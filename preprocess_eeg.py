#!/usr/bin/env -S  python  #
# -*- coding: utf-8 -*-
# Samuel Louviot, PhD 
#
# MIT License
#
# Copyright (c) 2023 Samuel Louviot
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Useful functions to preprocess eeg files recorded with EGI system during ant experiments.
These functions are dealing with subtelties of: the EGI system, the ant task, the lab naming conventions (and their variations).
It has been developped when exploring and analyzing data in adult_TBI_attention project.
"""
# ==============================================================================
#                             STANDARD LIBRARY IMPORTS
# ==============================================================================
import re
import os

# ==============================================================================
#                       IMPORTS CUSTOM MODULES AND PACKAGES
# ==============================================================================
import mne                                    # pip install mne or conda install mne
from usefull_functions import prompt_message  # Think about coding setup.py
import numpy as np                            # pip install numpy or conda install numpy
import pyprep                                 # pip install pyprep or conda install -c conda-forge pyprep

# ==============================================================================
#                               CLASS AND FUNCTIONS
# ==============================================================================
def ref_naming(raw,desired_ref_name='Cz'):
    """ref_naming.
    For files recorded with EGI system the name of the reference is not the
    same. Sometime it's VertexReference, sometime it's VREF and sometime there
    is no reference at all. This little function deal with any reference name
    and rename it as desired. If it doesn't exist, create
    a flat channel named with the desired ref name (because ref channel in EGI system are flat)

    Args:
        raw (:obj:mne.io.Raw): raw data
        desired_ref_name (str, optional): desired name of the reference. Defaults to 'Cz'.
    
    Returns:
        raw (:obj:mne.io.Raw object): raw data with the reference renamed in place.

    """
    # Get the name of the reference
    ref_chan_name = [
        ch_name

        for ch_name in raw.info["ch_names"]

        if re.search("ref", ch_name, re.IGNORECASE)
    ]

    # Rename the reference to "Cz"

    if ref_chan_name:
        raw.rename_channels({ref_chan_name[0]:desired_ref_name})

    # If the reference is not there, create a flat channel
    elif not ref_chan_name and desired_ref_name not in raw.info["ch_names"]:
        flat_data = np.zeros((1, raw.n_times))
        # create an info object for the flat channel
        flat_info = mne.create_info(
            ch_names=[desired_ref_name], sfreq=raw.info["sfreq"], ch_types=["eeg"]
        )
        # create a rawArray object for the flat channel
        flat_chan = mne.io.RawArray(data=flat_data, info=flat_info)
        # add the flat channel to the raw object
        raw.add_channels([flat_chan])

def read_raw_eeg(filename, preload=False):
    """read_raw_eeg.
    Read eeg data from different format and standardize channel names by removing potential trailing spaces.
    Format allowed are:
    - egi (.mff, .RAW)
    - bdf (.bdf)
    - edf (.edf)
    - fif (.fif)
    - eeglab (.set)
    
    Args:
        filename (str): path to the file
        preload (bool, optional): If True, the data will be preloaded into memory. Defaults to False.
    
    Returns:
        raw: mne.Raw object
        status: (str) 'ok' if the file is read correctly, 'corrupted' if the file is corrupted
    """
    if os.path.exists(filename):
        extension = os.path.splitext(filename)[1]
        if extension == '.mff' or extension == '.RAW':
            reader = getattr(mne.io, 'read_raw_egi')
        elif extension == '.bdf':
            reader = getattr(mne.io, 'read_raw_bdf')
        elif extension == '.edf':
            reader = getattr(mne.io, 'read_raw_edf')
        elif extension == '.fif':
            reader = getattr(mne.io, 'read_raw_fif')
        elif extension == '.set':
            reader = getattr(mne.io, 'read_raw_eeglab')
        try:
            raw = reader(filename, preload=preload, exclude=None)
            # Remove spaces in channel names
            del_space = [i.strip() for i in raw.info["ch_names"]]
            chan_mapping = dict()

            for old_chan, new_chan in zip(raw.info["ch_names"], del_space):
                chan_mapping[old_chan] = new_chan
            raw.rename_channels(chan_mapping)
            ref_naming(raw)
            status = "ok"
        except:
            print(f"File {filename} is corrupted")
            status = "corrupted"
            raw = []

        return raw, status
    else:
        print(f"File {filename} doesn't exist")
        status = "nonexistent"
        return status

def get_stim_chan_name(raw, selection = None):
    """get_stim_chan_name.
    Get in alphabetical order triggers and their values from the stim channels.

    Args:
        raw (mne.Raw object): raw data
        selection (list, optional): list of channels to select. Defaults to None.
            If channels in selection are not stim channels, they will be ignored.
    
    Returns:
        stim_chan_names (list): list of stim channels names
        stim_chan_data (np.array): array of stim channels data
    """
    # Get stimulation channel names
    if selection:
        selection.sort()
        temp = raw.copy().pick('stim').info["ch_names"]
        stim_chan_names = [ch for ch in temp if ch in selection]
    else:
        stim_chan_names = raw.copy().pick('stim').info["ch_names"]
        stim_chan_names.sort()
        
    stim_chan_data = raw.copy().pick(stim_chan_names).get_data()
    

    return stim_chan_names, stim_chan_data

def run_eog_projection(
    raw,
    eog_electrodes=["E8", "E25", "E21", "E14", "E22", "E9"],
):
    """run_eog_projection.

    Parameters
    ----------
    raw : (`mne.Raw` or `mne.Epochs` objects)
        data object obtained with mne
    eog_electrodes : (`list` of `str`)
        list of eog electrodes to use for the eog projection
    """
    existing_eog = [
        i for i in eog_electrodes if i in raw.info["ch_names"]
    ]  # in case of some electrodes have been dropped

    eog_projs, eog_events = mne.preprocessing.compute_proj_eog(
        raw, ch_name=existing_eog
    )
    raw.add_proj(eog_projs).apply_proj()

    return raw

def preprocess_eeg_ant(raw, montage_name="GSN-HydroCel-129"):
    """preprocess_eeg_ant.
    Run the pipeline to prepare the raw eeg data obtained during ant task and
    recorded with an egi system.

    Args:
        raw (mne.Raw object): raw data
        montage_name : (str) Default "GSN-HydroCel-129"
            provide the name of the montage of the cap used during the experiment.
            values has to follow the mne naming standard. For a list of acceptable
            values run `mne.channels.get_builtin_montages(descriptions=False).
            To get a list of tuples (montage name, description), pass descriptions
            to True.
    
    Returns:
        raw (mne.Raw object): raw data with annotations
        evt_count (dict): dictionnary with the number of events for each event type
        
    """
    stim_chan_names, _ = get_stim_chan_name(raw)

    if raw.preload:
        # Event_description ordered to match with channels order
        if any(['Co' in name for name in stim_chan_names]):
            stim_var = 'Co'
        else:
            stim_var = '00'
        event_desc = {
            'SESS':'session_start',
            'CELL':'CELL',
            'bgin':'trial_start',
            'PFxS':'fixation_practice',
            'PWar':'cue_practice',
            'PRAC':'practice_start',
            'TRSP':'trial_end',
            'FxS+':'fixation',
            'War+':'cue',
            f'{stim_var}01':'no_cue/congruent/up',
            f'{stim_var}02':'no_cue/incongruent/up',
            f'{stim_var}04':'no_cue/congruent/down',
            f'{stim_var}05':'no_cue/incongruent/down',
            f'{stim_var}07':'center/congruent/up',
            f'{stim_var}08':'center/incongruent/up',
            f'{stim_var}10':'center/congruent/down',
            f'{stim_var}11':'center/incongruent/down',
            f'{stim_var}19':'spatial/congruent/up',
            f'{stim_var}20':'spatial/incongruent/up',
            f'{stim_var}22':'spatial/congruent/down',
            f'{stim_var}23':'spatial/incongruent/down',
            'ESTR':'session_start',
            'EEND':'session_end',
            'resp':'response',
            'RT':'response',
            'TSTR':'trial_start',
            'TEND':'trial_end',
            'eyec':'eyes_closed',
            'eyeo':'eyes_open'
        } 

        # Because stim values are not the same across file I needed to find a
        # walkaround by selecting only the desired channel and pick their value and then
        # remap the value to name the events.

        if 'STI 014' in stim_chan_names:
            stim_chan_names.pop(stim_chan_names.index('STI 014'))
            
        # Arange StimChan to get always the right value for the right event
        ordered_stim_chan_dict = {
            evt:event_desc[evt] for evt in stim_chan_names
        }
        # Find events based on stim channels
        # The old EGI system have only 0 or 1 on each channel so it can't be sorted
        # using stim values. Events values have to be assigned by hand

        # Delete spaces in event ID because there are blank spaces randomly!!
        raw.event_id = {
            key.strip(): value for key, value in raw.event_id.items()
        }
        
        idx = [raw.info["ch_names"].index(i) for i in ordered_stim_chan_dict.keys()]
        # Dictionnary with stim channels paramters
        chan_stim_dict = {
            stim_chan_names[i]: {
                "event": event_desc[stim_chan_names[i]],
                "ch_index": idx[i],
                "id_value": raw.event_id.get(stim_chan_names[i]),
            }

            for i in range(len(stim_chan_names))
        }
        raw.event_id = {
            key: value + max(raw.event_id.values()) + 1

            for key, value in raw.event_id.items()
        }

        for val, key in enumerate(chan_stim_dict.keys()):
            raw.event_id[key] = val + 1
            ch_index = chan_stim_dict[key]["ch_index"]
            chan_stim_dict[key]["id_value"] = val + 1
            raw._data[ch_index, np.nonzero(raw._data[ch_index, :])] = val + 1

        # Sanity check
        prompt_message(
            "- stim id reassignement -",
            line_separator=" ",
            separator_position="both",
            message_alignment="left",
        )

        for key in chan_stim_dict.keys():
            print(
                "".join(
                    [
                        f"channel name: {key}",
                        4 * " ",
                        f"id: {chan_stim_dict[key]['id_value']}",
                        4 * " ",
                        f"event:{chan_stim_dict[key]['event']}",
                        4 * " ",
                        f"raw_val: {int(np.unique(raw._data[chan_stim_dict[key]['ch_index']])[1])}",
                    ]
                )
            )

        events = mne.find_events(
            raw,
            output="onset",
            shortest_event=0,
            consecutive=False,
            stim_channel=stim_chan_names,
            verbose=False,
        )

        prompt_message(
            "- check id and events before anotating -",
            line_separator=" ",
            separator_position="both",
            message_alignment="left",
        )

        event_dict = {
            chan_stim_dict[key]["id_value"]: chan_stim_dict[key]["event"]

            for key in chan_stim_dict.keys()
        }

        # Generate annotation from events
        prompt_message(
            "- annotating -",
            line_separator=" ",
            separator_position="both",
            message_alignment="left",
        )
        annot_from_event = mne.annotations_from_events(
            events,
            event_desc=event_dict,
            sfreq=raw.info["sfreq"],
            orig_time=raw.info["meas_date"],
        )

        # Apply annotations to raw
        raw.set_annotations(raw.annotations + annot_from_event)

        montage = mne.channels.make_standard_montage(montage_name)
        raw.set_montage(
            montage, match_case=True, match_alias=False, on_missing="ignore"
        )

        annot = raw.annotations.description
        evt_count = {
            evt_name:np.where(annot == evt_name)[0].shape[0]
            for evt_name in event_desc.values()
        }
        return raw, evt_count

    else:
        prompt_message(
            """ Data are not preloaded into memory, please read the data by
                using the option preload = True""",
            all_cap=False,
            line_separator="*",
        )

def calculate_rt(raw):
    """calculate_rt.
    Calculate the response time for each trial.

    Args:
        raw (mne.Raw object): raw data

    Returns:
        rt (list): list of response time for each trial
    """
    
    # Event_description ordered to match with channels order
    if any(['Co' in name for name in raw.info["ch_names"]]):
        stim_var = 'Co'
        response_name = 'RT'
    else:
        stim_var = '00'
        response_name = 'resp'
        
    desired = [
        f'{stim_var}01',
        f'{stim_var}02',
        f'{stim_var}04',
        f'{stim_var}05',
        f'{stim_var}07',
        f'{stim_var}08',
        f'{stim_var}10',
        f'{stim_var}11',
        f'{stim_var}19',
        f'{stim_var}20',
        f'{stim_var}22',
        f'{stim_var}23',
        response_name
    ]
    selection = [ch for ch in raw.info["ch_names"] if ch in desired]
    stim_chan_names, stim_chan_data = get_stim_chan_name(raw, selection = selection)
    response_chan_data = stim_chan_data[stim_chan_names.index(response_name),:]
    stim_chan_data = np.delete(stim_chan_data, stim_chan_names.index(response_name), axis=0)
    response_chan_data[response_chan_data > 0] = 1
    stim_chan_data[stim_chan_data > 0] = 1
    stim_chan_data = np.squeeze(np.sum(stim_chan_data, axis=0))
    index = np.nonzero(stim_chan_data)[0]
    rt = list()
    cooldown = int(round(0.200*raw.info["sfreq"]))
    limit = int(round(1.7*raw.info["sfreq"]))
    for i in index:
        if any(response_chan_data[i+cooldown:i+limit].astype(bool)):
            rt.append(np.nonzero(response_chan_data[i+cooldown:i+limit,])[0][0]/raw.info["sfreq"])
        else:
            rt.append('N/A')
    return rt

def preprocess_eeg_resting(raw, montage_name="GSN-HydroCel-129"):
    """Detect eyes open and eyes close events and annotate them.

    Args:
        raw (:obj: mne.io.raw): An instance of mne.io.Raw.
        montage_name (str, optional): provide the name of the montage of the cap used during the experiment.
            values has to follow the mne naming standard. For a list of acceptable
            values run `mne.channels.get_builtin_montages(descriptions=False).
            To get a list of tuples (montage name, description), pass descriptions
            to True. Defaults to "GSN-HydroCel-129".

    Returns:
        raw (:obj: mne.io.raw): An instance of mne.io.Raw with annotations
        spec (dict): dictionnary with the number of events for each event type
    """
    EVENT_DURATION = 90
    prompt_message(
        "- looking for eyes open/close events -",
        line_separator=" ",
        separator_position="both",
        message_alignment="left",
    )
    montage_name = "GSN-HydroCel-129"
    sampling_frequency = raw.info["sfreq"]
    evts_name = ["eyeo", "eyec"]
    sorted_evt = dict(eyeo=np.array([]), eyec=np.array([]))

    # Check the existence of channels (membership of 2 lists)
    check_membership = [item in raw.info["ch_names"] for item in evts_name]

    eyes_open_exists = check_membership[0]  # "eyeo" in raw.info["ch_names"]
    eyes_close_exists = check_membership[1]  # "eyec" in raw.info["ch_names"]

    time_epoch = dict()
    if all(check_membership):
        prompt_message(
            "Eyes close and eyes open events found",
            line_separator="",
            separator_position="above",
            message_alignment="left",
        )
        # Check existence of events in the channels

        for evt in evts_name:
            output_evt = mne.find_events(raw, stim_channel=[evt], verbose="CRITICAL")
            sorted_evt[evt] = output_evt
            time_epoch[evt] = np.diff(output_evt[:, 0],n=1) / sampling_frequency
            print(f"{evt} time: {time_epoch[evt]}")

        # FORCE EVT ID TO 1 FOR EYES OPEN AND 2 FOR EYES CLOSE BECAUSE IT IS
        # NOT THE SAME VALUE ACCROSS FILES!!!!!!

        if sorted_evt["eyeo"].any():
            sorted_evt["eyeo"][:,2] = 1

        if sorted_evt["eyec"].any():
            sorted_evt["eyec"][:,2] = 2

        resting_events = np.array([sorted_evt["eyeo"][0],
                                                 sorted_evt["eyec"][0]])
        MARGIN = 5
        onset_seconds = resting_events[:, 0] / sampling_frequency

        prompt_message(
            "- annotating -",
            line_separator=" ",
            separator_position="both",
            message_alignment="left",
        )
        events_annotations = mne.Annotations(
            onset=onset_seconds + MARGIN,
            duration=EVENT_DURATION,
            description=["eyes_open", "eyes_close"],
            orig_time=raw.info["meas_date"],
        )

        raw.set_annotations(raw.annotations + events_annotations)
        montage = mne.channels.make_standard_montage(montage_name)
        # Apply montage to raw object
        raw.set_montage(
            montage, match_case=True, match_alias=False, on_missing="ignore"
        )

    if not eyes_open_exists:
        prompt_message(
            "WARNING: eyes open event not found. Can't define an exact eyes open event.",
            line_separator="*",
            separator_position="above",
            message_alignment="left",
        )

    if not eyes_close_exists:
        prompt_message(
            "WARNING: eyes close event not found. Can't define an exact eyes close event.",
            line_separator="*",
            separator_position="above",
            message_alignment="left",
        )
    spec = dict(
        trigger_eyeo_nb = len(sorted_evt["eyeo"]),
        trigger_eyec_nb = len(sorted_evt["eyec"]),
        time_epoch_eyeo=str(time_epoch.get("eyeo")),
        time_epoch_eyec=str(time_epoch.get("eyec")),
        event_duration_eyeo = EVENT_DURATION,
        event_duration_eyec = EVENT_DURATION,
    )
    return raw, spec

def extract_resting_block(raw):
    """
    Extract resting blocks from raw data based on annotations
    
    Args:
        raw (:obj: mne.io.raw): An instance of mne.io.Raw.
    
    Returns:
        raw_concatenated (:obj: mne.io.raw): An instance mne.io.Raw resulting from a concatenation of extracted blocks
    """

    prompt_message(
        "- extracting resting blocks -",
        line_separator=" ",
        separator_position="both",
        message_alignment="left",
    )

    if (
        "eyes_open" in raw.annotations.description
        and "eyes_close" in raw.annotations.description
    ):
        annot = raw.annotations
        resting_evt_idx = list(
            np.where(
                np.logical_or(
                    annot.description == "eyes_close", annot.description == "eyes_open"
                )
            )[0]
        )
        resting_blocks = list()
        margin_time = 0.3

        for idx in resting_evt_idx:
            tmin_crop = annot[idx]["onset"] - margin_time
            tmax_crop = annot[idx]["onset"] + annot[idx]["duration"] + margin_time
            resting_blocks.append(raw.copy().crop(tmin=tmin_crop, tmax=tmax_crop))
        raw_concatenated = mne.concatenate_raws(resting_blocks)

        return raw_concatenated
    else:
        prompt_message(
            """One or all annotations of resting states are missing,
            please check if 'eyes_open' and/or 'eyes_close' exist""",
            all_cap=False,
            line_separator="*",
            separator_position="above",
            message_alignment="left",
        )

def extract_ant_blocks(raw):
    """extract_ant_blocks.
    Extract ant blocks from raw data based on annotations

    Args:
        raw (:obj: mne.io.Raw): An instance of mne.io.Raw.

    Returns:
        raw_concatenated (:obj: mne.io.Raw): An instance mne.io.Raw resulting from a concatenation of extracted blocks
        pulses_stats (dict): dictionnary with the number of events for each event type
    """
    prompt_message(
        "- extracting blocks -",
        line_separator=" ",
        separator_position="both",
        message_alignment="left",
    )
    stim_chan_names, stim_chan_data = get_stim_chan_name(raw)

    if any(['Co' in name for name in stim_chan_names]):
        stim_var = 'Co'
    else:
        stim_var = '00'

    desired_elec = [
        'bgin',
        'FxS+',
        'War+',
        f'{stim_var}01',
        f'{stim_var}02',
        f'{stim_var}04',
        f'{stim_var}05',
        f'{stim_var}07',
        f'{stim_var}08',
        f'{stim_var}10',
        f'{stim_var}11',
        f'{stim_var}19',
        f'{stim_var}20',
        f'{stim_var}22',
        f'{stim_var}23',
        'resp',
        'RT',
        'TSTR',
        'TEND']
    
    ch_index = [stim_chan_names.index(ch) for ch in stim_chan_names if ch in desired_elec]
    new_stim_chan_data = stim_chan_data[ch_index,:]
    
    prompt_message(
        "filtering",
        line_separator="",
        separator_position="below",
        message_alignment="left",
    )
    raw.filter(1, None, verbose="CRITICAL")

    # Merge all stimuli together to have the ant_block shape discretized (all stim)
    new_stim_chan_data[new_stim_chan_data > 0] = 1
    target_data = new_stim_chan_data[:-1,:]
    response_data = new_stim_chan_data[-1,:]
    ant_block_discrete= np.sum(target_data, axis=0)

    bin_block = np.zeros(len(raw))

    # Get when a target trigger is on
    idx = np.squeeze(np.nonzero(ant_block_discrete))

    # How many seconds before and after the trigger to consider as a same event
    threshold = 3 * int(raw.info["sfreq"])

    # Create blocks
    for i in idx:
        bin_block[i-threshold:i+threshold] = 1

    # Detect edges
    edges = np.diff(bin_block)
    rising = np.where(edges == 1)[0]
    falling = np.where(edges == -1)[0]

    # Deal when data start or finish without block returning to 0 (to close to 
    # the edge of the data)
    if falling[-1] < rising[-1]:
        # If the last block is not finished, go back to 4 samples to finish the
        # block
        falling = np.append(falling, len(bin_block) - 4)
    if rising[0] > falling[0]:
        rising = np.insert(rising, 0, 0)

    idx_block_duration = [
        stop - start for start, stop in zip(rising,falling)
    ]

    onsets = [raw.times[start] for start in rising]
    duration = [raw.times[duration_idx] for duration_idx in idx_block_duration]

    # Count number of targets and responses

    if len(onsets) > 3:
        onsets.pop(0)
        duration.pop(0)
        rising = np.delete(rising, 0)
        falling = np.delete(falling, 0)

    nb_targets, nb_responses = [], []

    for r,f in zip(rising, falling):
        nb_targets.append(sum(ant_block_discrete[r:f]))
        nb_responses.append(sum(response_data[r:f]))

    np_total_targets = np.sum(nb_targets)
    nb_total_responses = np.sum(nb_responses)

    pulses_stats = {
        "nb_blocks": len(onsets),
        "triggers_target_nb": str(nb_targets),
        "triggers_response_nb": str(nb_responses),
        "triggers_target_nb_total": np_total_targets,
        "triggers_response_nb_total": nb_total_responses
    }

    block_annotations = mne.Annotations(
        onsets,
        duration,
        [f"block_{n}" for n in range(1, len(onsets) + 1)],
        orig_time=raw.info["meas_date"],
    )
    raw.set_annotations(raw.annotations + block_annotations)

    for i in range(len(rising)):
        print(f"Block{i+1}: {nb_targets[i]} trials | {nb_responses[i]} responses")
    blocks = raw.crop_by_annotations(annotations=block_annotations)
    raw_concatenated = mne.concatenate_raws(blocks)

    return raw_concatenated, pulses_stats

def run_prep(raw, montage_name="GSN-HydroCel-129"):
    """run_prep.
    Run the prep pipeline on the raw data

    Args:
        raw (:obj: mne.io.Raw): An instance of mne.io.Raw.
        montage_name (str, optional): provide the name of the montage of the cap used during the experiment.
            values has to follow the mne naming standard. For a list of acceptable
            values run `mne.channels.get_builtin_montages(descriptions=False).
            To get a list of tuples (montage name, description), pass descriptions to True. Default "GSN-HydroCel-129"

    Returns:
        prep (:obj: pyprep.prep): An instance of pyprep.prep
    """
    prompt_message(
        "- running prep -",
        line_separator=" ",
        separator_position="both",
        message_alignment="left",
    )
    montage = mne.channels.make_standard_montage(
        montage_name
    )  # to modify to be adaptable to any egi nets

    if raw.info["sfreq"] > 250:
        prompt_message(
            "Downsampling to 250 Hz",
            all_cap=False,
            line_separator="",
            message_alignment="left",
        )

        raw.filter(l_freq=1, h_freq=125, verbose="CRITICAL")
        raw.resample(250, verbose="CRITICAL")

    if not raw.preload:
        raw.load_data()

    prep_params = {
        "ref_chs": "eeg",
        "reref_chs": "eeg",
        "line_freqs": np.arange(60, raw.info["sfreq"] / 2, 60),
    }
    prompt_message(
        "Executing PrepPipeline",
        line_separator="",
        message_alignment="left",
        separator_position="both",
    )
    prep = pyprep.PrepPipeline(raw, prep_params, montage, channel_wise=True)
    prep.fit()

    return prep
