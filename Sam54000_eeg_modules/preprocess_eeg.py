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
import pickle
from copy import copy
from pathlib import Path

# ==============================================================================
#                       IMPORTS CUSTOM MODULES AND PACKAGES
# ==============================================================================
import mne                                    # pip install mne or conda install mne
import numpy as np                            # pip install numpy or conda install numpy
import pyprep                                 # pip install pyprep or conda install -c conda-forge pyprep
import mne_bids
import pandas as pd

from matplotlib import pyplot as plt
from pytz import timezone
from mne.io import read_raw_fif
from mne.preprocessing import ICA, read_ica
from .cli import prompt_message
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
        raw (:obj:`mne.io.Raw`): raw data
        desired_ref_name (str, optional): desired name of the reference. Defaults to 'Cz'.
    
    Returns:
        raw (:obj:`mne.io.Raw`): raw data with the reference renamed in place.

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
        raw = []
        return raw, status

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
            f'{stim_var}04':'no_cue/congruent/dwn',
            f'{stim_var}05':'no_cue/incongruent/dwn',
            f'{stim_var}07':'center/congruent/up',
            f'{stim_var}08':'center/incongruent/up',
            f'{stim_var}10':'center/congruent/dwn',
            f'{stim_var}11':'center/incongruent/dwn',
            f'{stim_var}19':'spatial/congruent/up',
            f'{stim_var}20':'spatial/incongruent/up',
            f'{stim_var}22':'spatial/congruent/dwn',
            f'{stim_var}23':'spatial/incongruent/dwn',
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
            evt:event_desc.get(evt) for evt in stim_chan_names
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
                "event": event_desc.get(stim_chan_names[i]),
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
            ch_index = chan_stim_dict.get(key)["ch_index"]
            chan_stim_dict.get(key)["id_value"] = val + 1
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
                        f"id: {chan_stim_dict.get(key).get('id_value')}",
                        4 * " ",
                        f"event:{chan_stim_dict.get(key).get('event')}",
                        4 * " ",
                        f"raw_val: {np.unique(np.flatnonzero(raw._data[chan_stim_dict.get(key).get('ch_index')]))}",
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
            chan_stim_dict.get(key)["id_value"]: chan_stim_dict.get(key)["event"]

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
        raw (:obj:`mne.io.raw`): An instance of mne.io.Raw.
        montage_name (str, optional): provide the name of the montage of the cap used during the experiment.
            values has to follow the mne naming standard. For a list of acceptable
            values run `mne.channels.get_builtin_montages(descriptions=False).
            To get a list of tuples (montage name, description), pass descriptions
            to True. Defaults to "GSN-HydroCel-129".

    Returns:
        raw (:obj:`mne.io.raw`): An instance of mne.io.Raw with annotations
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

def extract_resting_block(raw, return_blocks=False):
    """
    Extract resting blocks from raw data based on annotations
    
    Args:
        raw (:obj:`mne.io.raw`): An instance of mne.io.Raw.
    
    Returns:
        raw_concatenated (:obj:`mne.io.raw`): An instance mne.io.Raw resulting from a concatenation of extracted blocks
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
        
        if return_blocks:
            return resting_blocks
        else: 
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
        raw (:obj:`mne.io.Raw`): An instance of mne.io.Raw.

    Returns:
        raw_concatenated (:obj:`mne.io.Raw`): An instance mne.io.Raw resulting from a concatenation of extracted blocks
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
        'RT'
        ]

    ch_index = [stim_chan_names.index(ch) for ch in stim_chan_names if ch.strip() in desired_elec]
    new_stim_chan_data = stim_chan_data[ch_index,:]

    raw.filter(1, None, verbose="CRITICAL")

    # Merge all stimuli together to have the ant_block shape discretized (all stim)
    new_stim_chan_data[new_stim_chan_data > 0] = 1
    target_data = new_stim_chan_data[:-1,:]
    response_data = new_stim_chan_data[-1,:]
    ant_block_discrete= np.sum(target_data, axis=0)

    bin_block = np.zeros(len(raw))

    # Get when a target trigger is on
    index = np.where(ant_block_discrete == 1)[0]
    diff_samples = np.diff(index)
    sorted_diff_samples = np.sort(diff_samples, axis=0)
    threshold = sorted_diff_samples[-3]    # How many seconds before and after the trigger to consider as a same event

    # Create blocks
    for i in index:
        bin_block[i:i+threshold] = 1

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
        raw (:obj:`mne.io.Raw`): An instance of mne.io.Raw.
        montage_name (str, optional): provide the name of the montage of the cap used during the experiment.
            values has to follow the mne naming standard. For a list of acceptable
            values run `mne.channels.get_builtin_montages(descriptions=False).
            To get a list of tuples (montage name, description), pass descriptions to True. Default "GSN-HydroCel-129"

    Returns:
        prep (:obj:`pyprep.prep`): An instance of pyprep.prep
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

def apply_custom_baseline(data, base_data = None, time_window=None, method="mean"):
    """apply_custom_baseline.
    Apply a custom baseline to the data for a time-frequency analysis.

    Args:
        data (:obj:`mne.io.Raw` | :obj:`mne.Epochs`): An instance of mne.io.Raw or mne.Epochs.
        base_data (:obj:`mne.Epochs` | :obj:`mne.Evoked`): An instance of mne.Epochs or mne.Evoked that will be used as a baseline.
            Useful when data has to be referenced to a separated baseline condition such as resting state (which is separated from data)
        time_window (tuple, optional): The time window in seconds to use for the baseline. Defaults to None.
        method (str, optional): Method for applying the baseline. Can be:
            - "mean": Data are substracted by the mean calculated within the time window.
            - "ratio": Data are divided by the mean calculated within the time window.
            - "zscore": Data are substracted by the mean calculated within the time window 
                and divided by the standard deviation.
            - "logratio": log10 of data divided by the mean calculated within the time window

    Returns:
        baseline_applied (:obj:`mne.time_frequency.EpochsTFR`): The time frequency calculated 
            with the baseline applied.
    """
    
    d = data.get_data()

    if not base_data:
        base_data = data.copy()

    if time_window is None:
        bld = base_data.copy().get_data()
    else:
        bld = base_data.copy().crop(*time_window).get_data()

    m = np.mean(bld, axis=-1, keepdims=True)

    if method == "mean":
        d -= m
    elif method == "ratio":
        d /= m
    elif method == "zscore":
        d -= m
        d /= np.std(data._data, axis=-1, keepdims=True)
    elif method == "logratio":
        d = np.log10(d / m)

    baseline_applied = mne.time_frequency.EpochsTFR(data.info,
                                 d,
                                 data.times,
                                 data.freqs,
                                 method="multitaper",
                                 events=data.events,
                                 event_id=data.event_id,
                                 metadata=data.metadata)

    return baseline_applied

def run_ica(raw,n_components = 0.999, vizualisation=False, apply=False, interpolate = False, **kwargs):
    """run_ica.
    Special wrapper for Independant Component Analysis (ICA) from mne-python.

    Args:
        raw (:obj:`mne.io.Raw`): Instance of mne.io.Raw.
        n_components (int | float, optional): The number of components to compute. 
            If `n_components` is `float` computes the maximum number of 
            component that explain `n_components` the standard deviation. For more
            details see: https://mne.tools/stable/generated/mne.preprocessing.ICA.html#mne.preprocessing.ICA
            Defaults to 0.999.
        vizualisation (bool, optional): If True, plot the sources. Defaults to False.
        apply (bool, optional): If True, apply the ICA to the raw data. Defaults to False.
        interpolate (bool, optional): If True, interpolate bad channels before running ICA. Defaults to False.

    """

    if interpolate:
        raw.interpolate_bads(reset_bads=True)

    ica = ICA(
        n_components=n_components,
        max_iter="auto",
        method="picard",
        random_state=97,
        verbose="CRITICAL",
    )
    ica.fit(raw)

    if vizualisation:
        ica.plot_sources(raw, block=True)

    if apply:
        ica.apply(raw)

        return ica, raw

    return ica

def read_eeg_and_ica(raw_filename, ica_filename=None, interpolate=False):
    """read_data_and_ica.
    Read the raw data and apply the ICA if an ICA object is provided.
    
    Args: 
        raw_filename (str | path-like): The filename of the raw data should be in '.fif' format
            because ICA is usually applied on preprocessed data previously done with mne-python
            (bad segments annotated, prep pipeline done).
        ica_filename (str | path-like | None, optional): The filename of the ICA object.
            Defaults to None.
    
    Returns:
        raw (:obj:`mne.io.Raw`): A mne.io.Raw instance with ICA applied.
    
    """
    raw = read_raw_fif(raw_filename, preload=True)
    if interpolate:
        raw.interpolate_bads(reset_bads=True)

    if ica_filename is not None:
        ica = read_ica(ica_filename)
        ica.apply(raw)

    return raw, ica

class epochs_stats:
    """epochs_stats.

    This class give the descriptive statistics of an mne.Epochs object.
    The descriptive statistics are calculated on the signal amplitude along the
    time dimension within each epoch. The time along which the statistics are
    calculated is defined between the onset time of the fixation cross and the
    cue onset time. The onset time of the fixation cross varies from one trial
    to another. The cue onset time is always the same.

    Args:
        data (:obj:`mne.Epochs`): Instance of mne.Epochs. It needs to contain metadata read
            from eprime file. The metadata must contain the column 'PreTargetDuration'.
        metadata (:obj:`pandas.DataFrame`): Metadata read from eprime file.

    Attributes:
        average (:obj:`numpy.ndarray`): The average calculated
        std (:obj:`numpy.ndarray`): The standard deviation calculated
        slope (:obj:`numpy.ndarray`): The slope calculated
        maximum (:obj:`numpy.ndarray`): Maximum value within the time window
        minimum (:obj:`numpy.ndarray`): Minimum value within the time window
        median (:obj:`numpy.ndarray`): The median calculated within the time window
        time_max (:obj:`numpy.ndarray`): Time in seconds of the maximum value within the time window
        time_min (:obj:`numpy.ndarray`): Time in seconds of the minimum value within the time window
        fixation_duration (tuple): Time in seconds between the fixation cross and the cue
        event_names (tuple): the name of the event

    """

    def __init__(self, data, metadata=None):
        self.data = data

        # Create empty arrays (size (nb of epochs, nb of channels, 1)) using
        # list comprehension

        (mean, std, slope, maximum, minimum, median,) = [
            np.empty(
                (
                    np.shape(data.get_data())[0],
                    np.shape(data.get_data())[1],
                    1,
                )
            )

            for i in range(6)
        ]

        # Create other empty array with different size (nb of epochs, 1)
        fixation_duration, idx_max, idx_min = [
            np.empty((np.shape(data.get_data())[0])) for i in range(3)
        ]

        # Create empty lists
        event_names, chunk_list, nb_epochs = list(), list(), list()

        # Loop through individual epochs

        for i in range(len(data)):
            """A more pythonic way to iterate the loop would be to use enumerate.
            However if I use `enumerate` epoch will be a numpy array whereas if I do epochs[i] it returns
            individual epochs object. For the process it is more interesting to get the Epochs object because I can get
            all the informations."""

            # Get the event names
            event_names.append(list(data[i].event_id.keys())[0])

        # Count the number of epochs for each condition

        for condition in data.event_id.keys():
            nb_epochs.append(len(data[condition]))

        self.channels = data.info["ch_names"]
        self.nb_epochs = nb_epochs
        self.mean = mean
        self.std = std
        self.slope = slope
        self.maximum = maximum
        self.time_max= idx_max * data.info["sfreq"]
        self.minimum = minimum
        self.time_min = idx_min * data.info["sfreq"]
        self.median = median
        self.fixation_duration = fixation_duration
        self.event_names = event_names
        self.idx_epochs = [i for i in range(len(data))]

        if metadata is not None:
            self.reaction_time = metadata["ReactionTime"].values
        self.time_series = chunk_list

    def pick_channel(self, channel_name):
        # TODO try to make a method like mne "pick_channels" where a list of
        # channel name can be passed
        """
        Return a copy of the epochs_stats object with only the selected channel.

       Args: 
            channel_name (str): The name of the channel to select.

        Returns
            epochs_stats (:obj:`epochs_stats`): The modified instance with statistics for only
                the selected channels.
        """

        # Get the index of the selected channel
        channel_idx = self.channels.index(channel_name)

        # Create a copy of the current object
        new_stats = copy(self)

        # Replace the data with the selected channel data
        new_stats.average = self.average[:, channel_idx, :]
        new_stats.std = self.std[:, channel_idx, :]
        new_stats.slope = self.slope[:, channel_idx]
        new_stats.maximum = self.maximum[:, channel_idx, :]
        new_stats.minimum = self.minimum[:, channel_idx, :]
        new_stats.median = self.median[:, channel_idx, :]
        new_stats.time_series = [i[:, channel_idx, :] for i in self.time_series]

        # Return the new object

        return new_stats

    def save(self, filename):
        """save.
        Save the object as a pickle file.

        Args:
            filename (str | path-like): Filename of the saved file.
        """
        with open(filename, "wb") as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

    # TODO def plot()

def parse_extrema(overall_metadata, value_to_parse = 'reaction_time'):
    """parse_extrema.
    Calculate the first and last quartile of data.
    
    Args:
        overall_data (:obj:`pandas.DataFrame`): Dataframe from which the extrema will be calculated.
            Should be a dataframe with the columns 'subject', 'value', 'reaction_time' and 'condition'.
        value_to_parse (str, optional): The column name of the value to parse. Defaults to 'reaction_time'.
    
    """
    quartile = round(0.25*np.median(overall_metadata.groupby('subject').value_counts()))
    if value_to_parse == 'reaction_time':
        qual = ['slowest','fastest']
    else:
        qual = ['highest','lowest']
    
    for h,q in zip(['head','tail'],qual):
        overall_metadata.loc[getattr(
            overall_metadata.sort_values(
                value_to_parse,
                ascending = False).groupby(['condition',
                                            'subject']),h)(quartile).index,value_to_parse+'qual'] = q

def apply_ssd(ssd, epochs):
    """apply_ssd.
    Apply the ssd object to the epochs object.
    
    Args:
        ssd (:obj:`mne.decoding.SSD`): The ssd object.
        epochs (:obj:`mne.Epochs`): The epochs object.
    
    Returns:
        trans_epoch (:obj:`mne.Epochs`): The transformed epochs object.
    """
    epochs_data = epochs.get_data()
    filtered_data = ssd.apply(epochs_data)
    tmin = epochs.tmin

    trans_epoch = mne.EpochsArray(data=filtered_data, info=epochs.info, tmin=tmin,
                                 events=epochs.events, event_id=epochs.event_id)
    trans_epoch.metadata = epochs.metadata

    return trans_epoch

def create_ssd_object(instance, band="theta", saving_filename=None, save=False):
    """create_ssd_object.
    Create the spectro spatial object estimated within the frequency band
    of interest.

    Args:
        instance (:obj:`mne.raw` | :obj:`mne.epochs`): The instance from which the ssd object will be estimated.
        band (str): Frequency band of interest.
        saving_filename (str | None): If save = `True`, the filename under which the ssd object will be saved. Defaults to `None`.
        save (bool): Indiciate if the object will be saved under the name of saving_filename. Defaults to `False`.

    Returns:
        ssd (:obj:`mne.decoding.SSD`): The ssd object estimated from the instance.
    """
    info = instance.info
    data = instance.get_data()
    frequency_band = dict(
        delta=(1, 3),
        theta=(4, 8),
        alpha=(8, 12),
        beta=(12, 30),
        gamma=(30, 100),
    )
    # freq_noise = (frequency_band[band][0]-1,frequency_band[1]+1)

    filt_parameter_signal = {
        "l_freq": frequency_band[band][0],
        "h_freq": frequency_band[band][1],
        "l_trans_bandwidth": 1,
        "h_trans_bandwidth": 1,
    }
    filt_parameter_noise = {
        "l_freq": frequency_band[band][0] - 1,
        "h_freq": frequency_band[band][1] + 1,
        "l_trans_bandwidth": 1,
        "h_trans_bandwidth": 1,
    }

    ssd = mne.decoding.SSD(
        n_components=4,
        info=info,
        rank="full",
        filt_params_signal=filt_parameter_signal,
        filt_params_noise=filt_parameter_noise,
        reg="oas",
    )
    ssd.fit(data)

    if save:
        with open(saving_filename, "wb") as outp:
            pickle.dump(ssd, outp, pickle.HIGHEST_PROTOCOL)

    return ssd

def plot_bad(root = None, subject = 'HC001', session = '01', task = 'ant', group = 'HC', montage_name = 'GSN-HydroCel-129'):
    
    """plot_bad
    Plot bad electrodes after pyprep and visual inspection

    Args:
        root (Posixpath, optional): Path to look. Defaults to None.
        subject (str, optional): Subject ID. Defaults to 'HC001'.
        session (str, optional): Session ID. Defaults to '01'.
        task (str, optional): Task name. Defaults to 'ant'.
        group (str, optional): Group ID. Defaults to 'HC'.
        montage_name (str, optional): Montage. Defaults to 'GSN-HydroCel-129'.
    
    Returns:
        fig (:obj:`matplotlib.figure.Figure`): The figure object.
        ax (:obj:`matplotlib.axes.Axes`): The axes object.
    """
    BIDS_path = mne_bids.BIDSPath(root = root,
                                 subject = subject,
                                 suffix='eeg',
                                 datatype='eeg',
                                 description='reviewed',
                                 session=session,
                                 extension='.fif',
                                 task=task
                                 )
    raw = mne_bids.read_raw_bids(BIDS_path)
    csv_file = Path(BIDS_path.root,'metadata_reviewed',f'metadata_review_{group}.csv')
    df = pd.read_csv(csv_file)
    d = df[df['subject'] == subject]
    nb_bad = df[f'auto_bad_elec_nb_{task}'][0] + df[f'bad_after_review_nb_{task}'][0]
    percentage = round(nb_bad * 100 / raw.info['nchan'],2)
    fig, ax = plt.subplots()
    string = [f'Number: {nb_bad}',
          f'Percentage: {percentage} %']
    plt.text(-0.1, 0.12,'\n'.join(string), fontsize=11)
    raw.plot_sensors(axes = ax)
    return fig, ax

def create_epochs(raw, tmin = -0.8, tmax=0.8, beh_data = None, locking = "stim", **kwargs):
    """create_epochs.
    Create epochs from raw data based on annotations during an ANT task.

    Args:
        raw (:obj: mne.io.Raw): The raw data processed with the prep pipeline.
        tmin (float): time in seconds to start the epoch.
        tmax (float): time in seconds to end the epoch.
        beh_data (:obj: pandas.DataFrame): The behavioral data in a pandas dataframe.
        locking (str, optional): Wether to lock the vizualisation on the stim "stim" or response
            "resp". Defaults to "stim".

    Returns:
        None
    """

    conditions = {"stim": 
        [
        "no_cue/congruent/up",
        "no_cue/incongruent/up",
        "no_cue/congruent/dwn",
        "no_cue/incongruent/dwn",
        "center/congruent/up",
        "center/incongruent/up",
        "center/congruent/dwn",
        "center/incongruent/dwn",
        "spatial/congruent/up",
        "spatial/incongruent/up",
        "spatial/congruent/dwn",
        "spatial/incongruent/dwn",
    ],
        "resp": ["response"]
    }
    
    if locking not in ["stim", "resp"]:
        raise ValueError("locking must be either 'stim' or 'resp'")
    if locking == "resp":
        beh_data = None
    
    # Create events from annotations
    events_from_annot = mne.events_from_annotations(raw)
    
    # Get the event id for each condition
    evt_dict = {key: events_from_annot[1][key] for key in conditions.get(locking)}
    
    # Get the event indices
    unfiltered_evt = events_from_annot[0]
    mask = np.isin(unfiltered_evt[:, 2], list(evt_dict.values()))
    evt = unfiltered_evt[mask]
    
    # Create the epochs
    epochs = mne.Epochs(
        raw,
        events=evt,
        event_id=evt_dict,
        tmin=tmin,
        tmax=tmax,
        preload=True,
        **kwargs,
    )

    # Harmonize IDS for cross subject comparison
    # Add 13 to all event IDs in the epoch
    epochs.events[:, 2] += max(list(epochs.event_id.values()))
    
    # Create a dictionary that maps original event IDs to new IDs
    event_id_map = {}
    for key, value in epochs.event_id.items():
        event_id_map[key] = value + max(list(epochs.event_id.values()))
        
    # Update the event_id dictionary with the new IDs
    epochs.event_id = event_id_map
    
    # Replace the original event IDs in the epoch with the new IDs
    for new_id, condition in enumerate(conditions.get(locking)):
        epochs.events[:, 2][
            np.where(epochs.events[:, 2] == epochs.event_id[condition])
        ] = (new_id + 1)
        epochs.event_id[condition] = new_id + 1
        
    # Get dataframe from annotations
    annotations = epochs.annotations  # get the annotations object
    annotations_df = annotations.to_data_frame()  # export annotation into dataframe
    annotations_df = annotations_df[annotations_df["description"].isin(conditions.get(locking))]
    annotations_df = annotations_df.reset_index()
    DeltaTime = list()
    timestamps = annotations_df["onset"].to_list()  # extract the timestamps
    for timestamp in timestamps:
        DT = timestamp.to_pydatetime()
        DT = timezone("UTC").localize(DT)
        DT = (DT - raw.info["meas_date"]).total_seconds()
        DeltaTime.append(DT)
    annotations_df["DeltaTime"] = DeltaTime
    
    # get dropped epochs index
    droped_epoch = [index for index, string in enumerate(epochs.drop_log) if string]
    
    # in case we want to add the metadata from the behavioral file
    if beh_data is not None:
        df = beh_data.dataframe
        df["description"] = (
            df["CueType"] + "/" + df["FlankerType"] + "/" + df["TargetType"]
        )
        df = df.reset_index()
        
        # HARMONIZE THESE GODDAMN EVENT BECAUSE SOMETIME THERE ARE SOME MISSING IN THE RAW!!!!!
        if len(df) > len(annotations_df):
            a = list(df["description"])
            b = list(annotations_df["description"])
            dataframe_a = left = df
            dataframe_b = right = annotations_df

        elif len(df) < len(annotations_df):
            a = list(annotations_df["description"])
            b = list(df["description"])
            dataframe_a = right = annotations_df
            dataframe_b = left = df

        elif len(df) == len(annotations_df):
            dataframe_a = left = df
            dataframe_b = right = annotations_df

        if len(dataframe_a) != len(dataframe_b):
            idx = list()
            while len(a) > len(b):
                for i in range(len(a)):
                    if a[i] != b[i]:
                        idx.append(i)
                        a.pop(i)
                        break
                    
            dataframe_a.drop(index=idx, inplace=True)
            dataframe_a.reset_index(inplace=True)
            
        left.reset_index(inplace=True)
        right.reset_index(inplace=True)
        metadata = left.join(right, lsuffix="_beh")
        metadata = metadata.drop(index=droped_epoch, axis=0)
        metadata = metadata.drop(columns=["index_beh", "level_0"])
        metadata.reset_index(inplace=True)
    else:
        metadata = annotations_df.reset_index()
        metadata = metadata.drop(index=droped_epoch, axis=0)
    epochs.metadata = metadata
    return epochs
