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
Useful functions for EEG analysis.
"""
# ==============================================================================
#                             STANDARD LIBRARY IMPORTS
# ==============================================================================
import pickle
from copy import copy

# ==============================================================================
#                         IMPORTS CUSTOM MODULES AND PACKAGES
# ==============================================================================
import mne
import numpy as np

from mne.io import read_raw_fif
from mne.preprocessing import ICA, read_ica

# ==============================================================================
#                               CLASS AND FUNCTIONS
# ==============================================================================
def apply_custom_baseline(data, time_window=None, method="mean"):
    """apply_custom_baseline.
    Apply a custom baseline to the data for a time-frequency analysis.

    Args:
        data (:obj:`mne.io.Raw` | :obj:`mne.Epochs`): An instance of mne.io.Raw or mne.Epochs.
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
    if time_window is None:
        bld = data.copy().get_data()
    else:
        bld = data.copy().crop(*time_window).get_data()
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

def prompt_message(
    message,
    all_cap=True,
    line_separator="=",
    separator_position="both",
    message_alignment="centered",
):
    """prompt_message.
    Print messages on the prompt with nice formating

    Args:
        message (str): The message to print
        all_cap (bool, optional): Wether to pring the message in all capital letters. 
            Defaults to True.
        line_separator (str, optional): Choose the character that is use to delimitate the messages. 
            Defaults to "=".
        separator_position (str, optional): The position of the separator can be either "both", "above", "below". 
            Defaults to "both".
        message_alignment (str, optional): How to align the text either "centered", "left", "right". 
            Defaults to "centered".
    """
    if all_cap:
        message = message.upper()

    if separator_position == "both":
        upper_separator, lower_separator = (
            line_separator * 80 + "\n",
            "\n" + line_separator * 80,
        )
    elif separator_position == "above":
        upper_separator, lower_separator = line_separator * 80 + "\n", ""
    elif separator_position == "below":
        upper_separator, lower_separator = "", line_separator * 80 + "\n"

    if message_alignment == "left":
        position = ""
    elif message_alignment == "centered":
        position = " " * (40 - int(np.floor(len(message) * (2**-1))))
    elif message_alignment == "right":
        position = " " * (80 - len(message) - 1)

    print("\n" + upper_separator + position + message + lower_separator)

def run_ica(raw,n_components = 0.999, vizualisation=False, apply=False, **kwargs):
    """run_ica.
    Special wrapper for Independant Component Analysis (ICA) from mne-python.

    Args:
        raw (:obj:`mne.io.Raw`): Instance of mne.io.Raw.
        n_components (int | float, optional): The number of components to compute. 
            If `n_components` is `float` computes the maximum number of 
            component that explain `n_components` the standard deviation. For more
            details see: https://mne.tools/stable/generated/mne.preprocessing.ICA.html#mne.preprocessing.ICA
            Defaults to 0.999.

    """

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

def read_data_and_ica(raw_filename, ica_filename=None):
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
    raw.interpolate_bads(reset_bads=True)

    if ica_filename is not None:
        ica = read_ica(ica_filename)
        ica.apply(raw)

    return raw

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

            # Onset time of the fixation cross compared to the target (t=0)

            if metadata is None:
                fixation_cross_onset = data.tmin
                CUE_ONSET_TIME = data.tmax

            elif metadata is not None:
                fixation_cross_onset = -metadata["PreTargetDuration"].iloc[i] / 1000
                CUE_ONSET_TIME = -0.5

            fixation_duration[i] = abs(CUE_ONSET_TIME - fixation_cross_onset)

            # Extract data from fixation onset to cue onset for each epoch
            chunk = data[i].copy().crop(tmin=fixation_cross_onset, tmax=CUE_ONSET_TIME)

            chunk_array = chunk.get_data()

            # Get the signal only within the desired time window
            chunk_list.append(chunk_array)

            # Calculate the descriptive statistics
            mean[i, :, :] = chunk_array.mean(axis=2, keepdims=True)
            std[i, :, :] = chunk_array.std(axis=2, keepdims=True)
            median[i, :, :] = np.median(chunk_array, axis=2, keepdims=True)
            idx_max[i] = np.unravel_index(chunk_array.argmax(), chunk_array.shape)[2]
            idx_min[i] = np.unravel_index(chunk_array.argmin(), chunk_array.shape)[2]
            maximum[i, :, :] = np.transpose(chunk_array[:, :, int(idx_max[i])])
            minimum[i, :, :] = np.transpose(chunk_array[:, :, int(idx_min[i])])

            # Calculate the slope across each channel

            for channel in range(np.shape(chunk)[0]):
                slope[i, channel] = np.polyfit(
                    chunk.times, np.squeeze(chunk_array[:, channel, :]), deg=1
                )[0]

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
