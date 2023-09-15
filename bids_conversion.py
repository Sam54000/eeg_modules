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
This code is used to convert informations in GeoScan files (in an xml format) into BIDS readable files.
It generates a tsv file containing the electrodes coordinates and a json file containing the coordinates system.

An EEG folder in a BIDS format should include the following files:
    - A sidecar JSON file that contains general information about the EEG recording.
        It has to be named according to the following convention:
        sub-<participant_label>[_ses-<session_label>]_task-<task_label>[_acq-<label>][_rec-<label>][_run-<index>]_eeg.json
    - A table file that contains the channel information. A channel is described as any stream of data which can be TTL triggers,
        respiratory belt, EEG electrodes, etc. It has to be named according to the following convention:
        sub-<participant_label>[_ses-<session_label>]_task-<task_label>[_acq-<label>][_rec-<label>][_run-<index>]_channels.tsv
    - A channels description file that give additional information about the entry of the channels table file.
        It has to be named according to the following convention:
        sub-<participant_label>[_ses-<session_label>]_task-<task_label>[_acq-<label>][_rec-<label>][_run-<index>]_channels.json
    - An electrodes description file table that contains the EEG electrodes coordinates.
        It has to be named according to the following convention:
        sub-<participant_label>[_ses-<session_label>]_task-<task_label>[_acq-<label>][_rec-<label>][_run-<index>]_electrodes.tsv
    - An electrodes coordinates system description file that contains the EEG electrodes coordinates system.
        It has to be named according to the following convention:
        sub-<participant_label>[_ses-<session_label>]_task-<task_label>[_acq-<label>][_rec-<label>][_run-<index>]_coordsystem.json
    - An events file that contains in a table the events information. It has to be named according to the following convention:
        sub-<participant_label>[_ses-<session_label>]_task-<task_label>[_acq-<label>][_rec-<label>][_run-<index>]_events.tsv
    - An events description file that give additional information about the entry of the events table file.
        It has to be named according to the following convention:
        sub-<participant_label>[_ses-<session_label>]_task-<task_label>[_acq-<label>][_rec-<label>][_run-<index>]_events.json
    - An EEG file that contains the EEG data. It should be in a standardized format such as .edf 
        (see: https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html#eeg-recording-data). 
        It has to be named according to the following convention:
        sub-<participant_label>[_ses-<session_label>]_task-<task_label>[_acq-<label>][_rec-<label>][_run-<index>]_eeg.edf
"""
# ==============================================================================
#                       IMPORTS STANDARD MODULES AND PACKAGES
# ==============================================================================
import json
import os
import argparse
import sys
import xml.etree.ElementTree as ET

from pathlib import Path

# ==============================================================================
#                       IMPORTS CUSTOM MODULES AND PACKAGES
# ==============================================================================
import pandas as pd #pip install pandas or conda install pandas
import mne #pip install mne or conda install mne
import mne_bids #pip install mne-bids or conda install mne-bids

sys.path.append('/Users/samuel/codes/modules') #Think about making a setup.py file
from preprocess_ant import *

# ==============================================================================
#                               CLASS AND FUNCTIONS
# ==============================================================================
class Convertor:
    """Convertor
    This class is used to convert a raw EGI eeg files into a set of BIDS compatible files.
    
    Args:
        BIDSpath (:obj: `mne_bids.BIDSPath`): Path object generated from mne_bids.BIDSPath constructor. It has to point to the saving EEG folder.
        eeg_filename (str, path-like): Path to the raw EEG file.
        electrodes_loc_xml (str, path-like): Path to the xml file containing the electrodes coordinate if the file exists.
        
    Attributes:
        BIDSpath (:obj: `mne_bids.BIDSPath`): Path object generated from mne_bids.BIDSPath constructor. It has to point to the saving EEG folder.
        eeg_filename (str, path-like): Path to the raw EEG file.
        electrodes_loc_xml (str, path-like): Path to the xml file containing the electrodes coordinate if the file exists.
        BIDSentities_from_eeg (dict): Dictionary containing the BIDS entities extracted from the eeg_filename.
        raw (:obj: `mne.io.Raw`): Raw object containing the EEG data.
    """
    def __init__(self,BIDSpath,electrodes_loc_xml=None,eeg_filename=None):
        self.electrodes_loc_xml=electrodes_loc_xml
        self.BIDSpath = BIDSpath
        fileparts = os.path.splitext(eeg_filename)
        self.BIDSentities_from_eeg= mne_bids.get_entities_from_fname(fileparts[0])
        self.BIDSpath.mkdir()
        self.raw,_= read_raw_eeg(eeg_filename, preload=True)

        if self.BIDSentities_from_eeg['task'].lower() == 'ant':
            self.raw,_ = preprocess_eeg_ant(self.raw, montage_name="GSN-HydroCel-129")
        elif self.BIDSentities_from_eeg['task'].lower() == 'rest':
            self.raw,_ = preprocess_eeg_resting(self.raw, montage_name="GSN-HydroCel-129")

    def electrodes_description(self,
                               elec_type = 'sponge',
                               ref_name = 'Cz'):
        """generate the electrodes description tsv and coordsystem json files
        
        Args:
            elec_type (str, optional): Type of the electrodes. Defaults to 'sponge'.
            ref_name (str, optional): Name of the reference electrode. Defaults to 'Cz'.

        Returns:
            self(:obj: Convertor): instance of the Convertor class
        """
        if self.electrodes_loc_xml is None:
            print('No electrodes_loc_xml file specified')
            return
        else:
            # Read and parse elements in the xml file
            tree = ET.parse(self.electrodes_loc_xml)
            root = tree.getroot()
            
            # Create a dataframe containing the electrodes coordinates
            electrodes= pd.DataFrame(columns=['name', 'x', 'y', 'z', 'type'])
            
            # Enter the coordinates in the dataframe for each electrodes found
            for elec in range(len(root[0][1])):
                if root[0][1][elec][2].text != '2':
                    if root[0][1][elec][2].text == '0':
                        electrodes.loc[elec, 'name'] = f'E{root[0][1][elec][1].text}'
                        electrodes.loc[elec, 'type'] = elec_type
                    elif root[0][1][elec][2].text == '1':
                        electrodes.loc[elec, 'name'] = ref_name
                        electrodes.loc[elec, 'type'] = elec_type
                    electrodes.loc[elec, 'x'] = root[0][1][elec][3].text
                    electrodes.loc[elec, 'y'] = root[0][1][elec][4].text
                    electrodes.loc[elec, 'z'] = root[0][1][elec][5].text
                elif root[0][1][elec][2].text == '2':
                    if "right" in root[0][1][elec][0].text.lower():
                        rpa = [float(i.text) for i in root[0][1][elec][3:6]]
                    elif "left" in root[0][1][elec][0].text.lower():
                        lpa = [float(i.text) for i in root[0][1][elec][3:6]]
                    elif "nasion" in root[0][1][elec][0].text.lower():
                        nas = [float(i.text) for i in root[0][1][elec][3:6]]

            self.BIDSpath.update(suffix = 'electrodes', extension = '.tsv')
            electrodes.to_csv(self.BIDSpath.fpath)
            print_string = f'electrodes description saved in:'
            spaces = ' ' * (39 - len(print_string))
            print(f'{print_string}{spaces}{self.BIDSpath.fpath}')

            self.BIDSpath.update(suffix = 'eeg', extension = '.edf')
            coordsystem = {'IntentedFor': str(self.BIDSpath.fpath),
            'EEGCoordinateSystem': 'CapTrak',
            'EEGCoordinateUnits': 'cm',
            'EEGCoordinateSystemDescription': 'RAS orientation: Origin halfway between LPA and RPA, positive x-axis towards RPA, positive y-axis orthogonal to x-axis through Nasion,  z-axis orthogonal to xy-plane, pointing in superior direction',
            'FiducialDescription': 'Electrodes and fiducials were digitized with a GeoScan system from EGI',
            'FiducialCoordinates': {'NAS': rpa, 
                                    'LPA': lpa,
                                    'RPA': nas},
            'FiducialsCoordinateSystem': 'CapTrak',
            'FiducialsCoordinateUnits': 'cm',
            'FiducialCoordinateSystemDescription': 'The three fiducials are NAS for nasion the most anterior point on the skull, LPA for left preauricular the most lateral point on the left ear, and RPA for right preauricular the most lateral point on the right ear. Each fiducial coordinates are expressed with a list of 3 values representing the coordinates in x,y and z, respectively. The coordinate system is CapTrak, which is a right-handed Cartesian coordinate system with the origin halfway between LPA and RPA, positive x-axis towards RPA, positive y-axis orthogonal to x-axis through Nasion,  z-axis orthogonal to xy-plane, pointing in superior direction'} 
            json_obj = json.dumps(coordsystem)
            self.BIDSpath.update(suffix = 'coordsystem', extension = '.json')
            with open(self.BIDSpath.fpath, "w") as outfile:
                outfile.write(json_obj)
            
            print_string = f'electrodes coordinate system saved in:'
            spaces = ' ' * (39 - len(print_string))
            print(f'{print_string}{spaces}{self.BIDSpath.fpath}')
            return self
    
    def channel_description(self):
        """generate the channel description tsv file
        
        Returns:
            self(:obj: Convertor): instance of the Convertor class
        """
        self.BIDSpath.update(**self.BIDSentities_from_eeg, suffix = 'channels', extension = '.tsv')
        
        channel_description = pd.DataFrame(columns=['name', 'type', 'units', 'description', 'sampling_frequency', 'low_cutoff', 'high_cutoff','status'])
        for i, ch in enumerate(self.raw.info['chs']):
            name = str(ch['ch_name'])
            temp = str(ch['unit'])
            unit = temp[temp.index('FIFF'):-1].split('_')[-1]
            temp = str(ch['kind'])
            ch_type = temp[temp.index('FIFF'):-1].split('_')[1]
            description = 'n/a'
            stim_vars = ['Co','00']
            for stim_var in stim_vars:
                special_channels = {
                    'RespThorax': {'ch_type':'RESP',
                                'description':'respiration'},
                    'RespAbdomen': {'ch_type':'RESP',
                                'description':'respiration'},
                    'ECG': {'ch_type':'ECG',
                                'description':'electrocardiogram'},
                    'VREF': {'ch_type':'REF',
                                'description':'reference placed on Cz'},
                    'SESS': {'ch_type':'TRIG',
                                'description':'start of the session'},
                    'CELL': {'ch_type':'TRIG',
                                'description':'unknown'},
                    'bgin': {'ch_type':'TRIG',
                                'description':'start of a trial'},
                    'resp': {'ch_type':'TRIG',
                                'description':'response from the participant when they pressed a button'},
                    'PFxS': {'ch_type':'TRIG',
                                'description':'Fixation cross onset during practice trials'},
                    'PWar': {'ch_type':'TRIG',
                                'description':'cue onset during practice trials'},
                    'PRAC': {'ch_type':'TRIG',
                                'description':'start of a practice trial'},
                    'TRSP': {'ch_type':'TRIG',
                                'description':'end of a trial'},
                    'FxS+': {'ch_type':'TRIG',
                                'description':'Fixation cross onset'},
                    'War+': {'ch_type':'TRIG',
                                'description':'cue onset'},
                    f'{stim_var}01': {'ch_type':'TRIG',
                                'description':'no_cue/congruent/up'},
                    f'{stim_var}02': {'ch_type':'TRIG',
                                'description':'no_cue/incongruent/up'},
                    f'{stim_var}04': {'ch_type':'TRIG',
                                'description':'no_cue/congruent/down'},
                    f'{stim_var}05': {'ch_type':'TRIG',
                                'description':'no_cue/incongruent/down'},
                    f'{stim_var}07': {'ch_type':'TRIG',
                                'description':'center/congruent/up'},
                    f'{stim_var}08': {'ch_type':'TRIG',
                                'description':'center/incongruent/up'},
                    f'{stim_var}10': {'ch_type':'TRIG',
                                'description':'center/congruent/down'},
                    f'{stim_var}11': {'ch_type':'TRIG',
                                'description':'center/incongruent/down'},
                    f'{stim_var}19': {'ch_type':'TRIG',
                                'description':'spatial/congruent/up'},
                    f'{stim_var}20': {'ch_type':'TRIG',
                                'description':'spatial/incongruent/up'},
                    f'{stim_var}22': {'ch_type':'TRIG',
                                'description':'spatial/congruent/down'},
                    f'{stim_var}23': {'ch_type':'TRIG',
                                'description':'spatial/incongruent/down'},
                    'EEND': {'ch_type':'TRIG',
                                'description':'end of the session'},
                    'ESTR': {'ch_type':'TRIG',
                                'description':'start of the session'},
                    'RT':   {'ch_type':'TRIG',
                                'description':'response from the participant when they pressed a button'},
                    'TEND': {'ch_type':'TRIG',
                                'description':'end of a trial'},
                    'TSTR': {'ch_type':'TRIG',
                                'description':'start of a trial'},
                    'eyec': {'ch_type':'TRIG',
                                'description':'eye closed'},
                    'eyeo': {'ch_type':'TRIG',
                                'description':'eye open'}
                }

                if name.strip() in special_channels.keys():
                    ch_type = special_channels[name.strip()]['ch_type']
                    description = special_channels[name.strip()]['description']
                
            sampling_frequency =self.raw.info['sfreq']
            low_cutoff =self.raw.info.get('highpass')
            high_cutoff =self.raw.info.get('lowpass')

            if name in self.raw.info['bads']:
                status = 'bad'
            elif 'stim' in name.lower() or 'resp' in name.lower() or 'ecg' in name.lower():
                status = 'n/a'
            else:
                status = 'good'
            channel_description.loc[i] = [str(name), ch_type, unit, description, sampling_frequency, low_cutoff, high_cutoff, status]
        channel_description.to_csv(self.BIDSpath.fpath)
        print_string = f'channels description saved in:'
        spaces = ' ' * (39 - len(print_string))
        print(f'{print_string}{spaces}{self.BIDSpath.fpath}')
        
        return self

    def generate_events(self):
        """generate the events tsv and json files

        Returns:
            self (:obj: Convertor): instance of the Convertor class
        """
        self.BIDSpath.update(**self.BIDSentities_from_eeg, suffix = 'events', extension = '.tsv')
        annot = self.raw.annotations
        events = pd.DataFrame(columns=['onset', 'duration', 'trial_type', 'response_time'])
        if 'rest' in self.BIDSentities_from_eeg['task'].lower():
            mapping = [
                'eyes_open',
                'eyes_close',
            ]
        else:
            rt = calculate_rt(self.raw)
            events['response_time'] = rt
            mapping = [
                'no_cue/congruent/up',
                'no_cue/incongruent/up',
                'no_cue/congruent/down',
                'no_cue/incongruent/down',
                'center/congruent/up',
                'center/incongruent/up',
                'center/congruent/down',
                'center/incongruent/down',
                'spatial/congruent/up',
                'spatial/incongruent/up',
                'spatial/congruent/down',
                'spatial/incongruent/down',
            ]

        onset = list()
        duration = list()
        trial_type = list()
        for i, desc in enumerate(annot.description):
            if desc in mapping:
                onset.append(annot.onset[i])
                duration.append(annot.duration[i])
                trial_type.append(desc)

        events['onset'] = onset
        events['duration'] = duration
        events['trial_type'] = trial_type
        
        events.to_csv(self.BIDSpath.fpath)
        print_string = f'events saved in:'
        spaces = ' ' * (39 - len(print_string))
        print(f'{print_string}{spaces}{self.BIDSpath.fpath}')
        
        
        
        description = {
            'onset': 'onset of the event in seconds from the beginning of the recording raw.info["meas_date"]',
            'duration': 'duration of the event in seconds',
            'trial_type':{
                "LongName": "Attention Network Task (Flanker Task)",
                "Description": '''A set of five arrows appears above or under a fixation cross located in the middle of the screen. 
                Participants were instructed to answer, as fast as possible, which direction the middle arrow is pointing to by pressing the appropriate response button. 
                The arrows above or under the fixation cross are either pointing in the same direction (congruent condition) 
                or in different directions (incongruent condition). 
                A cue appears 500ms before the stimulus either in the middle of the screen (center condition) or in place of where the stimulus will appear (spatial condition).
                ''',
                "Levels": {
                    'no_cue/congruent/up': 'no cue, stimulus with congruent flankers appearing above the fixation cross',
                    'no_cue/incongruent/up': 'no cue, stimulus with incongruent flankers appearing above the fixation cross',
                    'no_cue/congruent/down': 'no cue, stimulus with congruent flankers appearing under the fixation cross',
                    'no_cue/incongruent/down': 'no cue, stimulus with incongruent flankers appearing under the fixation cross',
                    'center/congruent/up': 'center cue, stimulus with congruent flankers appearing above the fixation cross',
                    'center/incongruent/up': 'center cue, stimulus with incongruent flankers appearing above the fixation cross',
                    'center/congruent/down': 'center cue, stimulus with congruent flankers appearing under the fixation cross',
                    'center/incongruent/down': 'center cue, stimulus with incongruent flankers appearing under the fixation cross',
                    'spatial/congruent/up': 'spatial cue, stimulus with congruent flankers appearing above the fixation cross',
                    'spatial/incongruent/up': 'spatial cue, stimulus with incongruent flankers appearing above the fixation cross',
                    'spatial/congruent/down': 'spatial cue, stimulus with congruent flankers appearing under the fixation cross',
                    'spatial/incongruent/down': 'spatial cue, stimulus with incongruent flankers appearing under the fixation cross',
                },
            "response_time": "Time in seconds between the onset of the stimulus and when the participant pressed the button",
            }
        }
        self.BIDSpath.update(extension='.json')
        json_obj = json.dumps(description, indent=None)
        with open(f'{self.BIDSpath.fpath}', "w") as outfile:
            outfile.write(json_obj) 
        print_string = f'events description saved in:'
        spaces = ' ' * (39 - len(print_string))
        print(f'{print_string}{spaces}{self.BIDSpath.fpath}')
        
        return self
    
    def generate_sidecar_json(self,
                              task_description = None,
                              instructions = None,
                              cogatlas = None,
                              manufacturer = 'EGI',
                              device = 'NET Amp 400',
                              device_serial_number = None,
                              software_versions = None,
                              powerline_frequency = 60,
                              eeg_placement_scheme = '5-10',
                              ref_name = 'Cz',
                              eeg_ground = 'placed on CPz',
                              recording_type = 'continuous'
                              ):
        """generate the sidecar json file

        Args:
            task_description (str, optional): Longer description of the task. Defaults to None.
            instructions (str, optional): Instructions that has been given to the participant. 
                It can be the prompt before the task. Defaults to None.
            cogatlas ('str', optional): URL of the CogAtlas for the task. Defaults to None.
            manufacturer (str, optional): Manufacturer of the EEG amplifier. Defaults to 'EGI'.
            device (str, optional): Model name of the amplifier. Defaults to 'NET Amp 400'.
            device_serial_number (str, optional): Serial number of the amplifier. Defaults to None.
            software_versions (str, optional): Version of the recording software used. Defaults to None.
            powerline_frequency (int, optional): Powerline frequency. Defaults to 60 (USA). Put 50 if data
                was recorded in Europe.
            eeg_placement_scheme (str, optional): EEG system must be one of:
                '10-20','10-10','5-10','5-5', 'egi', 'biosemi'. Defaults to '5-10'.
            eeg_reference (str, optional): Reference electrode(s). Defaults to 'single electrode placed on Cz'.
            eeg_ground (str, optional): Location of the ground electrode for the common mode rejection. 
                Defaults to 'placed on CPz'.
            recording_type (str, optional): The type of the recording must be one of: 
                'continuous', 'epoched', 'discontinuous'. Defaults to 'continuous'.

        Returns:
            self(:obj: Convertor): instance of the Convertor class
        """

        if eeg_placement_scheme not in ['10-20','10-10','5-10','5-5', 'egi', 'biosemi']:
            raise ValueError('eeg_placement_scheme must be one of: 10-20, 10-10, 5-10, 5-5, egi or biosemi')
        if recording_type not in ['continuous', 'epoched', 'discontinuous']:
            raise ValueError('recording_type must be one of: continuous, epoched, discontinuous')
        
        eeg, ecg, eog, emg, misc, trig = 0, 0, 0, 0, 0, 0
        for ch in self.raw.info['chs']:
            ch_type = str(ch['kind'])[str(ch['kind']).index('FIF'):-1].split('_')[1]
            if ch_type == 'EEG':
                eeg+=1
            elif ch_type == 'ECG':
                ecg+=1
            elif ch_type == 'EOG':
                eog+=1
            elif ch_type == 'EMG':
                emg+=1
            elif ch_type == 'STIM':
                trig+=1
            else:
                misc+=1
        rec_duration = self.raw.get_data().shape[1]/ self.raw.info['sfreq']
        self.description = {
          "Experimenter": self.raw.info['experimenter'],
          "TaskName":self.BIDSpath.task,
          "TaskDescription":task_description,
          "Instructions":instructions,
          "CogAtlas":cogatlas,
          "InstitutionName":"Weill Cornell Medicine",
          "InstitutionAddress":"525 East 68th Street New York, NY 10065 United States",
          "InstitutionDepartmentName":"Radiology",
          "Manufacturer":manufacturer,
          "ManufacturersModelName":device,
          "DeviceSerialNumber":device_serial_number,
          "SoftwareVersions":software_versions,
          "CapManufacturer":manufacturer,
          "CapManufacturersModelName": self.raw.info['device_info']['type'],
          "EEGChannelCount":eeg,
          "EOGChannelCount":eog,
          "ECGChannelCount":ecg,
          "EMGChannelCount":emg,
          "MiscChannelCount":misc,
          "TriggerChannelCount":trig,
          "PowerLineFrequency":powerline_frequency,
          "SamplingFrequency": self.raw.info['sfreq'],
          "EEGPlacementScheme":eeg_placement_scheme,
          "EEGReference":ref_name,
          "EEGGround":eeg_ground,
          "SoftwareFilters":{
            "Anti-aliasing filter":{
                        "half-amplitude cutoff (Hz)": self.raw.info['lowpass'], 
                        "Roll-off": "6dB/Octave"
                        },
            "High-pass filter": self.raw.info['highpass']
          },
          "HardwareFilters":{
            "ADC's decimation filter (hardware bandwidth limit)":{
              "-3dB cutoff point (Hz)":2000,
            }
          },
          "RecordingDuration":rec_duration,
          "RecordingType":"continuous", 
        }
        if recording_type == 'epoched':
          self.description['EpochLength'] = self.raw.times[1] - self.raw.times[0]
        self.BIDSpath.update(suffix = 'eeg', extension='.json')
        json_obj = json.dumps(self.description)
        with open(f'{self.BIDSpath.fpath}', "w") as outfile:
            outfile.write(json_obj)
        
        print_string = f'sidecar JSON saved in:'
        spaces = ' ' * (39 - len(print_string))
        print(f'{print_string}{spaces}{self.BIDSpath.fpath}')
        
        return self

    def convert_eeg(self):
        self.BIDSpath.update(suffix = 'eeg', extension = '.edf')
        mne.export.export_raw(self.BIDSpath.fpath,self.raw, overwrite=False)
        
        print_string = f'EEG data saved in:'
        spaces = ' ' * (39 - len(print_string))
        print(f'{print_string}{spaces}{self.BIDSpath.fpath}')

        return self