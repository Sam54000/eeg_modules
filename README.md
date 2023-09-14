# eprime_conversion.py
Read, parse and convert eprime files into a csv file.
Example of usage:
``` python
import eprime_conversion

# Creating the eprime_data object
filename = '~/proj-project/source/sub-subject01/ses-session01/eeg/eprime_file.txt'
data = eprime_conversion.eprime_data(filename)

# Export the dataframe into a csv format
output_filename = '~/proj-project/source/sub-subject01/ses-session01/eeg/output_file.csv'
data.export2csv()
```

# preprocess_eeg.py
Useful functions to preprocess eeg files recorded with EGI system during ant experiments. These functions are dealing with subtleties of: the EGI system, the ant task, the lab naming conventions (and their variations).
It has been developed when exploring and analyzing data in adult_TBI_attention project.

## ref_naming(raw,desired_ref_name='Cz')
For files recorded with EGI system the name of the reference is not the
same. Sometime it's VertexReference, sometime it's VREF and sometime there
is no reference at all. This little function deal with any reference name
and rename it as desired. If it doesn't exist, create
a flat channel named with the desired ref name (because ref channel in EGI system are flat)

<b><u>Args:</u></b>
	raw ([mne.io.Raw](https://mne.tools/stable/generated/mne.io.Raw.html#mne.io.Raw) object): raw data
    desired_ref_name (str, optional): desired name of the reference. Defaults to 'Cz'.

<b><u>Returns:</u></b>
	raw ([mne.io.Raw](https://mne.tools/stable/generated/mne.io.Raw.html#mne.io.Raw) object): raw data modified in place

## read_raw_eeg(filename, preload=False)
Read eeg data from different format and standardize channel names by removing potential trailing spaces.
Allowed format are:
    * egi (.mff, .RAW)
    * bdf (.bdf)
    * edf (.edf)
    * fif (.fif)
    * eeglab (.set)

<b><u>Args:</u></b>
    filename (str): path to the file
    preload (bool, optional): If True, the data will be preloaded into memory. Defaults to False.
    
<b><u>Returns:</u></b>
    raw ([mne.io.Raw](https://mne.tools/stable/generated/mne.io.Raw.html#mne.io.Raw) object): raw data
    status (str): 'ok' if the file is read correctly, 'corrupted' if the file is corrupted, 'nonexistent' if the file doesn't exists.

Example of usage:
``` python
filename = '~/proj-project/source/sub-subject01/ses-session01/eeg/eeg_file.mff'
raw = read_raw_eeg(filename)
```
## get_stim_chan_name(raw, selection = None)
Get triggers in the alphabetical order and their values from stim channels.

<b><u>Args:</b></u>
    raw ([mne.io.Raw](https://mne.tools/stable/generated/mne.io.Raw.html#mne.io.Raw) object): raw data
    selection (list, optional): list of channels to select. Defaults to None. If channels in selection are not stim channels, they will be ignored.
    
<b><u>Returns:</u></b>
	stim_chan_names (list): list of stim channels names
	stim_chan_data (np.array): array of stim channels data
## preprocess_eeg_ant(raw, montage_name="GSN-HydroCel-129")
Run a pipeline to prepare the raw eeg data obtained during an ant task.

Args:
    raw ([mne.io.Raw](https://mne.tools/stable/generated/mne.io.Raw.html#mne.io.Raw) object): raw data
    montage_name : (str) Default "GSN-HydroCel-129"
        provide the name of the montage of the cap used during the experiment. values has to follow the mne naming standard. For a list of acceptable values run `mne.channels.get_builtin_montages(descriptions=False)`. To get a list of tuples (montage name, description), pass descriptions to True.

Returns:
    raw ([mne.io.Raw](https://mne.tools/stable/generated/mne.io.Raw.html#mne.io.Raw) object): raw data with annotations
    evt_count (dict): dictionnary with the number of events for each event type
## calculate_rt(raw)
## preprocess_eeg_resting(raw, montage_name="GSN-HydroCel-129")
## extract_resting_block(raw)
## extract_ant_blocks(raw)
## run_prep(raw, montage_name="GSN-HydroCel-129")
