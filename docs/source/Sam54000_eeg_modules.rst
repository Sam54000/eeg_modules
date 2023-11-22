The `eeg_modules` package
=======================

.. contents:: :local:

Summary
-------

.. automodule:: Sam54000_eeg_modules

Classes
-------

- `Sam54000_eeg_modules.bids_conversion.Convertor` for converting EEG files and electrode location into BIDS format.
- `Sam54000_eeg_modules.eprime_conversion.eprime_data` for converting E-Prime files into a pandas dataframe.

Global functions and constants
------------------------------

.. autofunction:: Sam54000_eeg_modules.preprocess_eeg.read_raw_eeg
.. autofunction:: Sam54000_eeg_modules.preprocess_eeg.get_stim_chan_name