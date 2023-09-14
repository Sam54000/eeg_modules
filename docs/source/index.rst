EEG modules documentation
===================================

# Introduction
Here are useful modules for EEG analysis using MNE-Python, mne_bids, and pyprep.
These modules have been developped in the specific case of the adult_TBI_attention
project at Sudhin Shah's lab at Weill Cornell Medicine. This package derives 
from the analysis of EEG data during an attention network task (ANT)[^1][^2][^3][^4] recorded in
a large cohort of TBI patients and controls. The EEG data were recorded useg a 128-channel
EGI system (HydroCel Geodesic Sensor Net). Here are also referenced modules for
converting the data from EGI to BIDS format. Feel free to take example from these modules
to adapt them to your own needs.
# Installation
To install the package:
1. Clone the repository on your computer
   `git clone https://github.com/Sam54000/eeg_modules.git`
2. Go to the created `eeg_modules` directory on your computer
3. Install the package using pip
   `python -m pip install -e .`




.. note::

   This project is under active development.

Contents
--------

.. toctree::

   usage
   api
