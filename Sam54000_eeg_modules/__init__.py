"""
Here are some simple and useful modules for EEG data analysis and BIDS formatting.
"""
from .bids_conversion import Convertor
from .cli import create_form, console_menu, query_yes_no, prompt_message
from .preprocess_eeg import *
from .eprime_conversion import eprime_data
__all__ = ('Convertor', 
           'create_form', 
           'console_menu', 
           'query_yes_no', 
           'prompt_message', 
           'ref_naming',
           'read_raw_eeg',
           'calculate_rt',
           'preprocess_eeg_resting',
           'preprocess_eeg_ant',
           'extract_resting_block',
           'extract_ant_block',
           'run_prep',
           'apply_custom_baseline',
           'run_ica',
           'read_eeg_and_ica',
           'epochs_stats',
           'parse_extrema',
           'create_ssd_object',
           'apply_ssd',
           'eprime_data',
           )
__version__ = "0.0.20"