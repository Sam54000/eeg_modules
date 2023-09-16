#!/usr/bin/env python3

from simple_term_menu import TerminalMenu
from mne.channels import get_builtin_montages, get_builtin_ch_adjacencies

def main():
    options = get_builtin_montages() + get_builtin_ch_adjacencies()
    terminal_menu = TerminalMenu(options, title="Select the montage used for the EEG experiment (use arros to navigate and Enter to select):)")
    menu_entry_index = terminal_menu.show()
    print(f"You have selected {options[menu_entry_index]}!")

if __name__ == "__main__":
    main()