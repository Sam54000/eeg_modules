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
This is to format the Command Line Interface (CLI) of the package.
"""
# ==============================================================================
#                       IMPORTS STANDARD MODULES AND PACKAGES
# ==============================================================================
import json
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime
import curses
# ==============================================================================
#                         IMPORTS CUSTOM MODULES AND PACKAGES
# ==============================================================================

from simple_term_menu import TerminalMenu
import numpy as np

# ==============================================================================
#                           CLASS AND FUNCTIONS DEFINITIONS
# ==============================================================================

def create_form(stdscr, form_data = None):
    # Initialize the curses library
    curses.curs_set(1)
    stdscr.clear()
    stdscr.refresh()

    # Create an ordered dictionary to store form entries

    # Set the initial cursor position
    current_row = 0

    while True:
        stdscr.clear()
        stdscr.addstr(0, 0, "Press 'q' to quit")

        # Display the form entries
        for i, (field, value) in enumerate(form_data.items()):
            stdscr.addstr(i + 2, 0, f"{field}: {value}")

        # Highlight the currently selected field
        stdscr.addstr(current_row + 2, 0, f"{list(form_data.keys())[current_row]}: {form_data[list(form_data.keys())[current_row]]}", curses.A_BOLD)

        stdscr.refresh()
        curses.echo()

        # Get user input
        key = stdscr.getch()

        # Quit the form if 'q' is pressed
        if key == ord('q'):
            break

        # Move to the next field
        elif key == curses.KEY_DOWN:
            current_row = min(current_row + 1, len(form_data) - 1)

        # Move to the previous field
        elif key == curses.KEY_UP:
            current_row = max(current_row - 1, 0)

        # Edit the selected field
        elif key == curses.KEY_ENTER or key in [10, 13]:
            field = list(form_data.keys())[current_row]
            stdscr.addstr(current_row + 2, 0, f"{field}: ")
            stdscr.refresh()
            value = stdscr.getstr(current_row + 2, len(field) + 2).decode('utf-8')
            form_data[field] = value

    # When the user quits, return the form data as a dictionary
    return dict(form_data)

def console_menu(options = None, title = None):
    """consnole_menu
    This function is used to generate a console menu to select a set of options.

    Args:
        options (list, optional): List of options to select from. Defaults to None.

    Returns:
        output (str): name of the option selected
    """
    terminal_menu = TerminalMenu(options, title=title)
    menu_entry_index = terminal_menu.show()
    output = options[menu_entry_index]
    return output

def query_yes_no(question, default="no"):
    """Ask a yes/no question via raw_input() and return their answer.
    
    Args:
        question (str): A string that is presented to the user.
        default (str): The presumed answer if the user just hits <Enter>.
            It must be "yes" (the default), "no" or None (meaning
            an answer is required of the user).

    Returns:
        choice (str): The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":"yes", "no":"no"}
    if default == None:
        prompt = " [yes/no] "
    elif default == "yes":
        prompt = " [YES/no] "
    elif default == "no":
        prompt = " [yes/NO] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while 1:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return default
        elif choice in valid.keys():
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no'")
    return choice

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