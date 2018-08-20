from __future__ import (absolute_import, division, print_function)

import os

allowed_instruments = ["EMU", "MUSR", "CHRONUS", "HIFI"]
allowed_extensions = ["nxs"]


def filter_for_extensions(extensions):
    str_list = ["*." + str(ext) for ext in extensions]
    return "Files (" + ", ".join(str_list) + ")"


def get_instrument_directory(instrument):
    """
    If instrument is supported, returns the directory name on the ISIS network
    in which its data can be found
    """
    if instrument in allowed_instruments:
        instrument_directory = instrument
        if instrument == "CHRONUS":
            instrument_directory = "NDW1030"
        return instrument_directory
    else:
        return None


def check_file_exists(filename):
    return os.path.isfile(filename)


def get_current_run_filename(instrument):
    """
    If instrument is supported, attempts to find the file on the ISIS network which
    contains the data from its current (most up-to-date) run.
    """

    instrument_directory = get_instrument_directory(instrument)
    if instrument_directory is None:
        return ""

    autosave_file = os.sep + os.sep + instrument_directory + os.sep + "data" + os.sep + "autosave.run"
    autosave_points_to = ""
    if not check_file_exists(autosave_file):
        raise ValueError("Cannot find file : " + autosave_file)
    with open(autosave_file, 'r') as f:
        for line in f:
            if len(line.split('.')) == 2:
                autosave_points_to = line
    if autosave_points_to == "":
        current_run_filename = os.sep + os.sep + instrument_directory + os.sep + "data" \
                               + os.sep + instrument_directory + "auto_A.tmp"
    else:
        current_run_filename = os.sep + os.sep + instrument_directory + os.sep + "data" \
                               + os.sep + autosave_points_to
    return current_run_filename


def format_run_for_file(run):
    return "{0:08d}".format(run)


def file_path_for_instrument_and_run(instrument, run):
    """Returns the path to the data file for a given instrument/run"""
    base_dir = os.sep + os.sep + get_instrument_directory(instrument) + os.sep + "data"
    file_name = instrument + format_run_for_file(run) + ".nxs"
    return base_dir.lower() + os.sep + file_name


def parse_input_text_to_files(text, separator=','):
    files = text.split(separator)
    for filename in files:
        try:
            filestem = filename.split(".")[-1]
        except IndexError:
            raise ValueError("Please supply a valid list of files using separator" + separator)
        if filestem not in allowed_extensions:
            raise ValueError("Only " + "".join(allowed_extensions) + " files supported")
    return files


# Fill in digits to make num2 the same as num1
def fill_digits(num1, num2):
    extra_digits = len(str(num1)) - len(str(num2))
    return int(str(num1)[0:extra_digits] + str(num2))


def parse_range_string(text):
    range_list = text.split("-")
    if len(range_list) == 1:
        return [int(range_list[0])]
    if len(range_list) == 2:
        min_range = int(range_list[0])
        max_range = int(range_list[1])
        if min_range > max_range:
            # Assume the second number is truncated
            max_range = fill_digits(min_range, max_range)
        if max_range - min_range > 500:
            raise ValueError("Range from " + str(min_range) + " to " + str(max_range) + " is too large")
        return [min_range + i for i in range(max_range - min_range + 1)]


# return ordered and unique.
def parse_run_string(text):
    runs = []
    ranges = text.split(",")
    for runRange in ranges:
        runs += parse_range_string(runRange)
    return list(set(runs))


def parse_run_list_to_filenames(run_list, instrument):
    files = []
    for run in run_list:
        # pad to 8 zeros
        files += [instrument + '{:08d}'.format(run) + '.nxs']
    return files


def parse_user_input_to_files(input_text, extensions=allowed_extensions):
    input_list = input_text.split(";")
    filenames = []
    for text in input_list:
        # remove whitespace
        if os.path.splitext(text)[-1].lower() in extensions:
            filenames += [text]
    return filenames


def remove_duplicated_files_from_list(file_list):
    # first occurrence of a given file is taken
    files = [os.path.basename(full_path) for full_path in file_list]
    unique_files = [file_list[n] for n, file_name in enumerate(files) if file_name not in files[:n]]
    return unique_files
