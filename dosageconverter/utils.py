# Here, the helper functions come - does the file exist etc
import os
import logging
import gzip
import sys

def check_lengths(chromosome_numbers, genetics_filenames):
    if len(chromosome_numbers)!= len(genetics_filenames):
        print("List of chromosome numbers has to be the same as the number of genetics filenames")


def check_file_existence(filepath):
    if not os.path.exists(filepath):
        logging.warning("Wrong filename or nonexisting directory or a file: {}".format(filepath))
        sys.exit()


def read_file(file_extension, filepath):
    if file_extension == 'gz':
        read = gzip.open(filepath, 'rt')
        lenght = len(read.readlines())
        read.seek(0)
    else:
        read = open(filepath, 'r')
        lenght = len(read.readlines())
        read.seek(0)
    return read,lenght