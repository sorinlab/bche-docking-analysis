#!/usr/bin/env python
'''Selector script for Nick and Mariel.'''

import argparse
import os
from decimal import Decimal
from shutil import copy

# Function to check the existence of a directory.
# Used in conjunction with argparse to check that the given parameter
# dataset exist.


def valid_file(path):
    '''Function to check the validity of a user-provided path.'''
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(
            '\"{}\" does not exist (must be in the same directory or specify full path).' % path)
    return path

# Initialization of the argument parser.
PARSER = argparse.ArgumentParser(
    description='Selector script for Nick and Mariel. Selction of the \'l\' highest scores between N ICM dataset .log files.',
    epilog="Designed by Xavier Martinez on February 27th, 2016")

PARSER.add_argument('-c', type=valid_file, help='Path to store your results.' )

PARSER.add_argument('-l', type=int, default=10, help='Number of results that you want to see. Default is 10.')

PARSER.add_argument('log_file', type=valid_file, nargs="+",
                    help='The log file(s) to process. Can be 1 or more arguments.')
PARSER.add_argument('-v', action="store_true",
                    default=False, help='Display verbose output.')
args = PARSER.parse_args()

# Initialization of the variables that correspond to the arguments passed
# by the user.
log_arguments = args.log_file
VERBOSITY = args.v
num_Lines = args.l
out_Dir = args.c 

if VERBOSITY:
    print("-- Log files:")
    count = 1
    for log in log_arguments:
        print("\t[{}]\t\"{}\"".format(count, log))
        count += 1

results = []
for log in log_arguments:
    log_file_dir_path = os.path.dirname(os.path.realpath(log))
    with open(log, mode='r') as l:
        lines = [next(l).strip('\n') for x in range(50)]
        if VERBOSITY:
            print("-- First 50 lines of \"{}\".".format(log))
            count = 1
            for line in lines:
                print("\t[{}]\t {}".format(count, line))
                count += 1
        line_processing = [(log, line.split()[0], Decimal(line.split()[1]), log_file_dir_path)
                           for line in lines]
        for value in line_processing:
            results.append(value)
sorted_results = sorted(results, key=lambda item: item[2])
print("-- Results:")
print("\tSmallest entry is from \"{}\": {}\t{}".format(
    sorted_results[0][0], sorted_results[0][1], sorted_results[0][2]))
print("-- 50 Trials with smallest score:")
for x in range(num_Lines):
    print("\t[{}]\t{} \t{} \t{}".format(
        x + 1, sorted_results[x][0], sorted_results[x][1], sorted_results[x][2]))
    if out_Dir is not None:
        hash_count = 0
        while(os.path.exists('{}/{}{}.ob{}'.format(out_Dir, '#' * hash_count, sorted_results[x][1], '#' * hash_count))):
            hash_count += 1
        copy('{}/{}.ob'.format(log_file_dir_path, sorted_results[x][1]), '{}/{}{}.ob{}'.format(out_Dir, '#' * hash_count, sorted_results[x][1], '#' * hash_count))
