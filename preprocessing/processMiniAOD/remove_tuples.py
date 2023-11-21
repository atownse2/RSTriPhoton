import os
import sys

sys.path.append('../../')

from analysis.utils import sample_info as s

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('tuple_version', help='Version of the tuple to be removed')
parser.add_argument('--year', '-y', default='2018', help='Year of the tuple to be removed')
parser.add_argument('--test', '-t', help='Print only', action='store_true', default=False)

args = parser.parse_args()
tuple_version = args.tuple_version
year = args.year
test = args.test

for sample_type in s.sample_types:
    dir = s.get_tuple_dir(sample_type, tuple_version, year=year)
    if os.path.isdir(dir):
        print(f'Removing {dir}')
        if not test:
            os.system(f'rm -rf {dir}')