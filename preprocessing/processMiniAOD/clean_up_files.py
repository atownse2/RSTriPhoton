# Use uproot to open files and check if they have the branches we need
import uproot
import os

import signal

import json

import argparse

import sys
sys.path.append('../../')
import analysis.utils.sample_info as si

parser = argparse.ArgumentParser()
parser.add_argument('dType', type=str, help='data, signal, or mc')
parser.add_argument('--year', '-y', type=str, default='2018', help='year')
parser.add_argument('--tuple_version', '-v', default=si.current_tuple_version, type=str, help='tuple version')
parser.add_argument('--test', '-t', action='store_true', help='test mode') 
args = parser.parse_args()

dType = args.dType
year = args.year
tuple_version = args.tuple_version
test = args.test

if tuple_version == 'MiniAOD':
    branch = 'Events'
else:
    branch = 'flattener/tree'

# Check if tuple dir exists
if dType == 'data':
    sample_type = 'data'
elif dType == 'signal':
    sample_type = 'signal'
elif dType == 'GJets':
    sample_type = 'mc'
else:
    raise ValueError(f'Invalid dType {dType}')

tupledir = f'/hadoop/store/user/atownse2/RSTriPhoton/{sample_type}/{year}/{tuple_version}/'
if not os.path.exists(tupledir):
    print(f'Tuple dir does not exist: {tupledir}')
    exit(1)

cache_dir = '/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/preprocessing/processMiniAOD/cache/'
condor_dir = '/scratch365/atownse2/condor'
# Cache good files with last modified time
cache_name = cache_dir + 'clean_up_files.json'
if os.path.exists(cache_name):
    with open(cache_name) as f:
        good_files = json.load(f)
else:
    good_files = {}

# Save good files
def write_good_files():
    if test:
        return
    
    with open(cache_name, 'w') as f:
        json.dump(good_files, f, indent=4)

def signal_handler(sig, frame):
    print('Saving good files...')
    write_good_files()
    exit(0)

signal.signal(signal.SIGINT, signal_handler)

# Get list of files
filelist = [ f for f in os.listdir(tupledir) ]

for i, f in enumerate(filelist):
    try:
        file = tupledir + f
        if f in good_files and os.path.getmtime(file) == good_files[f]:
            continue

        print(f'Checking file {i+1}/{len(filelist)}')
        
        open_file = uproot.open(file)

        if len(open_file.keys()) == 0:
            print(f'File {file} has no branches')
            if test:
                break
            else:
                os.remove(file)
                continue
        elif open_file[branch].num_entries <= 10:
            print(f'File {file} has 10 or less events')
            if test:
                break
            else:
                os.remove(file)
        else:
            good_files[f] = os.path.getmtime(file)
        
        # Remove condor logs
        if not test:
            os.system(f"rm {condor_dir}/log/{f.replace('.root', '.log')}")
            os.system(f"rm {condor_dir}/out/{f.replace('.root', '.out')}")
            os.system(f"rm {condor_dir}/err/{f.replace('.root', '.err')}")

    except Exception as e:
        print(f'Error processing file {file}: {str(e)}')
        if test:
            break
        else:
            os.remove(file)
            continue


write_good_files()