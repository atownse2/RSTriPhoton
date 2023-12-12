import uproot
import os
import signal
import json
import argparse
import sys

top_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(top_dir)

import analysis.utils.sample_info as si

cache_dir = '/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/preprocessing/cache/'
condor_dir = '/scratch365/atownse2/condor'

# Cache good files with last modified time
cache_name = cache_dir + 'clean_up_files.json'

def clear_condor_logs(jobName):
    logfile = f'{condor_dir}/log/{jobName}.log'
    errfile = f'{condor_dir}/err/{jobName}.err'
    outfile = f'{condor_dir}/out/{jobName}.out'
    for file in [logfile, errfile, outfile]:
        if os.path.exists(file):
            os.remove(file)

def get_tuple_dict():
    if os.path.exists(cache_name):
        with open(cache_name) as f:
            tuple_dict = json.load(f)
    else:
        tuple_dict = {}
    return tuple_dict

def write_good_files(good_files,
                     tuple_version=si.get_tuple_version()):
    tuple_dict = get_tuple_dict()
    tuple_dict[tuple_version] = {'signal': good_files}

    with open(cache_name, 'w') as f:
        json.dump(tuple_dict, f, indent=4)

def get_good_files( group,
                    tuple_version=si.get_tuple_version(),
                    year='2018'
                    ):
    tuple_dict = get_tuple_dict()
    if tuple_version not in tuple_dict:
        tuple_dict[tuple_version] = {group: {}}

    good_files = {}
    tuple_dir = si.get_tuple_dir(group, tuple_version=tuple_version, year=year)
    for file, file_info in tuple_dict[tuple_version][group].items():
        file_path = tuple_dir + file
        if os.path.exists(file_path) and os.path.getmtime(file_path) == file_info['last_modified']:
            good_files[file] = file_info
    return good_files

def clean_up_files(group,
                   year='2018',
                   tuple_version=si.get_tuple_version(),
                   test=False):
    if tuple_version == 'MiniAOD':
        branch = 'Events'
    else:
        branch = 'flattener/tree'

    # Check if tuple dir exists
    if group != 'data' and group != 'signal':
        group = 'mc'

    tuple_dir = si.get_tuple_dir(group, tuple_version=tuple_version, year=year)
    if not os.path.exists(tuple_dir):
        print(f'Tuple dir does not exist: {tuple_dir}')
        exit(1)

    good_files = get_good_files(group, tuple_version=tuple_version, year=year)

    # Get list of files
    filelist = [f for f in os.listdir(tuple_dir)]

    for i, f in enumerate(filelist):
        file = tuple_dir + f
        if f in good_files and os.path.getmtime(file) == good_files[f]['last_modified']:
            continue
        try:
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
                file_info = {'n_events': open_file[branch].num_entries,
                             'last_modified': os.path.getmtime(file)}
                good_files[f] = file_info

        except Exception as e:
            print(f'Error processing file {file}: {str(e)}')

    for file in filelist:
        if file in good_files:
            clear_condor_logs(file.replace('.root',''))
        else:
            print(f'File {file} is not good')
            if not test:
                os.remove(tuple_dir + file)

    write_good_files(good_files, tuple_version=tuple_version)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('group', type=str, help='data, signal, or mc')
    parser.add_argument('--year', '-y', type=str, default='2018', help='year')
    parser.add_argument('--tuple_version', '-v', default=si.get_tuple_version(), type=str, help='tuple version')
    parser.add_argument('--test', '-t', action='store_true', help='test mode')
    parser.add_argument('--delete_cache', action='store_true', help='delete cache')
    args = parser.parse_args()

    if args.delete_cache:
        os.remove(cache_name)
    
    clean_up_files(args.group, args.year, args.tuple_version, args.test)