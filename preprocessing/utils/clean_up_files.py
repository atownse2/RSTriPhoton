import os
import json
import argparse
import sys

import signal

import concurrent.futures

USER = os.environ['USER']
top_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(top_dir)

from preprocessing.utils import tuple_info as ti
cache_dir = f'{top_dir}/preprocessing/cache/'
condor_dir = f'/scratch365/{USER}/condor'

# Cache good files with last modified time
cache_name = cache_dir + 'tuple_info.json'

def clear_condor_logs(jobName):
    logfile = f'{condor_dir}/log/{jobName}.log'
    errfile = f'{condor_dir}/err/{jobName}.err'
    outfile = f'{condor_dir}/out/{jobName}.out'
    for file in [logfile, errfile, outfile]:
        if os.path.exists(file):
            os.remove(file)

def count_events(filepath):
    '''Returns the number of events in the file'''
    filename = filepath.split('/')[-1]

    if 'MiniAOD' in filename:
        tree = 'Events'
    else:
        tree = 'flattener/tree'

    import uproot
    try:
        with uproot.open(filepath) as file:
            return file[tree].num_entries # type: ignore
    except:
        return None

def get_good_files():
    '''Returns a dict of good files'''
    if not os.path.exists(cache_name):
        return {}

    with open(cache_name) as f:
        return json.load(f)

def cache_good_files(good_files):
    with open(cache_name, 'w') as f:
        json.dump(good_files, f, indent=4)

def run_cleaner(dType, tuple_version, **kwargs):
    cleaner = FileCleaner()
    signal.signal(signal.SIGINT, cleaner.signal_handler)
    cleaner.clean_up_files(dType, tuple_version, **kwargs)

class FileCleaner:
    def __init__(self):
        self.good_files = get_good_files()
    
    def signal_handler(self, sig, frame):
        print('Exiting, saving good files')
        cache_good_files(self.good_files)
        sys.exit(0)

    def check_file(self, filepath, file_dict, ask_delete=True):
        file = filepath.split('/')[-1]
        last_modified = os.path.getmtime(filepath)
        if file in file_dict and last_modified == file_dict[file]['last_modified']:
            return None

        n_events = count_events(filepath)
        if n_events == None:
            if ask_delete:
                input(f'Error opening {filepath}. Press enter to delete it.')
            os.remove(file)
            return None
    
        # Clear condor logs and update file info
        clear_condor_logs(file.replace('.root', ''))
        file_info = {'n_events': n_events,
                    'last_modified': last_modified}
        return file_info
        
    def clean_up_files(self,
            dType,
            tuple_version,
            datasets=None, # For cleaning datasets not in samples.yml
            ask_delete=True,
            test=False):

        # Check if tuple dir exists
        tuple_dir = ti.get_tuple_dir(dType, tuple_version)
        if not os.path.exists(tuple_dir):
            print(f'Tuple dir does not exist: {tuple_dir}')
            exit(1)

        if dType not in self.good_files:
            self.good_files[dType] = {}

        if datasets is None:
            datasets = ti.get_datasets(dType)

        for dataset in datasets:
            print(f'Checking {dataset}')
            if dataset not in self.good_files[dType]:
                self.good_files[dType][dataset] = {}
            if tuple_version not in self.good_files[dType][dataset]:
                self.good_files[dType][dataset][tuple_version] = {'files': {}}

            tuple_info = self.good_files[dType][dataset][tuple_version]

            file_dict = tuple_info['files']
            all_files = [f for f in os.listdir(tuple_dir) if dataset in f]

            with concurrent.futures.ThreadPoolExecutor() as executor:
                futures = [executor.submit(self.check_file, f'{tuple_dir}/{file}', file_dict, ask_delete) for file in all_files]
                for future, file in zip(concurrent.futures.as_completed(futures), all_files):
                    file_info = future.result()
                    if file_info is not None:
                        file_dict[file] = file_info
                        print(f'Added {file} to good files with {file_info["n_events"]} events')

            tuple_info['n_files'] = len(file_dict.keys())
            tuple_info['n_events'] = sum([info['n_events'] for info in file_dict.values()])
            tuple_info['files'] = file_dict

            self.good_files[dType][dataset][tuple_version] = tuple_info

        cache_good_files(self.good_files)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dType', type=str, help='data, signal, or GJets')
    parser.add_argument('tuple_version', type=str, help='tuple version')
    parser.add_argument('--test', '-t', action='store_true', help='test mode')
    parser.add_argument('--delete_cache', action='store_true', help='delete cache')
    args = parser.parse_args()

    if args.delete_cache:
        os.remove(cache_name)
    
    run_cleaner(args.dType, args.tuple_version, test=args.test)