import os
import sys

top_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(top_dir)

from preprocessing.utils import signal_info as signal_info

dTypes = ['data', 'GJets', 'signal']
years = ['2018']

local_storage = '/hadoop/store/user/atownse2/RSTriPhoton'

# XRootD redirectors
hadoop_redirector = "root://deepthought.crc.nd.edu/"
nd_redirector = "root://ndcms.crc.nd.edu/"

tuple_versions = {'MiniAOD': 'MiniAOD',
                 'current': 'FlatAODv4',
                 'new': 'FlatAODv5',
                 'old': 'FlatAODv3',}

def get_tuple_version(version='current'):
    if version in tuple_versions:
        return tuple_versions[version]
    else:
        raise ValueError(f'Invalid version {version}')

def get_trigger_index(trigger_name, year):
    '''Returns the index of the trigger'''
    trigger_file = f'{top_dir}/preprocessing/utils/triggers/triggerNames_{year}.txt'
    with open(trigger_file, 'r') as tfile:
        trigger_names = tfile.read().splitlines()
    return trigger_names.index(trigger_name)

def get_tuple_dir(dType,
                  tuple_version):
    return f'{local_storage}/{dType}/{tuple_version}'

def get_sample_info():
    import yaml
    with open(f'{top_dir}/samples.yml') as f:
        sample_info = yaml.load(f, Loader=yaml.FullLoader)
    return sample_info

def write_sample_info(sample_info):
    '''Writes the sample_info to the samples.yml file'''
    import yaml
    with open(f'{top_dir}/samples.yml', 'w') as f:
        yaml.dump(sample_info, f)

def update_local_samples(
    dType,
    tuple_version,
    **signal_kwargs):

    if tuple_version == 'MiniAOD':
        print('Updating MiniAOD not supported')

    sample_info = get_sample_info()

    if dType == 'signal':
        datasets = signal_info.signal_datasets(years, **signal_kwargs)
    else:
        datasets = sample_info[dType]

    # Clean files
    from preprocessing.utils import clean_up_files as clean
    clean.run_cleaner(dType, tuple_version, datasets=datasets)

    for dataset, dataset_info in datasets.items():

        tuple_dir = get_tuple_dir(dType, tuple_version=tuple_version, storage=storage)
        filetag = f'{dataset}_{tuple_version}'

        # Check if the files actually exist
        if not os.path.isdir(tuple_dir) or len([f for f in os.listdir(tuple_dir) if filetag in f]) == 0:
            continue

        dataset_info[tuple_version] = f'{storage}:{tuple_dir}/{filetag}'
    
    sample_info.update({dType: datasets})
    write_sample_info(sample_info)
    return sample_info

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('tuple_version', type=str, help='tuple version')
    parser.add_argument('--storage', '-s', type=str, default='hadoop', help='vast or hadoop')
    parser.add_argument('--mass_grid_version', '-m', type=str, default='current', help='mass grid version')

    args = parser.parse_args()

    tuple_version = args.tuple_version
    storage = args.storage
    mass_grid_version = args.mass_grid_version

    for dType in dTypes:
        update_local_samples(dType, **vars(args))