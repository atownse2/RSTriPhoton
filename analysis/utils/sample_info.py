#Creation and manipulation of the sample information
import os
import sys

import json

import preprocessing.utils.signal_info as signal_info
import preprocessing.utils.tuple_info as ti

top_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(top_dir)

meta_dir = f'{top_dir}/analysis/metadata'

# XRootD redirectors
hadoop_redirector = "root://deepthought.crc.nd.edu/"
nd_redirector = "root://ndcms.crc.nd.edu/"

tuple_version = 'NanoAODv9'

# General I/O functions
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

def get_dType(dataset):
    '''Returns the dType of the dataset'''
    sample_info = get_sample_info()

    for dType in sample_info.keys():
        datasets = sample_info[dType].keys()
        if dataset in datasets:
            return dType

# def get_dataset(dType, era, signal_point=None):
#     '''Returns the dataset for the given dType and year'''

#     if dType == 'signal':
#         assert signal_point != None, 'Must specify signal point for signal dataset'

#     sample_info = get_sample_info()
#     for dataset in sample_info[dType].keys():
#         if dType == 'signal' and 
#         if era in dataset:
#             return dataset
    


def get_datasets(dType):
    '''Returns a list of datasets for the given dType'''
    sample_info = get_sample_info()
    return sample_info[dType].keys()

def xrd_list_dir(directory, redirector):
    import XRootD.client
    status, listing = XRootD.client.FileSystem(redirector).dirlist(directory)
    if not status.ok:
        raise RuntimeError(f"Failed to list directory {directory} on {redirector}, make sure you have a proxy.")
    return [f.name for f in listing]

def das_list_dir(dataset):
    filelist_dir = f'{top_dir}/analysis/cache/filelists'
    if not os.path.exists(filelist_dir):
        os.makedirs(filelist_dir)

    cache = f'{filelist_dir}/{dataset[1:].replace("/","_")}.txt'
    if not os.path.exists(cache):
        os.system(f'dasgoclient -query="file dataset={dataset}" > {cache}')
    
    with open(cache) as f:
        return f.read().split('\n')[:-1]


def list_dir(access):
    '''Returns a list of files in the given directory'''

    access_method, path = access.split(':')
    file_dir = os.path.dirname(path)

    if access_method == 'vast':
        return [f for f in os.listdir(file_dir)]
    elif access_method == 'hadoop':
        return xrd_list_dir(file_dir.replace("/hadoop",""), hadoop_redirector)
    elif access_method == 'das':
        return [nd_redirector+f for f in das_list_dir(path)]

def get_filelist(
        dataset_or_dType,
        tuple_version=tuple_version,
        sample_info=get_sample_info(),
        storage='vast',):
    '''Returns a list of files for the given dataset or dType'''

    dType = None
    if dataset_or_dType in sample_info.keys():
        dType = dataset_or_dType
        datasets = sample_info[dataset_or_dType].keys()
    else:
        dType = get_dType(dataset_or_dType)
        if dataset_or_dType in sample_info[dType].keys():
            datasets = [dataset_or_dType]
        else:
            print(f'{dataset_or_dType} not found in sample_info')
            return None
    
    all_files = []
    for dataset in datasets:
        access = sample_info[dType][dataset][tuple_version]

        listed_storage = access.split(':')[0]
        if storage == 'vast' and listed_storage == 'hadoop':
            access = access.replace('hadoop:/hadoop/store/user/atownse2/', 'vast:/project01/ndcms/atownse2/')
        elif storage == 'hadoop' and listed_storage == 'vast':
            access = access.replace('vast:/project01/ndcms/atownse2/', 'hadoop:/hadoop/store/user/atownse2/')

        all_files += list_dir(access)

    return all_files

def get_filesets(dType,
                year=None,
                tuple_version=tuple_version,
                sample_info=get_sample_info(),
                storage='vast'):
    '''Returns a dictionary of datasets and their filelists'''

    datasets = sample_info[dType].keys()
    if year != None:
        datasets = [dataset for dataset in datasets if year in dataset]

    filesets = {
        dataset: {
            'files': get_filelist(
                dataset,
                tuple_version=tuple_version,
                storage=storage,
                sample_info=sample_info
                )
        } for dataset in datasets}
    
    return filesets


## Metadata functions
def get_xs(dataset):
    '''Returns the cross section of the dataset'''
    sample_info = get_sample_info()
    dType = get_dType(dataset)
    return sample_info[dType][dataset]['xs']
