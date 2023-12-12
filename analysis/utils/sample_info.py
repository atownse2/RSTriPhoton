#Creation and manipulation of the sample information

import os

import numpy as np

import re
import json
import yaml

top_dir = '/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton'
meta_dir = f'{top_dir}/analysis/metadata'

condor_dir = '/scratch365/atownse2/condor'
# Metadata
lumi_dict = { '2016': 35.92,
              '2017': 41.53,
              '2018': 59.74} #fb^-1

trigger_indicies = {'HLT_DoublePhoton70_v' : 0,
                    'HLT_DoublePhoton85_v' : 1,
                    'HLT_Photon200_v' : 2} # Should read this from the file...

sample_types = ['data', 'mc', 'signal']

mc_subsets = {'GJets': ['GJets_HT-40To100',
                        'GJets_HT-100To200',
                        'GJets_HT-200To400',
                        'GJets_HT-400To600',
                        'GJets_HT-600ToInf'],}

data_filetags = {'2016' : 'DoubleEG',
                 '2017' : 'DoubleEG',
                 '2018' : 'EGamma'}

signal_filetag = 'BkkToGRadionToGGG'

# Mass Grid
BKK_MASS = [180, 250, 500, 1000, 3000] # Removed 140, 160
MOES = [0.04, 0.02, 0.01, 0.005, 0.0025]

mass_grid = [ (M_BKK, np.round(MOE*(M_BKK/2), decimals=4)) for M_BKK in BKK_MASS for MOE in MOES]

# Adding more mass points but don't want to break workflows that use the old mass grid
new_BKK_MASS = [180, 250, 500, 1000, 1500, 2000, 2500, 3000]
new_mass_grid = [ (M_BKK, np.round(MOE*(M_BKK/2), decimals=4)) for M_BKK in new_BKK_MASS for MOE in MOES]

# XRootD redirectors
hadoop_redirector = "root://deepthought.crc.nd.edu/"
nd_redirector = "root://ndcms.crc.nd.edu/"
fnal_redirector = "root://cmsxrootd.fnal.gov/"

# Metadata functions
tuple_versions = {'MiniAOD': 'MiniAOD',
                 'current': 'FlatAODv4',
                 'new': 'FlatAODv5',
                 'old': 'FlatAODv3',}

def get_tuple_version(version='current'):
    if version in tuple_versions:
        return tuple_versions[version]
    else:
        raise ValueError(f'Invalid version {version}')

def get_tuple_dir(group,
                  tuple_version=get_tuple_version(),
                  year='2018'):
    """
    Group: data, mc, or signal
    Tuple cycle: new, current, old, or MiniAOD
    """
    return f'/hadoop/store/user/atownse2/RSTriPhoton/{group}/{year}/{tuple_version}/'      

def get_mass_grid(version='current'):
    if version == 'current':
        return mass_grid
    elif version == 'new':
        return new_mass_grid
    elif version == 'old':
        return mass_grid
    else:
        raise ValueError(f'Invalid version {version}')

def isMC(dataset):
    '''Returns true if the dataset is MC'''
    if 'data' in dataset:
        return False
    if any([filetag in dataset for filetag in data_filetags.values()]):
        return False
    return True

def isSignal(dataset):
    '''Returns true if the dataset is signal'''
    return signal_filetag in dataset

def get_xs(dataset):
    '''Returns the cross section of the dataset'''
    if isMC(dataset):
        if isSignal(dataset):
            return get_signal_cross_section(dataset)
        else:
            return load_xsdb_info(dataset)['xs']
    else:
        raise ValueError(f'Dataset {dataset} is not MC')

def get_trigger_index(trigger, era):
    """Reads triggers/trigger_names_<era>.txt and returns the index of the trigger"""
    with open(f'{meta_dir}/triggers/triggerNames_{era}.txt', 'r') as tfile:
        trigger_names = tfile.read().splitlines()
    return trigger_names.index(trigger)

def get_n_events(filepath):
    '''Returns the number of events in the file'''
    filename = filepath.split('/')[-1]

    if 'nanoAOD' in filename:
        tree = 'flattenerMatching/tree'
    elif 'MiniAOD' in filename:
        tree = 'Events'
    else:
        raise ValueError('Tier must be nanoAOD or MiniAOD')

    import uproot
    try:
        with uproot.open(filepath) as file:
            return file[tree].num_entries # type: ignore
    except:
        return 0

def get_signal_cross_section(filetag):
    xs_file = f'{top_dir}/analysis/metadata/xs/signal_xs.json'
    with open(xs_file, 'r') as f:
        xs = json.load(f)
    return xs[filetag]

def load_xsdb_info(dataset):
    '''Loads the xsdb info from the json file'''
    xsdb_path = f'{top_dir}/analysis/metadata/xs/'
    xsdb_jsons = [f for f in os.listdir(xsdb_path) if '.json' in f and dataset in f]
    if len(xsdb_jsons) == 0:
        raise ValueError(f'No xsdb json found for {dataset}')
    elif len(xsdb_jsons) > 1:
        raise ValueError(f'Multiple xsdb jsons found for {dataset}')

    with open(xsdb_path+xsdb_jsons[0], 'r') as f:
        info = json.load(f)

    xsdb_info = {}
    for campaign in info:
        if 'EXO-RunIIFall17MiniAODv2' in campaign['MCM']:
            xsdb_info['xs'] = float(campaign['cross_section'])

    return xsdb_info

# Signal functions
def get_signal_filetag(signal_point=None):
    '''Returns the filetag for the signal files
    signal_point : tuple(M_BKK, M_R)
    returns : {signal_filetag}_M1-<M_BKK>_R-<M_R>'''

    if signal_point is not None:
        M_BKK = signal_point[0]
        M_R = signal_point[1]
        if M_BKK/int(M_BKK) == 1:
            M_BKK = int(M_BKK)
        if int(M_R) != 0 and M_R/int(M_R) == 1:
            M_R = int(M_R)

        filetag = f'{signal_filetag}_M1-{M_BKK}_R0-{M_R}'.replace('.','p')
    else:
        filetag = signal_filetag

    return filetag

def m_moe_to_m_m(m_moe):
    '''Converts the m_moe to m_m'''
    M_BKK, MOE = m_moe
    M_R = (M_BKK/2)*MOE
    return (M_BKK, M_R)

def get_mass_point(fragment):
    '''Parse the mass point from the fragment'''
    if 'M' not in fragment or 'R0' not in fragment:
        raise ValueError('Fragment does not contain mass point')
    else:
        M1 = float(re.search(r'M1-(\d+p\d+|\d+)', fragment).group(1).replace('p', '.'))
        R = float(re.search(r'R0-(\d+p\d+|\d+)', fragment).group(1).replace('p', '.'))
        return M1, R

def get_all_signal_points(era='2018', in_mass_grid=True):
    '''Returns all signal points in MiniAOD tier,
    including those that are not in the mass grid
    era: era of dataset default is RunIISummer20UL18'''
    
    if in_mass_grid:
        return mass_grid
    else:
        points = set()
        filepaths = get_filelist('signal', era=era, tier='MiniAOD')
        for filepath in filepaths:
            fragment = filepath.split('/')[-1]
            M1, R = get_mass_point(fragment)
            points.add((M1, R))
        return points

def get_all_signal_filetags(era='2018', in_mass_grid=True):
    '''Returns a list of all signal filetags in mass grid'''
    return [get_signal_filetag(point) for point in get_all_signal_points(era=era, in_mass_grid=in_mass_grid)]

# General I/O functions

def get_das_datasets(group):

    with open('../metadata/samples.yml') as f:
        sample_info = yaml.load(f, Loader=yaml.FullLoader)

    datasets = [ key for key, value in sample_info['samples'].items() if value['group']==group]
    return datasets

def get_filesets(dType,
                year='2018',
                tuple_version=get_tuple_version(),
                signal_in_mass_grid=True,
                useXRD=False, redirector="root://deepthought.crc.nd.edu/"):
    '''Returns a dictionary of datasets and their filelists'''

    datasets = {}
    if dType == 'signal':
        if signal_in_mass_grid:
            mass_points = mass_grid
        else:
            mass_points = get_all_signal_points(year)
        for mass_point in mass_points:
            filetag = get_signal_filetag(mass_point)
            datasets[filetag] = {'files': get_filelist('signal', subset=filetag, tuple_version=tuple_version, year=year, useXRD=useXRD, redirector=redirector)}
    elif dType == 'GJets':
        for dataset in mc_subsets['GJets']:
            datasets[dataset] = {'files': get_filelist('mc', subset=dataset, tuple_version=tuple_version, year=year, useXRD=useXRD, redirector=redirector)}
    elif dType == 'data':
        datasets['data'] = {'files': get_filelist('data', tuple_version, year=year, useXRD=useXRD, redirector=redirector)}
    else:
        raise ValueError(f'Invalid dType {dType}')

    return datasets

def get_filelist(group, subset=None,
                tuple_version=get_tuple_version(),
                signal_point = None,
                year='2018',
                useXRD=False, redirector="root://deepthought.crc.nd.edu/"):
    
    '''Returns a list of files for the given filetag'''

    dir = get_tuple_dir(group, tuple_version=tuple_version, year=year)

    filelist = [f for f in os.listdir(dir)]

    if subset is not None:
        filelist = [f for f in filelist if subset in f]

    if useXRD:
        dir = dir.replace('/hadoop', redirector)

    return [dir+f for f in filelist]

if __name__ == '__main__':
    #Testing

    #Test XRD access
    #redirector = 'root://cmsxrootd.fnal.gov/'
    #redirector = "root://deepthought.crc.nd.edu/"
    # redirector = "root://ndcms.crc.nd.edu/"

    # era = 'RunIISummer20UL18'

    #Test get_filelist
    # print(get_filelist('GJets', era=era))



    #update_signal_metadata(eras=['2016','2017','2018'], test=False)
    #update_bkg_metadata(eras=['2016','2017','2018'], bkg_list=mc_datasets, test=False)
    pass