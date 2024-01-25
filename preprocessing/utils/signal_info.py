#Creation and manipulation of the sample information

import os

import re
import json

import numpy as np

top_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from . import tuple_info as ti

signal_tag = 'BkkToGRadionToGGG'

# Metadata functions
def get_sample_info():
    import yaml
    with open(f'{top_dir}/samples.yml') as f:
        sample_info = yaml.load(f, Loader=yaml.FullLoader)
    return sample_info

def write_sample_info(sample_info):
    import yaml
    with open(f'{top_dir}/samples.yml', 'w') as f:
        yaml.dump(sample_info, f)


# Signal functions
def m_moe_to_m_m(m_moe):
    '''Converts the m_moe to m_m'''
    M_BKK = m_moe['M_BKK']
    MOE = m_moe['MOE']
    return {'M_BKK': M_BKK, 'M_R': (M_BKK/2)*MOE}

def get_mass_grid(BKKs, MOEs):
    return [{"M_BKK": M_BKK, "M_R": np.round(MOE*(M_BKK/2), decimals=4)} for M_BKK in BKKs for MOE in MOEs]

mass_grids = {
    'old' : get_mass_grid([180, 250, 500, 1000, 3000], [0.04, 0.02, 0.01, 0.005, 0.0025]),
    'current' : get_mass_grid([180, 250, 500, 1000, 1500, 2000, 2500, 3000], [0.04, 0.02, 0.01, 0.005, 0.0025]),
}

def get_signal_cross_section(filetag):
    xs_file = f'{top_dir}/analysis/metadata/xs/signal_xs.json'
    with open(xs_file, 'r') as f:
        xs = json.load(f)
    return xs[filetag]

def signal_point_tag(signal_point):
    M_BKK = signal_point['M_BKK']
    M_R = signal_point['M_R']
    if M_BKK/int(M_BKK) == 1:
        M_BKK = int(M_BKK)
    if int(M_R) != 0 and M_R/int(M_R) == 1:
        M_R = int(M_R)
    
    tag = f'{signal_tag}_M1-{M_BKK}_R0-{M_R}'.replace('.','p')
    return tag

def signal_dataset(signal_point, year):
    '''Returns the filetag for the signal files
    signal_point : tuple(M_BKK, M_R)
    returns : {signal_filetag}_M1-<M_BKK>_R-<M_R>'''
    dataset = f'{signal_point_tag(signal_point)}_{year}'
    return dataset

def dataset_to_mass(fragment):
    '''Parse the mass point from the fragment'''
    if 'M' not in fragment or 'R0' not in fragment:
        raise ValueError('Fragment does not contain mass point')
    else:
        M1 = float(re.search(r'M1-(\d+p\d+|\d+)', fragment).group(1).replace('p', '.'))
        R0 = float(re.search(r'R0-(\d+p\d+|\d+)', fragment).group(1).replace('p', '.'))
        return {'M_BKK': M1, 'M_R': R0}

def signal_datasets(years, mass_grid_version='old'):
    '''Returns a dictionary of signal datasets that is compatible with the samples.yml file'''

    signal_datasets = {}

    for year in years:
        for signal_point in mass_grids[mass_grid_version]:
            dataset = signal_dataset(signal_point, year)

            # Always include MiniAOD
            MiniAOD_dir = ti.get_tuple_dir('signal', tuple_version='MiniAOD')
            filetag = f'{dataset}_MiniAOD'

            # Check if the data exists
            if not os.path.isdir(MiniAOD_dir) or len([f for f in os.listdir(MiniAOD_dir) if filetag in f]) == 0:
                print("Warning: MiniAOD files do not exist for", dataset)
                continue

            signal_datasets.update({ dataset :{ 'MiniAOD':f'hadoop:{MiniAOD_dir}/{filetag}'}})

    return signal_datasets
