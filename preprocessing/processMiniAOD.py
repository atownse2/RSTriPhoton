import os
import sys

import htcondor
import argparse

top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(top_dir)

from preprocessing.utils import signal_info as si
from preprocessing.utils import tuple_info as ti
from preprocessing.utils import submit_jobs as sj

test_dir = f'{top_dir}/test/tuples'

filelist_dir = f'{top_dir}/preprocessing/cache/filelists'
ml_dir = f'{top_dir}/preprocessing/tools/CMSSW_10_6_19_patch2/src/'

# We can't write directly to Hadoop so we write to scratch365 first and then mv
mid_dir = f'/scratch365/{os.environ["USER"]}/tmp/'
if not os.path.isdir(mid_dir):
    os.makedirs(mid_dir)

## Functions
# Making filelists using das
def das_query(query, outputfile):
  print('Submitting query: ' + query)
  print(f'Query : {query}')
  os.system('dasgoclient --query "{}" > {}'.format(query, outputfile))

def get_MiniAOD_datasets(dType, years=si.years,
                    dataset=None,
                    mass_grid_version='current',):

    datasets = {}
    if dType == 'signal':
        for year in years:
            if dataset is None:
                for mass_point in si.get_mass_grid(version=mass_grid_version):
                    dataset = si.signal_dataset(mass_point, year)
                    datasets[dataset] = {}
            else:
                datasets[dataset] = {}
        
        input_dir = ti.get_tuple_dir('signal', 'MiniAOD')
        for dataset_name, info in datasets.items():
            files  = [f for f in os.listdir(input_dir) if dataset_name in f]
            info['input_files'] = [ f"file:{input_dir}/{f}" for f in files]

    else:
        # Read from sample yaml
        import yaml
        with open('samples.yml') as f:
            sample_info = yaml.load(f, Loader=yaml.FullLoader)
        
        for dataset_name, info in sample_info['dTypes'][dType].items():
            datasets[dataset_name] = {}
            
            # Look for filelist
            filelist = f'{filelist_dir}/{dataset_name}_filelist.txt'
            if not os.path.isfile(filelist):
                das_query(f"file dataset={info['db'].split(':')[-1]}", filelist)

            # Read filelist into list
            input_files = [ si.nd_redirector+f.replace('\n','') for f in open(filelist).readlines()]
            datasets[dataset_name]['input_files'] = input_files
    
    # Only process 10% of data for now
    if dType == 'data':
        for dataset_name, info in datasets.items():
            info['input_files'] = info['input_files'][:int(len(info['input_files'])*0.1)]


    return datasets

n_files_per_batch = {'data': 25,
                     'GJets': 25,
                     'signal': 10}

executable = f'{top_dir}/preprocessing/tools/scripts/run_ml.sh'
def processMiniAOD(
    dTypes,
    years=si.years,
    subset=None,
    mass_grid_version='current',
    tuple_version=si.get_tuple_version('new'),
    remake=False,
    condor=False,
    test=False
):
    '''Process MiniAOD files into flatTuple files'''

    maxEvents = -1

    if test:
        maxEvents = 10

    for dType in dTypes:
        output_dir = ti.get_tuple_dir(dType, tuple_version)
        datasets = get_MiniAOD_datasets(dType, years=years, subset=subset, mass_grid_version=mass_grid_version)

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        


        for dataset_name, info in datasets.items():
            print(f'Processing {dataset_name}')

            # Split into batches of size n
            input_files = info['input_files']
            n = n_files_per_batch[dType]
            input_file_batches = [input_files[i:i + n] for i in range(0, len(input_files), n)]

            # Process each batch
            for i, files in enumerate(input_file_batches):

                if test:
                    files = files[:2]
                    output_dir = test_dir

                output_file = f"{output_dir}/{dataset_name}_{i}.root"

                # Check if output file already exists
                if os.path.isfile(output_file) and not remake:
                    print(f'Output file already exists: {output_file}')
                    continue

                args = f"{ml_dir} {mid_dir} {','.join(files)} {output_file} {maxEvents} {dType!='data'}"
                if condor:
                    sj.submit_condor(executable, args, f'{dataset_name}_{i}')
                else:
                    sj.submit_local(executable, args)
                
                if test:
                    quit()



        

if __name__ == '__main__':
    ## Inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('--dTypes', nargs='+', default=['signal'], help='Data type to process')
    parser.add_argument('--years', nargs='+', default=['2018'], help='Years to process')
    parser.add_argument('--subset', default=None, help='Subset to process')
    parser.add_argument('--mass_grid_version', default='current', help='Version of mass grid')
    parser.add_argument('--tuple_version', default=si.get_tuple_version('new'), help='Version of tuple')
    parser.add_argument('--remake', action='store_true', help='Remake existing tuples')
    parser.add_argument('--condor', action='store_true', help='Submit jobs to condor')
    parser.add_argument('--test', action='store_true', help='Test mode')

    args = parser.parse_args()

    processMiniAOD(**vars(args))
