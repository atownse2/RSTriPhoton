import os
import sys

import htcondor
import argparse

sys.path.append('../../')
from analysis.utils import sample_info as si

condor_dir = '/scratch365/atownse2/condor'
script_dir = os.path.dirname(os.path.realpath(__file__))
test_dir = "/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/test/tuples"

## Inputs
parser = argparse.ArgumentParser()
parser.add_argument('sample_tag', help='Tag to identify subsets of MC/Data to process')
parser.add_argument('--m_moe', nargs=2, type=float, metavar=('M_BKK', 'MOE'), help='Specify a single point in the mass grid to process events for.')
parser.add_argument('--m_m' , nargs=2, type=float, metavar=('M_BKK', 'M_R'), help='Specify a single point in the mass grid to process events for.')
parser.add_argument('--nfiles', '-n', default=1, type=int, help="Number of files to process per subset of dataset")
parser.add_argument('--test', '-t', help='Run 10 events for one file and save to test dir', action='store_true', default=False)
parser.add_argument('--maxEvents', default=-1, type=int, help='Maximum number of events to process per file')
parser.add_argument('--year', '-e', default='2018', help='year to process. Ex: 2018')
parser.add_argument('--batch', '-b', help='Submit jobs to condor', action='store_true', default=False)
parser.add_argument('--remake', '-r', help='Remake filelists, remove old nanoAOD matching dataset tag and year', action='store_true', default=False)
parser.add_argument('--skim-bad-files', '-s', help='Remove bad files from filelists', action='store_true', default=False)
parser.add_argument('--n_max', default=-1, type=int, help="I want this many files total")
parser.add_argument('--print_only', help='Do not start or submit jobs', action='store_true', default=False)

args = parser.parse_args()
test = args.test
sample_tag = args.sample_tag
year = args.year
nfiles_per_subset = args.nfiles
batch = args.batch
remake = args.remake
skim = args.skim_bad_files
maxEvents = args.maxEvents


## Settings
tuple_version = si.current_tuple_version
useXRD = False # Necessary for das but not for signal

## Functions
# Making filelists using das
def das_query(query, outputfile):
  print('Submitting query: ' + query)
  print(f'Query : {query}')
  os.system('dasgoclient --query "{}" > {}'.format(query, outputfile))

# Clearing condor logs
def clear_logs():
    tag = sample_tag
    if sample_tag == 'data':
        if year == '2018':
            tag = 'EGamma'

    os.system(f'rm {condor_dir}/out/{tag}*.out')
    os.system(f'rm {condor_dir}/err/{tag}*.err')
    os.system(f'rm {condor_dir}/log/{tag}*.log')

# Submitting jobs
def submit_batch(inpath, outpath):
    tag = outpath.split('/')[-1].replace('.root','')
    process = htcondor.Submit({
        "executable": f"{script_dir}/run_ml.sh",
        "arguments": f"{inpath} {outpath} {maxEvents}",
        "output": f"{condor_dir}/out/{tag}.out",
        "error" : f"{condor_dir}/err/{tag}.err",
        "log"   : f"{condor_dir}/log/{tag}.log",
        'should_transfer_files': 'no',              
        "request_cpus": "1",
        "request_memory": "128MB",
        "request_disk": "128MB",})

    print("Submitting job for "+tag)
    schedd = htcondor.Schedd()
    if not args.print_only:
        submit_result = schedd.submit(process)

def submit_local(inpath, outpath):
    tag = outpath.split('/')[-1].replace('.root','')
    print("Starting job for " + tag)
    run = f'./run_ml.sh {inpath} {outpath} {maxEvents}'
    if not args.print_only:
        os.system(run)


# Find datasets to process
datasets = {}
isMC = True
if sample_tag == 'signal': # Maybe simplify this by adding the signal datasets to the sample yaml
    if args.m_moe:
        datasets[si.get_signal_filetag(si.m_moe_to_m_m(args.m_moe))] = {}
    elif args.m_m:
        datasets[si.get_signal_filetag(args.m_m)] = {}
    else:
        for file_tag in si.get_all_signal_filetags():
            datasets[file_tag] = {} 
else:
    # Read from sample yaml
    import yaml
    with open('samples.yml') as f:
        sample_info = yaml.load(f, Loader=yaml.FullLoader)

    if sample_tag == 'data':
        if year == '2018':
            sample_tag = 'EGamma'

    for dataset, info in sample_info["samples"].items():
        if sample_tag in dataset and year in dataset:
            datasets[dataset] = {}
            datasets[dataset]['das_dataset'] = info['db'].split(':')[-1]
            if info['group'] == 'data':
                isMC = False
            else:
                isMC = True


# Define output directory
if sample_tag == 'signal':
    sample_type = 'signal'
elif isMC:
    sample_type = 'mc'
else:
    sample_type = 'data'

outdir = si.get_tuple_dir(sample_type, tuple_version, year=year)

# If output directory doesn't exist, make it
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# Test settings
if test:
    nfiles_per_subset = 1
    maxEvents = 10
    outdir = test_dir

# Defining input and output files, XRD formmating can be done simpler if the deepthought redirector works for das as well
for dataset, info in datasets.items():
    if sample_tag == 'signal':
        # Signal is currently stored locally
        input_dir = si.get_tuple_dir('signal', 'MiniAOD', year=year)
        files  = [f for f in os.listdir(input_dir) if dataset in f]

        # Index flatTuple files by order in filelist
        input_files = [ f"{input_dir}/{f}" for f in files]
        output_files = [f'{outdir}/{f.replace("MiniAOD", tuple_version)}' for f in files]

        if useXRD:
            redirector = si.hadoop_redirector
            input_files = [ f.replace('/hadoop', redirector) for f in input_files ]
        else:
            input_files = [ f.replace('/hadoop', 'file:/hadoop') for f in input_files ]

        info['input_files'] = input_files
        info['output_files'] = output_files
    else:
        # Get from DAS
        redirector = si.nd_redirector
        #redirector = si.fnal_redirector

        # If filelist isn't already made, make it
        filelist = f'cache/filelists/{dataset}_filelist.txt'
        if not os.path.isfile(filelist):
            das_query(f"file dataset={info['das_dataset']}", filelist)

        # Read filelist into list
        input_files = [ redirector+f.replace('\n','') for f in open(filelist).readlines()]

        # Index flatTuple files by order in filelist
        output_files = [f'{outdir}/{dataset}_{tuple_version}_{i}.root' for i in range(len(input_files))]

        info['input_files'] = input_files
        info['output_files'] = output_files

## Only using 10% of the data for now
if not isMC:
    for dataset, info in datasets.items():
        info['input_files'] = info['input_files'][:int(len(info['input_files'])*0.1)]
        info['output_files'] = info['output_files'][:int(len(info['output_files'])*0.1)]

## Don't process files that already exist (unless remake is specified)
if not remake:
    for dataset, info in datasets.items():
        to_skip = [ i for i, f in enumerate(info['output_files']) if os.path.isfile(f)]
        n_existing = len(to_skip)
        print(f'Skipping {n_existing} files which are already processed from {dataset} of total {len(info["output_files"])} files')

        info['n_existing_files'] = n_existing
        info['input_files'] = [ f for i, f in enumerate(info['input_files']) if i not in to_skip]
        info['output_files'] = [ f for i, f in enumerate(info['output_files']) if i not in to_skip]

# Remove condor logs
# clear_logs()

# Submitting jobs
for dataset, info in datasets.items():
    n_existing_files = info['n_existing_files']
    print(f'Processing {min(len(info["input_files"]),nfiles_per_subset)} for {dataset}')
    for i, (MiniAOD, flatTuple) in enumerate(zip(info['input_files'], info['output_files'])):
        if i >= nfiles_per_subset:
            break
        if args.n_max > 0 and i+n_existing_files >= args.n_max:
            break

        if batch:
            submit_batch(MiniAOD, flatTuple)
        else:
            submit_local(MiniAOD, flatTuple)
        
        if test:
            quit()