import os

import dask
from dask.distributed import Client
from dask_jobqueue import HTCondorCluster

import awkward as ak

import numpy as np

import hist as bh

import pickle
import gzip


from . import histograms as hgm

from ..utils import sample_info as si
from ..utils.logger import Logger

l = Logger()


# Define paths
top_dir = '/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton'
out_dir = f'{top_dir}/output'

def get_outputs(dType, 
                    year='2018',
                    tuple_version=si.get_tuple_version(),
                    workflow='analyze',
                    remake=False,
                    scale=False,
                    cores_per_job=1,
                    n_jobs=1,
                    test=False,
                    write=True,
                    max_verbosity=0):
    '''
    Returns a dictionary of outputs for the given data type, year, and workflow.
    Histograms and cutflows are stored in a pkl file, while arrays are stored in a parquet file.
    '''

    if test:
        max_verbosity = 5
    l.set_max_verbosity(max_verbosity)

    pkl_file = f'{out_dir}/{dType}_{year}_{tuple_version}_{workflow}.pkl.gz'

    # If either file exists, load it
    if os.path.exists(pkl_file) and not remake:
        accumulator = {}

        if os.path.exists(pkl_file):
            l.log("Pkl file exists, loading it now")
            with gzip.open(pkl_file, 'rb') as f:
                accumulator.update(pickle.loads(f.read()))
        
        if scale and (dType == 'signal' or dType == 'GJets'):
            accumulator = hgm.scale_hists(accumulator, dType, year=year)

        return accumulator
    
    # Otherwise create them
    from coffea import processor
    from coffea.nanoevents import NanoEventsFactory, BaseSchema
    
    l.log("Histograms do not exist, creating them now")
    datasets = si.get_filesets(dType, tuple_version=tuple_version, year=year)

    if test:
        key = list(datasets.keys())[0]
        datasets = {key: datasets[key]}
        datasets[key]['files'] = datasets[key]['files'][:1]
        l.log(f'Test mode: only running on {key} with {datasets[key]["files"]}')
    tree = '/flattener/tree'

    executor = processor.iterative_executor()

    if workflow == 'analyze':
        from ..workflows import analyze as a
        a_processor = a.AnalysisProcessor
    elif workflow == 'explore':
        from ..workflows import explore as e
        a_processor = e.ExplorationProcessor
    elif workflow == 'optimize':
        raise NotImplementedError
    elif workflow == 'ruclu':
        from ..workflows import ruclu as r
        a_processor = r.RuCluProcessor
    
    run = processor.Runner(
        executor=executor,
        schema = BaseSchema)

    out = run(datasets, tree, processor_instance=a_processor())

    cluster = HTCondorCluster(
        cores=cores_per_job,
        memory='2GB',
        disk='2GB',
        log_directory=f'{top_dir}/test/logs',
    )

    cluster.scale(jobs=n_jobs)

    client = Client(cluster)

    accumulator = dask.compute(out)[0]
    
    # Save accumulator (only histograms and cutflows for now)
    if write:
        with gzip.open(pkl_file, 'wb') as f:
            f.write(pickle.dumps(accumulator))
    
    return accumulator