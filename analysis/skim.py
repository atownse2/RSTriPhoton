import sys
import os

import json

import pandas as pd
import numpy as np
import awkward as ak

import dask
from ndcctools.taskvine import DaskVine

top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(top_dir)

from coffea.nanoevents import NanoEventsFactory, BaseSchema

# from analysis.utils import sample_info as si
from analysis.utils import logger
from analysis.utils import sample_info as si

from analysis import selections as s

l = logger.Logger()

def skim(
        dType, analysis_region,
        years=["2018"],
        scaleout=None, # "DaskVine" or None
        storage="vast", # "vast" or "hadoop"
        verbosity=0,
        test=False):
    '''Loop over files in filenames, fill outputs'''

    l.set_max_verbosity(verbosity)

    # Get fileset
    fileset = si.get_filesets(dType, years=years, storage=storage)

    if test: # Only run over one file per dataset
        fileset = {d: {'files': dd['files'][:1]} for d, dd in fileset.items()}

    # Get outputs
    to_compute = { dataset: get_outputs(fileset['files'], analysis_region) for dataset, fileset in fileset.items() }

    # Compute outputs
    if scaleout == None:
        outputs = dask.compute(to_compute)[0]
    elif scaleout == "distributed":
        from distributed import Client
        with Client() as client:
            outputs = dask.compute(to_compute)[0]
    elif scaleout == "DaskVine":
        from ndcctools.taskvine import DaskVine
        print("Remember to start the factory: vine_factory -T condor --python-env=env.tar.gz -C vine_factory_config.json")
        m = DaskVine(name="triphoton_manager")
        outputs = dask.compute(
            to_compute,
            scheduler=m.get,
            resources={"cores": 1},
            resources_mode=None,
            lazy_transfers=True,
            )[0]
    else:
        raise ValueError("Invalid scheduler")

    write_outputs(outputs, analysis_region, year)

def get_output_filenames(dataset, analysis_region, extension):
    return f"{top_dir}/outputs/{dataset}_{analysis_region}.{extension}"

def get_outputs(files, analysis_region):
    events = NanoEventsFactory.from_root(
        files,
        schemaclass = BaseSchema,
        delayed = True,
    ).events()

    selection = s.EventSelection(events, analysis_region)

    pass_selections = selection.pass_selections
    cutflow = selection.cutflow

    outputs = { 'events' : events[pass_selections], 'cutflow' : cutflow.to_dict()}

    return outputs

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        return super().default(obj)

def write_outputs(outputs, analysis_region, year):
    for dataset, output in outputs.items():
        if output['events'] is not None:
            fout = get_output_filenames(dataset, analysis_region, year, 'parquet')
            ak.to_parquet(output['events'], fout)
        if output['cutflow'] is not None:
            fout = get_output_filenames(dataset, analysis_region, year, 'json')
            with open(fout, 'w') as f:
                json.dump(output['cutflow'], f, cls=NumpyEncoder, indent=4)

def load_outputs(dType, analysis_region, year, scale_mc=False):
    '''Load outputs from file'''
    
    datasets = si.get_datasets(dType)

    parquet_file = get_output_filenames(dType, analysis_region, year, "parquet")
    json_file = get_output_filenames(dType, analysis_region, year, "json")
    if not os.path.exists(parquet_file) or not os.path.exists(json_file):
        # Get input
        do_skims = input(f"Skims for {dType} {analysis_region} {year} do not exist. Do you want to run skims? (y/n) ")
        if do_skims == "y":
            skim(dType, analysis_region, year)
        else:
            return None, None
    
    # Load outputs
    outputs = pd.read_parquet(parquet_file)
    with open(json_file, "r") as f:
        cutflow = json.load(f)
    
    if scale_mc and dType != "data":
        raise NotImplementedError("Scale MC not implemented yet")

    return outputs, cutflow


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("dType", help="Data type to process: data or signal")
    parser.add_argument("analysis_region", help="Analysis region to process: None or preselection")
    parser.add_argument("-y", "--year", type=str, default="2018", help='years to run')
    parser.add_argument("-s", "--scaleout", type=str, default=None, help="Scaleout method None, DaskVine, or distributed")
    parser.add_argument("-o", "--storage", type=str, default="vast", help="Storage method vast or hadoop")
    parser.add_argument("-v", "--verbosity", type=int, default=0, help="Verbosity level")
    parser.add_argument("-t", "--test", action="store_true", help="Run one file per dataset")
    args = parser.parse_args()

    skim(
        args.dType, args.analysis_region, args.year,
        scaleout=args.scaleout,
        storage=args.storage,
        verbosity=args.verbosity,
        test=args.test,
    )