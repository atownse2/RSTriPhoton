import os
import sys
import yaml

top_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(top_dir)

# Make a das query to get the list of files
# dasgoclient -query="file dataset=/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD"

def das_query(query):
    return os.popen('dasgoclient -query="{}"'.format(query)).read().split('\n')[:-1]

def get_query(L1Name):
    return f"dataset dataset=/{L1Name}/*UL*MiniAODv2_NanoAODv9_v*/NANOAOD*"

def format_dataset(dataset):
    dataset_name = dataset.split('-')[0][1:].replace('/','_')
    return {dataset_name: {'NanoAODv9': f"das:{dataset}"}}

if __name__ == "__main__":

    samples = {}
    samples['data'] = {}
    for year, L1 in zip(['2017', '2018'],["DoubleEG", "EGamma"]):
        datasets = das_query(get_query(L1))
        for dataset in datasets:
            if year in dataset:
                samples['data'].update(format_dataset(dataset))
    
    with open(f'{top_dir}/samples.yml', 'w') as f:
        yaml.dump(samples, f, default_flow_style=False)