import os

import numpy as np

import hist as bh
import pickle

from coffea import processor
from coffea.nanoevents import NanoEventsFactory, BaseSchema
from coffea.analysis_tools import PackedSelection

from . import sample_info as si


from .logger import Logger
l = Logger()
# from ..workflows import optimize as o

# Define paths
top_dir = '/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton'
hist_dir = f'{top_dir}/hists'

# Define axes
dataset_axis = bh.axis.StrCategory([], name='dataset', label='Dataset', growth=True)
selection_axis = bh.axis.StrCategory([], name='selection', label='Selection', growth=True)

# Diphoton
diphoton_energy_axis = bh.axis.Regular(256, 0, 4000, name='diphoton_energy', label='diphoton energy [GeV]')
diphoton_mass_axis = bh.axis.Regular(1024, 0, 80, name='diphoton_mass', label='diphoton mass [GeV]')
diphoton_moe_axis = bh.axis.Regular(128, 0, 0.05, name='diphoton_moe', label='Diphoton $m/E$')
diphoton_pt_axis = bh.axis.Regular(256, 0, 4000, name='diphoton_pt', label='Diphoton $p_T$ [GeV]')
diphoton_eta_axis = bh.axis.Regular(32, -1.5, 1.5, name='diphoton_eta', label='Diphoton $\eta$')
diphoton_score_axis = bh.axis.Regular(128, 0, 1, name='diphoton_score', label='Diphoton classifier score')
diphoton_isolation_axis = bh.axis.Regular(128, 0, 1, name='diphoton_isolation', label='Diphoton isolation $E_{\gamma \gamma}/\sum_{pf|\Delta R<0.3} E_{pf}$')

# Photon
photon_energy_axis = bh.axis.Regular(256, 0, 4000, name='photon_energy', label='photon energy [GeV]')
photon_pt_axis = bh.axis.Regular(256, 0, 4000, name='photon_pt', label='photon $p_T$ [GeV]')
photon_eta_axis = bh.axis.Regular(32, -4, 4, name='photon_eta', label='photon $\eta$')
photon_sieie_axis = bh.axis.Regular(128, 0, 0.05, name='photon_sieie', label='$\sigma_{i\eta i\eta}$')

# Triphoton
delta_r_axis = bh.axis.Regular(128, 0, 7.5, name='delta_r', label='$\Delta R$')
delta_eta_axis = bh.axis.Regular(64, 0, 4, name='delta_eta', label='$| \Delta \eta |$')
ket_frac_axis = bh.axis.Regular(128, 0, 1, name='ket_frac', label='$| \\vec{p_{T_{\gamma}}} + \\vec{p_{T_{\gamma \gamma}}}| / \sum_\gamma{E}$')
energy_ratio_axis = bh.axis.Regular(128, 0, 2, name='energy_ratio', label='$E_{diphoton}/E_{\gamma}$')
alpha_axis = bh.axis.Regular(512, 0, 1, name='alpha', label='$m_{diphoton}/m_{triphoton}$')
triphoton_mass_axis = bh.axis.Regular(2048, 0, 4000, name='triphoton_mass', label='$m_{3\gamma} [GeV]$')

def get_axis(name, bins):
    if name == 'triphoton_mass':
        return bh.axis.Regular(bins, 0, 4000, name=name, label='$m_{3\gamma} [GeV]$')
    elif name == 'alpha':
        return bh.axis.Regular(bins, 0, 1, name=name, label='$m_{diphoton}/m_{triphoton}$')
    elif name == 'diphoton_mass':
        return bh.axis.Regular(bins, 0, 80, name=name, label='diphoton mass [GeV]')


# Event
pf_met_axis = bh.axis.Regular(128, 0, 500, name='pf_met', label='PF MET [GeV]')
jet_energy_frac_axis = bh.axis.Regular(128, 0, 1, name='jet_energy_frac', label='$\sum_{jets}{E}/\sum_{jets+\gamma}{E}$')

# Jet
jet_eta_axis = bh.axis.Regular(32, -4, 4, name='jet_eta', label='jet $\eta$')

def get_hists(dType, 
              year='2018',
              tuple_version=si.current_tuple_version,
              workflow='analyze',
              remake=False,
              scale=False,
              n_workers=None,
              test=False,
              max_verbosity=0):
    '''Returns a dictionary of histograms for the given data type'''

    if test:
        max_verbosity = 5
    l.set_max_verbosity(max_verbosity)

    l.log(f'Getting histograms for {dType} {year} {tuple_version} {workflow}')
    hist_file = f'{hist_dir}/{dType}_{year}_{tuple_version}_{workflow}.pkl'
    # Check if histograms already exist
    if os.path.exists(hist_file) and not remake:
        l.log("Histograms already exist, loading them now")
        with open(hist_file, 'rb') as f:
            hists = pickle.load(f)
        
        # Check if histograms are scaled
        if scale and (dType == 'signal' or dType == 'GJets'):
            hists = scale_hists(hists, dType, year=year)
    
        return hists
    
    # Otherwise create them
    l.log("Histograms do not exist, creating them now")
    datasets = si.get_datasets(dType, tuple_version=tuple_version, year=year)

    if test:
        key = list(datasets.keys())[0]
        datasets = {key: datasets[key]}
        datasets[key]['files'] = datasets[key]['files'][:1]
        l.log(f'Test mode: only running on {key} with {datasets[key]["files"]}')
    tree = '/flattener/tree'

    if dType == 'signal':
        n_workers = 2
    else:
        n_workers = 32

    executor = processor.futures_executor(workers=n_workers) #64 max
    #executor = processor.iterative_executor()

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

    hists = run(datasets, tree, processor_instance=a_processor())

    # Save histograms
    if not test:
        with open(hist_file, 'wb') as f:
            pickle.dump(hists, f)
    
    return hists

def scale_hists(accumulator, dType,
                year='2018'):

    cutflow = accumulator['cutflow']

    # Check if already scaled, there should not be a dataset axis
    for histname, hist in accumulator.items():
        if histname == 'cutflow':
            continue
        if 'dataset' not in hist.axes.name:
            return accumulator

    xs = {}
    if dType == 'signal':
        subsets = si.get_all_signal_filetags()
        for subset in subsets:
            xs[subset] = si.get_signal_cross_section(subset)['xs']
    else:
        subsets = si.mc_subsets[dType]
        for subset in subsets:
            xs[subset] = si.load_xsdb_info(subset)['xs']
    
    scale = {}
    for subset in subsets:
        n_events = cutflow[subset]['all events']
        scale[subset] = (1000*si.lumi_dict[year]*xs[subset])/n_events
    
    for histname, hist in accumulator.items():
        if histname == 'cutflow':
            continue

        for subset in subsets:
            hist[{'dataset': subset}] = hist[{'dataset': subset}].values()*scale[subset]
        
        if not dType == 'signal':
            accumulator[histname] = hist[{'dataset': sum}]

        # s = hist.stack('dataset')
        # for h in s:
        #     h *= scale[h.name]
        # accumulator[histname] = sum(s)
    
    return accumulator

def hist_to_root(hist, histname):
    import ROOT
    nbins = hist.axes[0].size
    xmin = hist.axes[0].edges[0]
    xmax = hist.axes[0].edges[-1]

    h = ROOT.TH1F(histname, histname, nbins, xmin, xmax)

    counts = hist.counts()
    variances = hist.variances()

    for i in range(nbins):
        h.SetBinContent(i+1, counts[i])
        if variances is not None:
            h.SetBinError(i+1, np.sqrt(variances[i]))

    return h