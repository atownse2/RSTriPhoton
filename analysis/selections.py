import sys

from collections import OrderedDict

# import awkward as ak
import dask_awkward as dak

# sys.path.append("../")
from .utils import sample_info as s
from . import objects as obj

from coffea.nanoevents.methods import nanoaod
from coffea.analysis_tools import PackedSelection

analysis_selections = {}

analysis_selections['preselection'] = {
    'photon_pt' : 60,
    'photon_ID' : 1, # At least loose ID
    'diphoton_pt' : 80,
    'diphoton_score' : 0.1,
    'delta_r' : 0.1,
}

analysis_selections['signal'] = analysis_selections['preselection']
analysis_selections['signal'].update({
    'ket_ratio' : 0.6, # Just explaining the concept of selection inheritance
    })

def get_selections(analysis_region, **modified_selections):

    if analysis_region is None:
        return {}
    elif analysis_region not in analysis_selections:
        raise ValueError(f"Analysis region {analysis_region} not recognized.")
    else:
        selections = analysis_selections[analysis_region]
        for key, value in modified_selections.items():
            selections[key] = value
        return selections

class Cutflow:
    def __init__(self):
        self.cutflow = OrderedDict()

    def __add__(self, other):
        for key in other.cutflow:
            if key not in self.cutflow:
                self.cutflow[key] = 0
            self.cutflow[key] += other.cutflow[key]
        return self

    def __getitem__(self, key):
        if key not in self.cutflow:
            self.cutflow[key] = 0
        return self.cutflow[key]

    def __setitem__(self, key, value):
        self.cutflow[key] = value

    def to_dict(self):
        return self.cutflow
    
    def from_dict(self, cutflow):
        self.cutflow = cutflow

class ObjectSelection:

    def __init__(self, objects, cutflow):

        self.cutflow = cutflow
        self.pass_selections = dak.ones_like(objects.pt, dtype=bool)
    
    def select(self, selection_name, selection):
        if selection_name is not None:
            self.cutflow[selection_name] = dak.sum(dak.num(selection)>0)

        if self.pass_selections is None:
            self.pass_selections = selection
        else:
            self.pass_selections = self.pass_selections & selection
        
        return self

class EventSelection:

    def __init__(
        self,
        events, 
        analysis_region,
        ):

        self.cutflow = Cutflow()
        self.pass_selections = dak.ones_like(events.run, dtype=bool)

        self.selections = get_selections(analysis_region)

        self.select_events(events)
        self.calculate_quantities(events)

    def select(self, selection_name, selection):
        if selection_name is not None:
            self.cutflow[selection_name] = dak.sum(selection)

        if self.pass_selections is None:
            self.pass_selections = selection
        else:
            self.pass_selections = self.pass_selections & selection
        
        return self

    def select_events(self, events):

        self.cutflow['total'] = dak.num(events, axis=0)
        
        ## Preselection
        photons = obj.get_photons(events)
        photonSelection = ObjectSelection(photons, self.cutflow)
        if 'photon_pt' in self.selections:
            pt_cut = self.selections['photon_pt']
            photonSelection.select(f'photon_pt_{pt_cut}', photons.pt > pt_cut)
        if 'photon_ID' in self.selections:
            ID_cut = self.selections['photon_ID']
            photonSelection.select(f'photon_{obj.photonID_tags[ID_cut]}ID', photons.cutBasedId >= ID_cut)
        photons = photons[photonSelection.pass_selections]

        diphotons = obj.get_diphotons(events)
        diphotonSelection = ObjectSelection(diphotons, self.cutflow).select(None, (diphotons.ruclu_moe > 0))
        if 'diphoton_pt' in self.selections:
            pt_cut = self.selections['diphoton_pt']
            diphotonSelection.select(f'diphoton_pt_{pt_cut}', diphotons.pt > pt_cut)
        if 'diphoton_score' in self.selections:
            score_cut = self.selections['diphoton_score']
            diphotonSelection.select(f'diphoton_score_{score_cut}', diphotons.dipho_score > score_cut)
        diphotons = diphotons[diphotonSelection.pass_selections]

        candidates = obj.get_triphoton_candidates(photons, diphotons)
        candidateSelection = ObjectSelection(candidates.triphoton, self.cutflow).select(None, candidates.triphoton.mass > 0)
        if 'delta_r' in self.selections:
            dr_cut = self.selections['delta_r']
            candidateSelection.select(f'delta_r_{dr_cut}', candidates.delta_r > dr_cut)
        candidates = candidates[candidateSelection.pass_selections]

        self.select('candidate_pass', dak.num(candidates) > 0)


        ## Other selections ...


        ## Store quantities in events
        events['candidates'] = candidates
    

    def calculate_quantities(self, events):
        pass
