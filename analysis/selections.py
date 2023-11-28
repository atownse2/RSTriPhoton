import sys


import awkward as ak

sys.path.append("../")
from analysis.utils import sample_info as s
from analysis import objects as obj

from coffea.nanoevents.methods import nanoaod
from coffea.analysis_tools import PackedSelection

photonID_tags = {1: 'loose',
                 2: 'medium',
                 3: 'tight',
                 0: 'no'}

class EventSelection:
    def __init__(self,
                 events,
                 trigger_names = ['HLT_DoublePhoton70_v',
                                  'HLT_DoublePhoton85_v',
                                  'HLT_Photon200_v'],
                 photon_pt_cut = 80,
                 photon_ID_cut = 1, # At least loose ID
                 diphoton_pt_cut = 80,
                 diphoton_score_cut = 0.9,
                 delta_r_cut = 0.15,
                 ket_frac_max = 0.5,
                 ):
        
        self.selections = PackedSelection()
        self.events = events

        if isinstance(trigger_names, str):
            self.trigger_names = [self.trigger_names]
        elif isinstance(trigger_names, list):
            self.trigger_names = trigger_names

        # Object selection
        self.photon_pt_cut = photon_pt_cut
        self.photon_ID_cut = photon_ID_cut

        self.diphoton_pt_cut = diphoton_pt_cut
        self.diphoton_score_cut = diphoton_score_cut
        
        # Candidate selection
        self.delta_r_cut = delta_r_cut
        self.ket_frac_max = ket_frac_max

    def apply_preselections(self):

        # Trigger selection
        for trigger_name in self.trigger_names:
            trigger_index = s.trigger_indicies[trigger_name]
            self.selections.add(trigger_name, self.events.triggers[:,trigger_index]==1)
        
        # Object selection
        photons = obj.get_photons(self.events)
        diphotons = obj.get_diphotons(self.events)

        good_photons = (photons.pt > self.photon_pt_cut) & \
                    (photons.cutBasedId >= self.photon_ID_cut)

        good_diphotons = (diphotons.pt > self.diphoton_pt_cut) & \
                        (diphotons.ruclu_moe > 0)

        self.selections.add(f'photon_pt{self.photon_pt_cut}_{photonID_tags[self.photon_ID_cut]}ID',
                      ak.num(good_photons)>0)
        self.selections.add(f'diphoton_pt{self.diphoton_pt_cut}',
                      ak.num(good_diphotons)>0)

        # Candidate selection
        candidates = obj.get_triphoton_candidates(photons[good_photons],
                                                diphotons[good_diphotons])

        # Delta R preselection cut
        candidates = candidates[candidates.delta_r > self.delta_r_cut]

        self.selections.add(f'candidate_dr{self.delta_r_cut}', ak.num(candidates)>0)

        candidates = candidates[candidates.ket_frac<self.ket_frac_max]
        candidates = candidates[ak.argsort(candidates.diphoton.dipho_score, axis=1, ascending=False)]
        self.selections.add(f'candidate_ket_frac{self.ket_frac_max}', ak.num(candidates)>0)

        # Pass preselection
        pass_preselection = self.selections.all(*self.selections.names)

        candidates = candidates[pass_preselection]
        candidate = candidates[:,0]
        
        self.events = self.events[pass_preselection]

        # Store candidates in events
        self.events['candidate'] = candidate
        self.events.behavior = nanoaod.behavior


