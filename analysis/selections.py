import sys


import awkward as ak

sys.path.append("../")
from analysis.utils import sample_info as s
from analysis import objects as obj

from coffea.nanoevents.methods import nanoaod



class EventSelection:
    def __init__(self,
                 trigger_names = ['HLT_DoublePhoton70_v',
                                  'HLT_DoublePhoton85_v',
                                  'HLT_Photon200_v'],
                 photon_energy_cut = 80,
                 photon_ID_cut = 0, # At least loose ID
                 diphoton_energy_cut = 80,
                 diphoton_score_cut = 0.9,
                 delta_r_cut = 0.15
                 ):

        if isinstance(trigger_names, str):
            self.trigger_names = [self.trigger_names]
        elif isinstance(trigger_names, list):
            self.trigger_names = trigger_names

        self.photon_energy_cut = photon_energy_cut
        self.photon_ID_cut = photon_ID_cut

        self.diphoton_energy_cut = diphoton_energy_cut
        self.diphoton_score_cut = diphoton_score_cut
        
        self.delta_r_cut = delta_r_cut

    def apply_preselection(self,
                           events,
                           cutflow={}
                           ):
        
        self.cut_string = f"""
        Trigger: OR of {self.trigger_names}
        Photon energy > {self.photon_energy_cut}
        Photon cut based ID > {self.photon_ID_cut}
        Diphoton energy > {self.diphoton_energy_cut}
        Diphoton moe > 0
        One pair of photon and diphoton with delta R > {self.delta_r_cut}
        """

        cutflow['all events'] = len(events)

        # Trigger selection

        pass_trigger = ak.zeros_like(events.triggers[:,0])
        trigger_indicies = [s.trigger_indicies[name] for name in self.trigger_names]
        for trigger_index in trigger_indicies:
            pass_trigger = pass_trigger | events.triggers[:,trigger_index]
        
        # Object selection
        photons = obj.get_photons(events)
        diphotons = obj.get_diphotons(events)

        good_photons = (photons.energy > self.photon_energy_cut) & \
                    (photons.cutBasedId > self.photon_ID_cut)

        good_diphotons = (diphotons.energy > self.diphoton_energy_cut) & \
                        (diphotons.ruclu_moe > 0)
        
        

        cutflow['a_good_photon'] = ak.sum(ak.any(good_photons, axis=1))
        cutflow['a_good_diphoton'] = ak.sum(ak.any(good_diphotons, axis=1))
        cutflow['a_good_photon_and_diphoton'] = ak.sum(ak.any(good_photons, axis=1) & \
                                                    ak.any(good_diphotons, axis=1))

        # Candidate selection
        candidates = obj.get_triphoton_candidates(photons[good_photons],
                                                diphotons[good_diphotons])

        # Delta R preselection cut
        candidates = candidates[candidates.delta_r > self.delta_r_cut]

        cutflow['a_candidate'] = ak.sum(ak.num(candidates)>0)

        # Pass preselection
        pass_preselection = pass_trigger & \
                            (ak.num(candidates) > 0)

        cutflow['pass_preselection'] = ak.sum(pass_preselection)

        # Store candidates in events
        events['candidates'] = candidates
        events.behavior = nanoaod.behavior

        return events[pass_preselection==1]
    
    

    