import sys
import awkward as ak

from coffea import processor
from coffea.analysis_tools import PackedSelection, Weights
from coffea.nanoevents.methods import candidate

import hist

top_dir = "/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/"
sys.path.append(top_dir)

import analysis.objects as o
import analysis.utils.histograms as h

photon_mva_cut = {'barrel': {80: 0.42, 90: -0.02},
                  'endcap': {80: 0.14, 90: -0.26}}

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self):
        pass
    
    @property
    def accumulator(self):

        storage = None #hist.storage.Weight() for when I need to add weights

        accumulator = {
            'diphoton_mass': hist.Hist(h.dataset_axis, h.diphoton_mass_axis, name='Counts', storage=storage),
            'diphoton_energy': hist.Hist(h.dataset_axis, h.diphoton_energy_axis, name='Counts', storage=storage),
            'diphoton_pt' : hist.Hist(h.dataset_axis, h.diphoton_pt_axis, name='Counts', storage=storage),
            'diphoton_eta': hist.Hist(h.dataset_axis, h.diphoton_eta_axis, name='Counts', storage=storage),
            'diphoton_isolation' : hist.Hist(h.dataset_axis, h.diphoton_isolation_axis, name='Counts', storage=storage),

            'photon_energy': hist.Hist(h.dataset_axis, h.photon_energy_axis, name='Counts', storage=storage),
            'photon_pt': hist.Hist(h.dataset_axis, h.photon_pt_axis, name='Counts', storage=storage),
            'photon_eta': hist.Hist(h.dataset_axis, h.photon_eta_axis, name='Counts', storage=storage),

            'pt_frac' : hist.Hist(h.dataset_axis, h.pt_frac_axis, name='Counts', storage=storage),
            'mass_ratio': hist.Hist(h.dataset_axis, h.mass_ratio_axis, name='Counts', storage=storage),
            'delta_r': hist.Hist(h.dataset_axis, h.delta_r_axis, name='Counts', storage=storage),
            'delta_eta' : hist.Hist(h.dataset_axis, h.delta_eta_axis, name='Counts', storage=storage),
            'triphoton_diphoton_mass': hist.Hist(h.dataset_axis, h.triphoton_mass_axis, h.diphoton_mass_axis, name='Counts', storage=storage),

            'pf_met' : hist.Hist(h.dataset_axis, h.pf_met_axis, name='Counts', storage=storage),
            'pc_hadronic_energy' : hist.Hist(h.dataset_axis, h.pc_hadron_energy_axis, name='Counts', storage=storage),

            'cutflow' : {}
        }
        return accumulator

    def process(self, events):
        dataset = events.metadata['dataset']
        #selection = PackedSelection()
        weight = 1 #Weights(len(events))
        output = self.accumulator

        # Cutflow
        cutflow = output['cutflow']
        cutflow['all events'] = len(events)

        # Trigger/Lumi/Filters
        #pass_triggers = events.triggers[:,self.trigger_index]
        trigger_DoublePhoton70 = events.triggers[:,0]
        trigger_DoublePhoton85 = events.triggers[:,1]
        trigger_Photon200 = events.triggers[:,2]

        cutflow['trigger_DoublePhoton70'] = ak.sum(trigger_DoublePhoton70)
        cutflow['trigger_DoublePhoton85'] = ak.sum(trigger_DoublePhoton85)
        cutflow['trigger_Photon200'] = ak.sum(trigger_Photon200)

        # Make objects
        photons = o.get_photons(events)
        diphotons = o.get_diphotons(events)

        #Make selections 2
        energy_cut = 80
        dipho_cut = 0.9

        good_diphotons = (diphotons.energy > energy_cut) & \
                        (diphotons.dipho_score > dipho_cut) & \
                        (diphotons.ruclu_moe > 0)
        a_good_diphoton = ak.any(good_diphotons, axis=1)

        cutflow['diphoton_pt80_dipho90'] = ak.sum(a_good_diphoton)

        good_photons = (photons.energy > energy_cut) & \
                        (photons.cutBasedId >= 1)

        a_good_photon = ak.any(good_photons, axis=1)

        cutflow['photon_pt80_loose'] = ak.sum(a_good_photon)
        cutflow['photon_pt80_loose_diphoton_pt80_dipho90'] = ak.sum(a_good_photon & a_good_diphoton)

        # Delta R cut
        combinations = ak.cartesian( {'photon': photons[good_photons],
                                    'diphoton': diphotons[good_diphotons]})

        delta_r = combinations.photon.delta_r(combinations.diphoton)
        combinations = combinations[delta_r > 1.5]
        cutflow['deltaR_1p5'] = ak.sum(ak.num(combinations)>0)

        # Effect of trigger
        pass_preselection = ak.num(combinations)>0
        cutflow['pass_preselection'] = ak.sum(pass_preselection)
        cutflow['pass_preselection_DoublePhoton70'] = ak.sum(pass_preselection & trigger_DoublePhoton70)
        cutflow['pass_preselection_DoublePhoton85'] = ak.sum(pass_preselection & trigger_DoublePhoton85)
        cutflow['pass_preselection_Photon200'] = ak.sum(pass_preselection & trigger_Photon200)

        # Choose DoublePhoton70 trigger for now
        combinations = combinations[trigger_DoublePhoton70==1]

        # Sort combinations by diphoton score
        combinations = combinations[ak.argsort(combinations.diphoton.dipho_score, ascending=False)]
        combination = combinations[ak.num(combinations) > 0][:,0]

        photon = combination.photon
        diphoton = combination.diphoton
        triphoton = photon + diphoton

        # Fill histograms
        output['diphoton_mass'].fill(dataset=dataset, diphoton_mass=diphoton.ruclu_mass, weight=weight)
        output['diphoton_energy'].fill(dataset=dataset, diphoton_energy=diphoton.energy, weight=weight)
        output['diphoton_pt'].fill(dataset=dataset, diphoton_pt=diphoton.pt, weight=weight)
        output['diphoton_eta'].fill(dataset=dataset, eta=diphoton.eta, weight=weight)
        output['photon_energy'].fill(dataset=dataset, photon_energy=photon.energy, weight=weight)
        output['photon_pt'].fill(dataset=dataset, photon_pt=photon.pt, weight=weight)
        output['photon_eta'].fill(dataset=dataset, eta=photon.eta, weight=weight)
        output['mass_ratio'].fill(dataset=dataset, mass_ratio=diphoton.ruclu_mass/triphoton.mass, weight=weight)
        output['dr'].fill(dataset=dataset, dr=ak.flatten(delta_r[delta_r>0.1]), weight=weight)
        output['triphoton_diphoton_mass'].fill(dataset=dataset, triphoton_mass=triphoton.mass, diphoton_mass=diphoton.ruclu_mass, weight=weight)

        return output
    
    def postprocess(self, accumulator):
        pass
