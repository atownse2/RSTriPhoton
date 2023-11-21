import sys
import awkward as ak

from coffea import processor
from coffea.analysis_tools import PackedSelection, Weights

import hist as bh

sys.path.append("../../")
import analysis.utils.sample_info as si
import analysis.objects as obj
import analysis.utils.histograms as h
import analysis.selections as selections

photon_mva_cut = {'barrel': {80: 0.42, 90: -0.02},
                  'endcap': {80: 0.14, 90: -0.26}}

class ExplorationProcessor(processor.ProcessorABC):
    def __init__(self):
        pass
    
    def process(self, events):
        dataset = events.metadata['dataset']
        isMC = si.isMC(dataset)
        isSignal = si.isSignal(dataset)

        # Define accumulators
        storage = None #bh.storage.Weight() for when I need to add weights

        output = {'cutflow': {dataset: {}}}
        cutflow = output['cutflow'][dataset]
        
        # Select Events
        EvtSel = selections.EventSelection() # Default cuts for preselection
        events = EvtSel.apply_preselection(events, cutflow)

        ## If difference ket_frac between top two candidates is less than 0.1, take the one with the highest diphoton score
        candidates = events.candidates[abs(events.candidates.ket_frac-events.candidates[:,0].ket_frac)<0.1]
        candidates = candidates[ak.argsort(candidates.diphoton.dipho_score, axis=1, ascending=False)]
        candidate = candidates[:,0]

        photon = candidate.photon
        diphoton = candidate.diphoton
        triphoton = candidate.triphoton

        ## Percent hadronic energy
        jets = obj.get_jets(events)

        jet_not_photon = jets.delta_r(photon) > 0.4
        jet_not_diphoton = jets.delta_r(diphoton) > 0.4

        jets = jets[jet_not_photon & jet_not_diphoton & (jets.pt>30)]

        sum_jet_energy = ak.sum(jets.energy, axis=1)
        jet_energy_frac = sum_jet_energy/(sum_jet_energy+photon.energy+diphoton.energy)

        # Histograms
        output['diphoton_mass'] = bh.Hist(h.dataset_axis, h.diphoton_mass_axis, name='Counts', storage=storage)
        output['diphoton_energy'] = bh.Hist(h.dataset_axis, h.diphoton_energy_axis, name='Counts', storage=storage)
        output['diphoton_pt'] = bh.Hist(h.dataset_axis, h.diphoton_pt_axis, name='Counts', storage=storage)
        output['diphoton_eta'] = bh.Hist(h.dataset_axis, h.diphoton_eta_axis, name='Counts', storage=storage)
        output['diphoton_score'] = bh.Hist(h.dataset_axis, h.diphoton_score_axis, name='Counts', storage=storage)
        output['diphoton_isolation'] = bh.Hist(h.dataset_axis, h.diphoton_isolation_axis, name='Counts', storage=storage)

        output['photon_energy'] = bh.Hist(h.dataset_axis, h.photon_energy_axis, name='Counts', storage=storage)
        output['photon_pt'] = bh.Hist(h.dataset_axis, h.photon_pt_axis, name='Counts', storage=storage)
        output['photon_eta'] = bh.Hist(h.dataset_axis, h.photon_eta_axis, name='Counts', storage=storage)

        output['ket_frac'] = bh.Hist(h.dataset_axis, h.ket_frac_axis, name='Counts', storage=storage)
        output['energy_ratio'] = bh.Hist(h.dataset_axis, h.energy_ratio_axis, name='Counts', storage=storage)
        output['diphoton_moe'] = bh.Hist(h.dataset_axis, h.diphoton_moe_axis, name='Counts', storage=storage)
        output['delta_r'] = bh.Hist(h.dataset_axis, h.delta_r_axis, name='Counts', storage=storage)
        output['delta_eta'] = bh.Hist(h.dataset_axis, h.delta_eta_axis, name='Counts', storage=storage)
        output['triphoton_diphoton_mass'] = bh.Hist(h.dataset_axis, h.triphoton_mass_axis, h.diphoton_mass_axis, name='Counts', storage=storage)
        output['triphoton_diphoton_0p9_mass'] = bh.Hist(h.dataset_axis, h.triphoton_mass_axis, h.diphoton_mass_axis, name='Counts', storage=storage)

        output['pf_met'] = bh.Hist(h.dataset_axis, h.pf_met_axis, name='Counts', storage=storage)
        output['jet_energy_frac'] = bh.Hist(h.dataset_axis, h.jet_energy_frac_axis, name='Counts', storage=storage)

        if isMC:
            output['diphoton_genMatch'] = bh.Hist(h.dataset_axis,
                                                    bh.axis.StrCategory([], name='genMatch', growth=True),
                                                    h.diphoton_mass_axis,
                                                    h.diphoton_isolation_axis,)
            
        output['diphoton_mass'].fill(dataset=dataset, diphoton_mass=diphoton.mass)
        output['diphoton_energy'].fill(dataset=dataset, diphoton_energy=diphoton.energy)
        output['diphoton_pt'].fill(dataset=dataset, diphoton_pt=diphoton.pt)
        output['diphoton_eta'].fill(dataset=dataset, diphoton_eta=diphoton.eta)
        output['diphoton_score'].fill(dataset=dataset, diphoton_score=diphoton.dipho_score)
        output['diphoton_isolation'].fill(dataset=dataset, diphoton_isolation=diphoton.pfIso)

        output['photon_energy'].fill(dataset=dataset, photon_energy=photon.energy)
        output['photon_pt'].fill(dataset=dataset, photon_pt=photon.pt)
        output['photon_eta'].fill(dataset=dataset, photon_eta=photon.eta)

        output['ket_frac'].fill(dataset=dataset, ket_frac=candidate.ket_frac)
        output['energy_ratio'].fill(dataset=dataset, energy_ratio=diphoton.energy/photon.energy)
        output['diphoton_moe'].fill(dataset=dataset, moe=2*diphoton.mass/triphoton.energy) # diphoton energy is expected to be half of resonance energy
        output['delta_r'].fill(dataset=dataset, delta_r=candidate.delta_r)
        output['delta_eta'].fill(dataset=dataset, delta_eta=candidate.delta_eta)
        output['triphoton_diphoton_mass'].fill(dataset=dataset, triphoton_mass=triphoton.mass, diphoton_mass=diphoton.mass)

        a_good_diphoton = diphoton.dipho_score > 0.9
        output['triphoton_diphoton_0p9_mass'].fill(dataset=dataset, triphoton_mass=triphoton[a_good_diphoton].mass, diphoton_mass=diphoton[a_good_diphoton].mass)

        output['pf_met'].fill(dataset=dataset, pf_met=events.met_energy[:,0])
        output['jet_energy_frac'].fill(dataset=dataset, jet_energy_frac=jet_energy_frac)

        # Apply tighter selections
        candidates = candidates[candidates.diphoton.dipho_score > 0.9]
        candidate = candidates[ak.num(candidates)>0][:,0]

        if isMC:
            gen_particles = obj.get_gen_particles(events)[ak.num(candidates)>0]
            gen_particles = gen_particles[gen_particles.pt>70]

            # Matching
            dr_match = 0.1

            if isSignal:
                gen_radions = gen_particles[(gen_particles.pdgId==obj.radion_pdgid) & (gen_particles.mother_pdgId==obj.bkk_pdgid)]
                gen_prompt_photons = gen_particles[(gen_particles.pdgId==22) & (gen_particles.mother_pdgId==obj.bkk_pdgid)]
                gen_others = gen_particles[(gen_particles.pdgId!=obj.radion_pdgid) & (gen_particles.pdgId!=22)]

                diphoton_matched_gen_radions = ak.any(candidate.diphoton.delta_r(gen_radions) < dr_match, axis=1)
                diphoton_matched_gen_prompt_photons = ak.any(candidate.diphoton.delta_r(gen_prompt_photons) < dr_match, axis=1)
                diphoton_matched_gen_others = ak.any(candidate.diphoton.delta_r(gen_others) < dr_match, axis=1)

                output['diphoton_genMatch'].fill(dataset=dataset,
                                                genMatch='radion',
                                                diphoton_mass=candidate.diphoton.mass[diphoton_matched_gen_radions],
                                                diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_radions])
                output['diphoton_genMatch'].fill(dataset=dataset,
                                                    genMatch='prompt_photon',
                                                    diphoton_mass=candidate.diphoton.mass[diphoton_matched_gen_prompt_photons],
                                                    diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_prompt_photons])
                output['diphoton_genMatch'].fill(dataset=dataset,
                                                    genMatch='others',
                                                    diphoton_mass=candidate.diphoton.mass[diphoton_matched_gen_others],
                                                    diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_others])
            else:
                gen_photons = gen_particles[(gen_particles.pdgId==22)]
                gen_electrons = gen_particles[(abs(gen_particles.pdgId)==11)]
                gen_others = gen_particles[(gen_particles.pdgId!=22) & (abs(gen_particles.pdgId)!=11)]

                diphoton_matched_gen_photons = ak.any(candidate.diphoton.delta_r(gen_photons) < dr_match, axis=1)
                diphoton_matched_gen_electrons = ak.any(candidate.diphoton.delta_r(gen_electrons) < dr_match, axis=1)
                diphoton_matched_gen_others = ak.any(candidate.diphoton.delta_r(gen_others) < dr_match, axis=1)

                output['diphoton_genMatch'].fill(dataset=dataset,
                                                    genMatch='photon',
                                                    diphoton_mass=candidate.diphoton.mass[diphoton_matched_gen_photons],
                                                    diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_photons])
                
                output['diphoton_genMatch'].fill(dataset=dataset,
                                                    genMatch='electron',
                                                    diphoton_mass=candidate.diphoton.mass[diphoton_matched_gen_electrons],
                                                    diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_electrons])
                
                output['diphoton_genMatch'].fill(dataset=dataset,
                                                    genMatch='others',
                                                    diphoton_mass=candidate.diphoton.mass[diphoton_matched_gen_others],
                                                    diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_others])

        return output

    def postprocess(self, accumulator):
        pass