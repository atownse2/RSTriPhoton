import sys

# import awkward as ak
import dask_awkward as dak

from coffea import processor
from coffea.nanoevents.methods import nanoaod
from coffea.analysis_tools import PackedSelection, Weights

# import hist as bh
import hist.dask as dh
import hist as bh

sys.path.append("../../")
import analysis.utils.sample_info as si

import analysis.processing.histograms as hg
import analysis.processing.objects as obj

class ExplorationProcessor(processor.ProcessorABC):
    def __init__(self,
                trigger_names = ['HLT_DoublePhoton70_v',
                                'HLT_DoublePhoton85_v',
                                'HLT_Photon200_v'],
                photon_pt_cut = 80,
                photon_ID_cut = 1, # At least loose ID
                diphoton_pt_cut = 80,
                diphoton_score_cut = 0.9,
                delta_r_cut = 0.15,
                ):
        
        self._trigger_names = trigger_names
        if isinstance(trigger_names, str):
            self.trigger_names = [self.trigger_names]

        self.photon_pt_cut = photon_pt_cut
        self.photon_ID_cut = photon_ID_cut
        self.photon_ID = obj.photonID_tags[self.photon_ID_cut]

        self.diphoton_pt_cut = diphoton_pt_cut
        self.diphoton_score_cut = diphoton_score_cut

        self.delta_r_cut = delta_r_cut

    def make_output(self):
        storage = None #dh.storage.Weight() for when I need to add weights

        output = {}

        # Cutflow
        output['cutflow'] = {self.dataset: {}}

        # Histograms
        output['diphoton_mass'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_mass_axis, name='Counts', storage=storage)
        output['diphoton_energy'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_energy_axis, name='Counts', storage=storage)
        output['diphoton_pt'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_pt_axis, name='Counts', storage=storage)
        output['diphoton_eta'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_eta_axis, name='Counts', storage=storage)
        output['diphoton_score'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_score_axis, name='Counts', storage=storage)
        output['diphoton_isolation'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_isolation_axis, name='Counts', storage=storage)

        output['photon_energy'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.photon_energy_axis, name='Counts', storage=storage)
        output['photon_pt'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.photon_pt_axis, name='Counts', storage=storage)
        output['photon_eta'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.photon_eta_axis, name='Counts', storage=storage)

        output['ket_frac'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.ket_frac_axis, name='Counts', storage=storage)
        output['energy_ratio'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.energy_ratio_axis, name='Counts', storage=storage)
        output['alpha'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.alpha_axis, name='Counts', storage=storage)
        output['delta_r'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.delta_r_axis, name='Counts', storage=storage)
        output['delta_eta'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.delta_eta_axis, name='Counts', storage=storage)

        output['triphoton_mass'] = dh.Hist(hg.dataset_axis, hg.selection_axis,
                                            hg.get_axis('triphoton_mass', 2048),
                                            name='Counts', storage=storage)
        output['diphoton_mass'] = dh.Hist(hg.dataset_axis, hg.selection_axis,
                                          hg.get_axis('diphoton_mass', 1024),
                                          name='Counts', storage=storage)

        # Reduce binning so that memory doesn't explode
        output['triphoton_diphoton_mass'] = dh.Hist(hg.dataset_axis, hg.selection_axis,
                                                    hg.get_axis('triphoton_mass', 512),
                                                    hg.get_axis('diphoton_mass', 256),
                                                    name='Counts', storage=storage)

        output['pf_met'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.pf_met_axis, name='Counts', storage=storage)
        output['jet_energy_frac'] = dh.Hist(hg.dataset_axis, hg.selection_axis, hg.jet_energy_frac_axis, name='Counts', storage=storage)

        if self.isMC:
            output['diphoton_genMatch'] = dh.Hist(hg.dataset_axis,
                                                  hg.selection_axis, 
                                                  bh.axis.StrCategory([], name='genMatch', growth=True),
                                                  hg.diphoton_isolation_axis,)

        return output

    def process(self, events):

        self.dataset = events.metadata['dataset']
        self.isMC = si.isMC(self.dataset)
        self.isSignal = si.isSignal(self.dataset)

        # Define accumulators
        output = self.make_output()
        cutflow = output['cutflow'][self.dataset]

        # Store cuts in histogram metadata
        selections = {}

        cutflow['all'] = dak.num(events, axis=0)
        # Do preselection

        ## Create objects (maybe not needed if schema is updated)
        photons = obj.get_photons(events)
        diphotons = obj.get_diphotons(events)

        pho_pass = (photons.pt > self.photon_pt_cut) & \
                    (photons.cutBasedId >= self.photon_ID_cut)

        dipho_pass = (diphotons.pt > self.diphoton_pt_cut) & \
                        (diphotons.ruclu_moe > 0)

        photons = photons[pho_pass]
        diphotons = diphotons[dipho_pass]

        ## Create candidates
        candidates = obj.get_triphoton_candidates(photons, diphotons)

        pass_dr = candidates.delta_r > self.delta_r_cut

        candidates = candidates[pass_dr]

        ## Cutflow
        cutflow[f'pho_pass_pt{self.photon_pt_cut}_{self.photon_ID}ID'] = dak.sum(dak.num(photons)>0) # Number of events with at least one photon passing
        cutflow[f'dipho_pass_pt{self.diphoton_pt_cut}'] = dak.sum(dak.num(diphotons)>0) # Number of events with at least one diphoton passing
        cutflow[f'pass_dr{self.delta_r_cut}'] = dak.sum(dak.num(candidates)>0) # Number of events with at least one candidate passing dr cut
        cutflow[f'pass_preselection'] = dak.sum(dak.num(candidates)>0) # Number of events with at least one candidate passing preselection

        ## Store preselection cuts in metadata
        selection = 'preselection'
        selections[selection] = f'photon_pt{self.photon_pt_cut}_{self.photon_ID}ID, diphoton_pt{self.diphoton_pt_cut}, candidate_dr{self.delta_r_cut}'

        ## Apply preselection
        pass_preselection = dak.num(candidates) > 0

        events = events[pass_preselection]
        candidates = candidates[pass_preselection] # Already sorted by diphoton score
        
        # Calculate quantities

        ## Get objects
        candidate = candidates[:,0] # 0 = Diphoton has highest diphoton score
        photon = candidate.photon
        diphoton = candidate.diphoton
        triphoton = candidate.triphoton

        ## Percent hadronic energy
        jets = obj.get_jets(events)

        jet_not_photon = jets.delta_r(photon) > 0.4
        jet_not_diphoton = jets.delta_r(diphoton) > 0.4

        jets = jets[jet_not_photon & jet_not_diphoton & (jets.pt>30)]

        sum_jet_energy = dak.sum(jets.energy, axis=1)
        jet_energy_frac = sum_jet_energy/(sum_jet_energy+photon.energy+diphoton.energy)

        # Fill histograms
        output['diphoton_mass'].fill(dataset=self.dataset, selection=selection, diphoton_mass=diphoton.mass)
        output['diphoton_energy'].fill(dataset=self.dataset, selection=selection, diphoton_energy=diphoton.energy)
        output['diphoton_pt'].fill(dataset=self.dataset, selection=selection, diphoton_pt=diphoton.pt)
        output['diphoton_eta'].fill(dataset=self.dataset, selection=selection, diphoton_eta=diphoton.eta)
        output['diphoton_score'].fill(dataset=self.dataset, selection=selection, diphoton_score=diphoton.dipho_score)
        output['diphoton_isolation'].fill(dataset=self.dataset, selection=selection, diphoton_isolation=diphoton.pfIso)

        output['photon_energy'].fill(dataset=self.dataset, selection=selection, photon_energy=photon.energy)
        output['photon_pt'].fill(dataset=self.dataset, selection=selection, photon_pt=photon.pt)
        output['photon_eta'].fill(dataset=self.dataset, selection=selection, photon_eta=photon.eta)

        output['ket_frac'].fill(dataset=self.dataset, selection=selection, ket_frac=candidate.ket_frac)
        output['energy_ratio'].fill(dataset=self.dataset, selection=selection, energy_ratio=diphoton.energy/photon.energy)
        output['alpha'].fill(dataset=self.dataset, selection=selection, alpha=diphoton.mass/triphoton.mass)
        output['delta_r'].fill(dataset=self.dataset, selection=selection, delta_r=candidate.delta_r)
        output['delta_eta'].fill(dataset=self.dataset, selection=selection, delta_eta=candidate.delta_eta)

        output['triphoton_mass'].fill(dataset=self.dataset, selection=selection, triphoton_mass=triphoton.mass)
        output['triphoton_diphoton_mass'].fill(dataset=self.dataset, selection=selection, triphoton_mass=triphoton.mass, diphoton_mass=diphoton.mass)

        output['pf_met'].fill(dataset=self.dataset, selection=selection, pf_met=events.met_energy[:,0])
        output['jet_energy_frac'].fill(dataset=self.dataset, selection=selection, jet_energy_frac=jet_energy_frac)

        # Apply additional selections

        ## Loose diphoton score cut
        pass_dipho_score = candidates.diphoton.dipho_score > self.diphoton_score_cut
        candidates = candidates[pass_dipho_score]

        ## Cutflow
        cutflow[f'candidate_dipho_score{self.diphoton_score_cut}'] = dak.sum(dak.num(candidates)>0) # Number of events with at least one candidate passing diphoton score cut

        ## Store cuts in histogram metadata
        selection = 'preselection + diphoton score'
        selections[selection] = selections['preselection'] + f', diphoton_score{self.diphoton_score_cut}'

        ## Apply additional selections
        mask = dak.num(candidates) > 0

        events = events[mask]
        candidates = candidates[mask]


        # Get objects
        candidate = candidates[:,0] # 0 = Diphoton has highest diphoton score
        photon = candidate.photon
        diphoton = candidate.diphoton
        triphoton = candidate.triphoton


        # Fill histograms
        output['diphoton_mass'].fill(dataset=self.dataset, selection=selection, diphoton_mass=diphoton.mass)
        output['triphoton_mass'].fill(dataset=self.dataset, selection=selection, triphoton_mass=triphoton.mass)
        output['triphoton_diphoton_mass'].fill(dataset=self.dataset, selection=selection, triphoton_mass=triphoton.mass, diphoton_mass=diphoton.mass)

        if self.isMC:
            gen_particles = obj.get_gen_particles(events)
            gen_particles = gen_particles[gen_particles.pt>70]

            # Matching
            dr_match = 0.1

            if self.isSignal:
                gen_radions = gen_particles[(gen_particles.pdgId==obj.radion_pdgid) & (gen_particles.mother_pdgId==obj.bkk_pdgid)]
                gen_prompt_photons = gen_particles[(gen_particles.pdgId==22) & (gen_particles.mother_pdgId==obj.bkk_pdgid)]
                gen_others = gen_particles[(gen_particles.pdgId!=obj.radion_pdgid) & (gen_particles.pdgId!=22)]

                diphoton_matched_gen_radions = dak.any(diphoton.delta_r(gen_radions) < dr_match, axis=1)
                diphoton_matched_gen_prompt_photons = dak.any(diphoton.delta_r(gen_prompt_photons) < dr_match, axis=1)
                diphoton_matched_gen_others = dak.any(diphoton.delta_r(gen_others) < dr_match, axis=1)

                output['diphoton_genMatch'].fill(dataset=self.dataset, selection=selection,
                                                genMatch='radion',
                                                diphoton_isolation=diphoton.pfIso[diphoton_matched_gen_radions])
                output['diphoton_genMatch'].fill(dataset=self.dataset, selection=selection,
                                                    genMatch='prompt_photon',
                                                    diphoton_isolation=diphoton.pfIso[diphoton_matched_gen_prompt_photons])
                output['diphoton_genMatch'].fill(dataset=self.dataset, selection=selection,
                                                    genMatch='others',
                                                    diphoton_isolation=diphoton.pfIso[diphoton_matched_gen_others])
            else:
                gen_photons = gen_particles[(gen_particles.pdgId==22)]
                gen_electrons = gen_particles[(abs(gen_particles.pdgId)==11)]
                gen_others = gen_particles[(gen_particles.pdgId!=22) & (abs(gen_particles.pdgId)!=11)]

                diphoton_matched_gen_photons = dak.any(diphoton.delta_r(gen_photons) < dr_match, axis=1)
                diphoton_matched_gen_electrons = dak.any(diphoton.delta_r(gen_electrons) < dr_match, axis=1)
                diphoton_matched_gen_others = dak.any(diphoton.delta_r(gen_others) < dr_match, axis=1)

                output['diphoton_genMatch'].fill(dataset=self.dataset, selection=selection,
                                                    genMatch='photon',
                                                    diphoton_isolation=diphoton.pfIso[diphoton_matched_gen_photons])
                
                output['diphoton_genMatch'].fill(dataset=self.dataset, selection=selection,
                                                    genMatch='electron',
                                                    diphoton_isolation=diphoton.pfIso[diphoton_matched_gen_electrons])
                
                output['diphoton_genMatch'].fill(dataset=self.dataset, selection=selection,
                                                    genMatch='others',
                                                    diphoton_isolation=diphoton.pfIso[diphoton_matched_gen_others])

        # Add selections to metadata
        for hname, hist in output.items():
            if isinstance(hist, dh.Hist) or isinstance(hist, bh.Hist):
                hist.metadata = selections

        return output

    def postprocess(self, accumulator):
        pass