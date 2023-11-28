import sys
import awkward as ak
# import dask_awkward as dak


from coffea import processor
from coffea.analysis_tools import PackedSelection, Weights

import hist as bh
# import dask_histogram as dh

sys.path.append("../../")
import analysis.utils.sample_info as si
import analysis.objects as obj
import analysis.utils.histograms as hg
import analysis.selections as sel

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
        selections = {}

        # Define Candidates
        
        # Select Events
        Event = sel.EventSelection(events) # Default cuts for preselection
        Event.apply_preselections()

        preselections = Event.selections #PackedSelection object
        events = Event.events

        cutflow['all events'] = len(events)
        for selection in preselections.names:
            cutflow[selection] = ak.sum(preselections.all(selection))
        cutflow['preselection'] = ak.sum(preselections.all(*preselections.names))

        selection = 'preselection'
        selections[selection] = ", ".join(preselections.names)

        candidate = events.candidate
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
        output['diphoton_mass'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_mass_axis, name='Counts', storage=storage)
        output['diphoton_energy'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_energy_axis, name='Counts', storage=storage)
        output['diphoton_pt'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_pt_axis, name='Counts', storage=storage)
        output['diphoton_eta'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_eta_axis, name='Counts', storage=storage)
        output['diphoton_score'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_score_axis, name='Counts', storage=storage)
        output['diphoton_isolation'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.diphoton_isolation_axis, name='Counts', storage=storage)

        output['photon_energy'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.photon_energy_axis, name='Counts', storage=storage)
        output['photon_pt'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.photon_pt_axis, name='Counts', storage=storage)
        output['photon_eta'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.photon_eta_axis, name='Counts', storage=storage)

        output['ket_frac'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.ket_frac_axis, name='Counts', storage=storage)
        output['energy_ratio'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.energy_ratio_axis, name='Counts', storage=storage)
        output['alpha'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.alpha_axis, name='Counts', storage=storage)
        output['delta_r'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.delta_r_axis, name='Counts', storage=storage)
        output['delta_eta'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.delta_eta_axis, name='Counts', storage=storage)

        output['triphoton_mass'] = bh.Hist(hg.dataset_axis, hg.selection_axis,
                                            hg.get_axis('triphoton_mass', 2048),
                                            name='Counts', storage=storage)
        output['diphoton_mass'] = bh.Hist(hg.dataset_axis, hg.selection_axis,
                                          hg.get_axis('diphoton_mass', 1024),
                                          name='Counts', storage=storage)

        # Reduce binning so that memory doesn't explode
        output['triphoton_diphoton_mass'] = bh.Hist(hg.dataset_axis, hg.selection_axis,
                                                    hg.get_axis('triphoton_mass', 512),
                                                    hg.get_axis('diphoton_mass', 256),
                                                    name='Counts', storage=storage)

        output['pf_met'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.pf_met_axis, name='Counts', storage=storage)
        output['jet_energy_frac'] = bh.Hist(hg.dataset_axis, hg.selection_axis, hg.jet_energy_frac_axis, name='Counts', storage=storage)

        if isMC:
            output['diphoton_genMatch'] = bh.Hist(hg.dataset_axis,
                                                  hg.selection_axis, 
                                                  bh.axis.StrCategory([], name='genMatch', growth=True),
                                                  hg.diphoton_isolation_axis,)

        ## Filling histograms
        output['diphoton_mass'].fill(dataset=dataset, selection=selection, diphoton_mass=diphoton.mass)
        output['diphoton_energy'].fill(dataset=dataset, selection=selection, diphoton_energy=diphoton.energy)
        output['diphoton_pt'].fill(dataset=dataset, selection=selection, diphoton_pt=diphoton.pt)
        output['diphoton_eta'].fill(dataset=dataset, selection=selection, diphoton_eta=diphoton.eta)
        output['diphoton_score'].fill(dataset=dataset, selection=selection, diphoton_score=diphoton.dipho_score)
        output['diphoton_isolation'].fill(dataset=dataset, selection=selection, diphoton_isolation=diphoton.pfIso)

        output['photon_energy'].fill(dataset=dataset, selection=selection, photon_energy=photon.energy)
        output['photon_pt'].fill(dataset=dataset, selection=selection, photon_pt=photon.pt)
        output['photon_eta'].fill(dataset=dataset, selection=selection, photon_eta=photon.eta)

        output['ket_frac'].fill(dataset=dataset, selection=selection, ket_frac=candidate.ket_frac)
        output['energy_ratio'].fill(dataset=dataset, selection=selection, energy_ratio=diphoton.energy/photon.energy)
        output['alpha'].fill(dataset=dataset, selection=selection, alpha=diphoton.mass/triphoton.mass)
        output['delta_r'].fill(dataset=dataset, selection=selection, delta_r=candidate.delta_r)
        output['delta_eta'].fill(dataset=dataset, selection=selection, delta_eta=candidate.delta_eta)

        output['triphoton_mass'].fill(dataset=dataset, selection=selection, triphoton_mass=triphoton.mass)
        output['triphoton_diphoton_mass'].fill(dataset=dataset, selection=selection, triphoton_mass=triphoton.mass, diphoton_mass=diphoton.mass)

        output['pf_met'].fill(dataset=dataset, selection=selection, pf_met=events.met_energy[:,0])
        output['jet_energy_frac'].fill(dataset=dataset, selection=selection, jet_energy_frac=jet_energy_frac)

        # Apply loose diphoton score cut
        dipho_score_min = 0.9
        events = events[events.candidate.diphoton.dipho_score > dipho_score_min]

        selection = 'preselection + diphoton score'
        selections[selection] = ", ".join(preselections.names) + f', diphoton_score{dipho_score_min}'

        candidate = events.candidate
        photon = candidate.photon
        diphoton = candidate.diphoton
        triphoton = candidate.triphoton

        output['diphoton_mass'].fill(dataset=dataset, selection=selection, diphoton_mass=candidate.diphoton.mass)
        output['triphoton_mass'].fill(dataset=dataset, selection=selection, triphoton_mass=candidate.triphoton.mass)

        output['triphoton_diphoton_mass'].fill(dataset=dataset, selection=selection, triphoton_mass=candidate.triphoton.mass, diphoton_mass=candidate.diphoton.mass)

        if isMC:
            gen_particles = obj.get_gen_particles(events)
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

                output['diphoton_genMatch'].fill(dataset=dataset, selection=selection,
                                                genMatch='radion',
                                                diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_radions])
                output['diphoton_genMatch'].fill(dataset=dataset, selection=selection,
                                                    genMatch='prompt_photon',
                                                    diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_prompt_photons])
                output['diphoton_genMatch'].fill(dataset=dataset, selection=selection,
                                                    genMatch='others',
                                                    diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_others])
            else:
                gen_photons = gen_particles[(gen_particles.pdgId==22)]
                gen_electrons = gen_particles[(abs(gen_particles.pdgId)==11)]
                gen_others = gen_particles[(gen_particles.pdgId!=22) & (abs(gen_particles.pdgId)!=11)]

                diphoton_matched_gen_photons = ak.any(candidate.diphoton.delta_r(gen_photons) < dr_match, axis=1)
                diphoton_matched_gen_electrons = ak.any(candidate.diphoton.delta_r(gen_electrons) < dr_match, axis=1)
                diphoton_matched_gen_others = ak.any(candidate.diphoton.delta_r(gen_others) < dr_match, axis=1)

                output['diphoton_genMatch'].fill(dataset=dataset, selection=selection,
                                                    genMatch='photon',
                                                    diphoton_mass=candidate.diphoton.mass[diphoton_matched_gen_photons],
                                                    diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_photons])
                
                output['diphoton_genMatch'].fill(dataset=dataset, selection=selection,
                                                    genMatch='electron',
                                                    diphoton_mass=candidate.diphoton.mass[diphoton_matched_gen_electrons],
                                                    diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_electrons])
                
                output['diphoton_genMatch'].fill(dataset=dataset, selection=selection,
                                                    genMatch='others',
                                                    diphoton_mass=candidate.diphoton.mass[diphoton_matched_gen_others],
                                                    diphoton_isolation=candidate.diphoton.pfIso[diphoton_matched_gen_others])

        # Add selections to metadata
        for hname, hist in output.items():
            if isinstance(hist, bh.Hist):
                hist.metadata = selections

        return output

    def postprocess(self, accumulator):
        pass