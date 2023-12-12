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


class RuCluProcessor(processor.ProcessorABC):
    def __init__(self):
        pass
    
    def process(self, events):
        dataset = events.metadata['dataset']
        isMC = si.isMC(dataset)
        isSignal = si.isSignal(dataset)

        # Define accumulators
        storage = None #bh.storage.Weight() for when I need to add weights

        output = {}

        output['diphoton_mass'] = bh.Hist(h.dataset_axis, h.diphoton_mass_axis, name='Counts', storage=storage)
        output['diphoton_energy'] = bh.Hist(h.dataset_axis, h.diphoton_energy_axis, name='Counts', storage=storage)
        output['diphoton_pt'] = bh.Hist(h.dataset_axis, h.diphoton_pt_axis, name='Counts', storage=storage)
        output['diphoton_eta'] = bh.Hist(h.dataset_axis, h.diphoton_eta_axis, name='Counts', storage=storage)
        output['diphoton_isolation_score'] = bh.Hist(h.dataset_axis,
                                                       h.diphoton_isolation_axis,
                                                       h.diphoton_score_axis, name='Counts', storage=storage)
        output['diphoton_isPhoton'] = bh.Hist(h.dataset_axis,
                                                bh.axis.Boolean(name='isPhoton', label='dr to pat pho < 0.15'),
                                                name = 'Counts', storage=storage)

        if isMC:
            output['diphoton_genMatch'] = bh.Hist(h.dataset_axis,
                                                    bh.axis.StrCategory([], name='genMatch', growth=True),
                                                    h.diphoton_score_axis,
                                                    h.diphoton_isolation_axis,)
        
        output['cutflow'] = {dataset: {}}
        cutflow = output['cutflow'][dataset]

        cutflow['all events'] = len(events)

        # Get objects
        diphotons = obj.get_diphotons(events)

        diphotons = diphotons[(diphotons.pt>80) &
                              (diphotons.ruclu_moe>0)]

        output['diphoton_mass'].fill(dataset=dataset, diphoton_mass=ak.flatten(diphotons.mass))
        output['diphoton_energy'].fill(dataset=dataset, diphoton_energy=ak.flatten(diphotons.energy))
        output['diphoton_pt'].fill(dataset=dataset, diphoton_pt=ak.flatten(diphotons.pt))
        output['diphoton_eta'].fill(dataset=dataset, diphoton_eta=ak.flatten(diphotons.eta))
        output['diphoton_isolation_score'].fill(dataset=dataset,
                                                diphoton_isolation=ak.flatten(diphotons.pfIso),
                                                diphoton_score=ak.flatten(diphotons.dipho_score))
        output['diphoton_isPhoton'].fill(dataset=dataset,
                                            isPhoton=ak.flatten(diphotons.isPhoton))
        


        if isMC:
            gen_particles = obj.get_gen_particles(events)
            gen_particles = gen_particles[gen_particles.pt>70]

            # Matching
            has_match, gen_matches = obj.deep_matching(diphotons, gen_particles, dr_max=0.1)
            diphotons = diphotons[has_match]

            if isSignal:
                is_radion = (gen_matches.pdgId==obj.radion_pdgid) & (gen_matches.mother_pdgId==obj.bkk_pdgid)
                is_prompt_photon = (gen_matches.pdgId==22) & (gen_matches.mother_pdgId==obj.bkk_pdgid)
                is_others = (gen_matches.pdgId!=obj.radion_pdgid) & (gen_matches.pdgId!=22)

                output['diphoton_genMatch'].fill(dataset=dataset,
                                                genMatch='radion',
                                                diphoton_score=ak.flatten(diphotons.dipho_score[is_radion]),
                                                diphoton_isolation=ak.flatten(diphotons.pfIso[is_radion]))
                output['diphoton_genMatch'].fill(dataset=dataset,
                                                    genMatch='prompt_photon',
                                                    diphoton_score=ak.flatten(diphotons.dipho_score[is_prompt_photon]),
                                                    diphoton_isolation=ak.flatten(diphotons.pfIso[is_prompt_photon]))
                output['diphoton_genMatch'].fill(dataset=dataset,
                                                    genMatch='others',
                                                    diphoton_score=ak.flatten(diphotons.dipho_score[is_others]),
                                                    diphoton_isolation=ak.flatten(diphotons.pfIso[is_others]))
            else:
                is_photon = (gen_matches.pdgId==22)
                is_electron = (abs(gen_matches.pdgId)==11)
                is_others = (gen_matches.pdgId!=22) & (abs(gen_matches.pdgId)!=11)

                output['diphoton_genMatch'].fill(dataset=dataset,
                                                    genMatch='photon',
                                                    diphoton_score=ak.flatten(diphotons.dipho_score[is_photon]),
                                                    diphoton_isolation=ak.flatten(diphotons.pfIso[is_photon]))
                
                output['diphoton_genMatch'].fill(dataset=dataset,
                                                    genMatch='electron',
                                                    diphoton_score=ak.flatten(diphotons.dipho_score[is_electron]),
                                                    diphoton_isolation=ak.flatten(diphotons.pfIso[is_electron]))
                
                output['diphoton_genMatch'].fill(dataset=dataset,
                                                    genMatch='others',
                                                    diphoton_score=ak.flatten(diphotons.dipho_score[is_others]),
                                                    diphoton_isolation=ak.flatten(diphotons.pfIso[is_others]))
        return output

    def postprocess(self, accumulator):
        pass