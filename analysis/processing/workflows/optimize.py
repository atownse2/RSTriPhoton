import sys

import dask_awkward as dak

import hist.dask as dh
import hist as bh

from coffea import processor

top_dir = "/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/"
sys.path.append(top_dir)

import analysis.utils.objects as obj
import analysis.utils.histograms as hgm
import analysis.utils.sample_info as si
import analysis.utils.fitting as fit

photon_mva_cut = {'barrel': {80: 0.42, 90: -0.02},
                  'endcap': {80: 0.14, 90: -0.26}}

# Exploring a grid of cuts on:
# - Photon ID
# - Diphoton Score
# - Diphoton Isolation
# - Delta R
# - KET fraction
# - Hadronic energy fraction

photon_id_axis = bh.axis.IntCategory([], name='photon_id', label='Photon ID', growth=True, flow=False)
diphoton_score_axis = bh.axis.Regular(4, 0.6, 1, name='diphoton_score', label='Diphoton Score', underflow=False, overflow=False)
diphoton_isolation_axis = bh.axis.Regular(4, 0.6, 1, name='diphoton_isolation', label='Diphoton Isolation', underflow=False, overflow=False)
delta_r_axis = bh.axis.Regular(2, 2, 3, name='delta_r', label='$\Delta R$', underflow=False)
ket_frac_axis = bh.axis.Regular(3, 0.1, 0.4, name='ket_frac', label='KET Fraction', overflow=False)
hadron_frac_axis = bh.axis.Regular(3, 0.1, 0.4, name='hadron_frac', label='Hadron Energy Fraction', overflow=False)

class OptimizationProcessor(processor.ProcessorABC):
    def __init__(self,
                trigger_names = ['HLT_DoublePhoton70_v',
                                'HLT_DoublePhoton85_v',
                                'HLT_Photon200_v'],
                photon_pt_cut = 80,
                photon_ID_cut = 1, # At least loose ID
                diphoton_pt_cut = 80,
                delta_r_cut = 0.15,
                mass_grid = si.get_mass_grid('old'), # Fix later
                ):
        
        self.mass_grid = mass_grid
        
        self._trigger_names = trigger_names
        if isinstance(trigger_names, str):
            self.trigger_names = [self.trigger_names]

        self.photon_pt_cut = photon_pt_cut
        self.photon_ID_cut = photon_ID_cut
        self.photon_ID = obj.photonID_tags[self.photon_ID_cut]

        self.diphoton_pt_cut = diphoton_pt_cut

        self.delta_r_cut = delta_r_cut
    
    def make_output(self):
        storage = None #dh.storage.Weight() for when I need to add weights

        output = {}

        # Cutflow
        output['cutflow'] = {self.dataset: {}}

        # Histograms
        output['in_signal_region'] = dh.Hist(bh.axis.StrCategory([], name='signal_point', label='Signal Point', growth=True, flow=False),
                                             photon_id_axis,
                                             diphoton_score_axis,
                                             diphoton_isolation_axis,
                                             delta_r_axis,
                                             ket_frac_axis,
                                             hadron_frac_axis,)

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

        ## Get quantities
        photon_id = photon.cutBasedId
        diphoton_score = diphoton.dipho_score
        diphoton_isolation = diphoton.ruclu_iso
        delta_r = candidate.delta_r
        ket_frac = candidate.ket_frac

        ### Percent hadronic energy
        jets = obj.get_jets(events)

        jet_not_photon = jets.delta_r(photon) > 0.4
        jet_not_diphoton = jets.delta_r(diphoton) > 0.4

        jets = jets[jet_not_photon & jet_not_diphoton & (jets.pt>30)]

        sum_jet_energy = dak.sum(jets.energy, axis=1)
        jet_energy_frac = sum_jet_energy/(sum_jet_energy+photon.energy+diphoton.energy)

        # Fill histograms
        for signal_point in self.mass_grid:
            signal_filetag = si.get_filetag(signal_point)
            in_ellipse = fit.in_ellipse(signal_point, triphoton.mass, diphoton.mass)
            output['in_signal_region'].fill(signal_point=signal_filetag,
                                            photon_id=photon_id[in_ellipse],
                                            diphoton_score=diphoton_score[in_ellipse],
                                            diphoton_isolation=diphoton_isolation[in_ellipse],
                                            delta_r=delta_r[in_ellipse],
                                            ket_frac=ket_frac[in_ellipse],
                                            hadron_frac=jet_energy_frac[in_ellipse],
                                            )

        return output
    
    def postprocess(self, accumulator):
        pass