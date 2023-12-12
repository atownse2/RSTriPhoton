import sys


# import awkward as ak
import dask_awkward as dak

sys.path.append("../")
from analysis.utils import sample_info as s
from analysis.utils import objects as obj

from coffea.nanoevents.methods import nanoaod
from coffea.analysis_tools import PackedSelection


def get_candidates(events):
    # Event selection
    photons = obj.get_photons(events)
    diphotons = obj.get_diphotons(events)

    pho_pass = (photons.pt > 80) & \
                (photons.cutBasedId >= 1)

    dipho_pass = (diphotons.pt > 80) & \
                (diphotons.dipho_score > 0.6) & \
                (diphotons.ruclu_moe > 0)

    photons = photons[pho_pass]
    diphotons = diphotons[dipho_pass]

    ## Create candidates
    candidates = obj.get_triphoton_candidates(photons, diphotons)

    pass_dr = candidates.delta_r > 0.1

    candidates = candidates[pass_dr]

    ## Apply preselection
    pass_preselection = dak.num(candidates) > 0

    events = events[pass_preselection]
    candidates = candidates[pass_preselection] # Already sorted by diphoton score
    candidate = candidates[:,0].compute()

    return candidate






# def apply_preselections(events,
#                         trigger_names = ['HLT_DoublePhoton70_v',
#                                          'HLT_DoublePhoton85_v',
#                                          'HLT_Photon200_v'],
#                         photon_pt_cut = 80,
#                         photon_ID_cut = 1, # At least loose ID
#                         diphoton_pt_cut = 80,
#                         diphoton_score_cut = 0.9,
#                         delta_r_cut = 0.15,
#                         ket_frac_max = 0.5,
#                         ):
    
#     selections = PackedSelection()

#     # Trigger selections
#     if isinstance(trigger_names, str):
#         trigger_names = [trigger_names]

#     for trigger_name in trigger_names:
#         trigger_index = s.trigger_indicies[trigger_name]
#         selections.add(trigger_name, events.triggers[:,trigger_index]==1)
    
#     # Object selection
#     photons = obj.get_photons(events)
#     diphotons = obj.get_diphotons(events)

#     pho_pass = (photons.pt > photon_pt_cut) & \
#                 (photons.cutBasedId >= photon_ID_cut)

#     dipho_pass = (diphotons.pt > diphoton_pt_cut) & \
#                     (diphotons.ruclu_moe > 0)

#     photons = photons[pho_pass]
#     diphotons = diphotons[dipho_pass]

#     selections.add(f'photon_pt{photon_pt_cut}_{photonID_tags[photon_ID_cut]}ID',
#                     dak.any(pho_pass, axis=1))
#     selections.add(f'diphoton_pt{diphoton_pt_cut}',
#                     dak.any(dipho_pass, axis=1))

#     # Candidate selection
#     candidates = obj.get_triphoton_candidates(photons, diphotons)

#     # Delta R preselection cut
#     pass_dr = candidates.delta_r > delta_r_cut
#     pass_ket_frac = candidates.ket_frac < ket_frac_max

#     selections.add(f'candidate_dr{delta_r_cut}', dak.any(pass_dr, axis=1))
#     selections.add(f'candidate_ket_frac{ket_frac_max}', dak.any(pass_ket_frac, axis=1))

#     candidates = candidates[pass_dr & pass_ket_frac]

#     # Pass preselection
#     pass_preselection = selections.all(*selections.names)
    
#     events = events[pass_preselection]
#     events['photons'] = photons[pass_preselection]
#     events['diphotons'] = diphotons[pass_preselection]
#     events['candidates'] = candidates[pass_preselection]
#     # events.behavior = nanoaod.behavior

#     return events, selections



