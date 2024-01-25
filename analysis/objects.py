import numpy as np

import dask_awkward as dak

from coffea.nanoevents.methods import candidate

radion_pdgid = 9000025
photon_pdgid = 22
bkk_pdgid = 9000121

photon_mva_cut = {'barrel': {80: 0.42, 90: -0.02},
                  'endcap': {80: 0.14, 90: -0.26}}

photonID_tags = {1: 'loose',
                 2: 'medium',
                 3: 'tight',
                 0: 'no'}

def get_photons(events):
    return dak.zip(
    {
        'energy': events.patpho_energy,
        'pt': events.patpho_pt,
        'eta': events.patpho_eta,
        'phi': events.patpho_phi,
        #'isEB' : events.patpho_eta < 1.4442,
        #'mvaID': events.patpho_mvaID,
        'cutBasedId': events.patpho_cutBased,
        #'hasPixelSeed' : events.patpho_hasPixelSeed,
        #'passElectronVeto' : events.patpho_passElectronVeto,
        'charge': dak.zeros_like(events.patpho_energy),
    },
        with_name='PtEtaPhiECandidate',
        behavior=candidate.behavior
    )

def get_diphotons(events):

    if not hasattr(events, 'ruclu_pt'):
        calculate_ruclu_pt(events)

    return dak.zip(
        {
            'pt': events.ruclu_pt,
            'eta': events.ruclu_eta,
            'phi': events.ruclu_phi,
            'mass' : events.ruclu_moe*events.ruclu_energy,
            'energy': events.ruclu_energy,
            'ruclu_moe' : events.ruclu_moe,
            'pfIso' : events.ruclu_pfIso,
            'dipho_score': events.ruclu_dipho,
            'monopho_score': events.ruclu_monopho,
            'hadron_score': events.ruclu_hadron,
            'isPhoton' : events.ruclu_isPhoton,
            'hasPixelSeed' : events.ruclu_hasPixelSeed,
            'passElectronVeto' : events.ruclu_passElectronVeto,
            'charge': dak.zeros_like(events.ruclu_energy),
        },
            with_name='PtEtaPhiMCandidate',
            behavior=candidate.behavior
        )

def get_jets(events):
    return dak.zip(
    {
        'energy': events.jet_energy,
        'pt': events.jet_pt,
        'eta': events.jet_eta,
        'phi': events.jet_phi,
        'charge': dak.zeros_like(events.jet_energy),
    },
        with_name='PtEtaPhiECandidate',
        behavior=candidate.behavior
    )

def get_gen_particles(events):
    return dak.zip(
    {
        'energy': events.genpart_energy,
        'mass' : events.genpart_mass,
        'pt': events.genpart_pt,
        'eta': events.genpart_eta,
        'phi': events.genpart_phi,
        'pdgId': events.genpart_pdgid,
        'mother_pdgId': events.genpart_motherpdgid,
        'charge': dak.zeros_like(events.genpart_energy),
    },
        with_name='PtEtaPhiMCandidate',
        behavior=candidate.behavior
    )

def get_triphoton_candidates(photons,diphotons):
    """
    Returns all combinations of photons and diphotons for every event.
    """
    combinations = dak.cartesian( [photons, diphotons])

    candidates = dak.zip({'photon': combinations['0'],
                          'diphoton': combinations['1']})

    candidates['triphoton'] = candidates.photon + candidates.diphoton    
    candidates['delta_r'] = candidates.photon.delta_r(candidates.diphoton)
    candidates['delta_eta'] = abs(candidates.photon.eta - candidates.diphoton.eta)
    candidates['ket_frac'] = candidates.triphoton.pt/candidates.triphoton.energy
    candidates['alpha'] = candidates.diphoton.mass/candidates.triphoton.mass

    # Sort by diphoton score
    candidates = candidates[dak.argsort(candidates.diphoton.dipho_score, axis=1, ascending=False)]

    return candidates

def deep_matching(obj0, obj1, dr_max=0.1):
    """Returns a boolean array with same shape as obj0 where 
    it is true if there is some obj1 within dr_max of obj0. Also
    returns the array of obj1, sorted by the percent difference in pt,
    that matches the shape of obj0[has_match]."""

    matching = dak.cartesian([obj0, obj1], nested=True)
    matching = matching[matching['0'].delta_r(matching['1'])<0.1]
    has_match = dak.num(matching['1'], axis=2) > 0

    # Sort axis 2 by the percent difference in pt
    matching = matching[dak.argsort(abs(matching['0'].pt - matching['1'].pt)/matching['0'].pt, axis=2, ascending=True)]
    matched_obj1 = matching['1'][has_match,0]

    return has_match, matched_obj1

# Functions for calculating kinematic variables
def calculate_ruclu_pt(events):
    E = events.ruclu_energy
    eta = events.ruclu_eta
    zpv = events.pvtx_z[:,0]

    events['ruclu_pt'] = E*np.sin(2*np.arctan(np.exp(-1*ZShift(eta,zpv))))

def ZShift( eta, zPV):
    R = 129
    theta = 2*np.arctan(np.exp(-1*abs(eta)))
    z = R / np.tan(theta) * np.sign(eta)
    zprime = z - zPV
    thetaprime = np.arctan(R/abs(zprime))
    etaprime = -np.sign(zprime)*np.log(np.tan(thetaprime/2))	  
    return etaprime