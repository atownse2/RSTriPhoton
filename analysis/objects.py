import numpy as np
import awkward as ak

from coffea.nanoevents.methods import candidate

radion_pdgid = 9000025
photon_pdgid = 22
bkk_pdgid = 9000121


def get_photons(events):
    return ak.zip(
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
        'charge': ak.zeros_like(events.patpho_energy),
    },
        with_name='PtEtaPhiECandidate',
        behavior=candidate.behavior
    )

def get_diphotons(events):

    if not hasattr(events, 'ruclu_pt'):
        calculate_ruclu_pt(events)

    return ak.zip(
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
            'charge': ak.zeros_like(events.ruclu_energy),
        },
            with_name='PtEtaPhiMCandidate',
            behavior=candidate.behavior
        )

def get_jets(events):
    return ak.zip(
    {
        'energy': events.jet_energy,
        'pt': events.jet_pt,
        'eta': events.jet_eta,
        'phi': events.jet_phi,
        'charge': ak.zeros_like(events.jet_energy),
    },
        with_name='PtEtaPhiECandidate',
        behavior=candidate.behavior
    )

def get_gen_particles(events):
    return ak.zip(
    {
        'energy': events.genpart_energy,
        'mass' : events.genpart_mass,
        'pt': events.genpart_pt,
        'eta': events.genpart_eta,
        'phi': events.genpart_phi,
        'pdgId': events.genpart_pdgid,
        'mother_pdgId': events.genpart_motherpdgid,
        'charge': ak.zeros_like(events.genpart_energy),
    },
        with_name='PtEtaPhiMCandidate',
        behavior=candidate.behavior
    )

def get_triphoton_candidates(photons, diphotons):
    """
    Returns all combinations of photons and diphotons for every event.
    """
    candidates = ak.cartesian( {'photon' : photons,
                          'diphoton': diphotons})
    
    candidates['delta_r'] = candidates.photon.delta_r(candidates.diphoton)
    candidates['delta_eta'] = abs(candidates.photon.eta - candidates.diphoton.eta)
    candidates['triphoton'] = candidates.photon + candidates.diphoton
    candidates['ket_frac'] = candidates.triphoton.pt/candidates.triphoton.energy
    candidates['mass_ratio'] = candidates.diphoton.mass/candidates.triphoton.mass

    # Sort by lowest ket fraction
    candidates = candidates[ak.argsort(candidates.ket_frac, ascending=True)]

    return candidates

def deep_matching(obj1, obj2, dr_max=0.1):
    """Returns a boolean array with same shape as obj1 where 
    it is true if there is some obj2 within dr_max of obj1."""
    matching = ak.cartesian({'obj1':obj1, 'obj2':obj2}, nested=True)
    matching = matching[matching.obj1.delta_r(matching.obj2)<0.1]

    # Sort axis 2 by the percent difference in pt
    matching = matching[ak.argsort(abs(matching.obj1.pt - matching.obj2.pt)/matching.obj1.pt, axis=2, ascending=True)]

    has_match = ak.num(matching.obj2, axis=2) > 0
    obj2 = matching.obj2[has_match,0]

    return has_match, obj2


# Functions for calculating kinematic variables
def calculate_ruclu_pt(events):
    events['ruclu_pt'] = pt(events.ruclu_energy, events.ruclu_eta, events.pvtx_z[:,0])

def ZShift( eta, zPV):
    R = 129
    theta = 2*np.arctan(np.exp(-1*abs(eta)))
    z = R / np.tan(theta) * np.sign(eta)
    zprime = z - zPV
    thetaprime = np.arctan(R/abs(zprime))
    etaprime = -np.sign(zprime)*np.log(np.tan(thetaprime/2))	  
    return etaprime

def pt(E, eta, zpv):
    """Calculate the transverse momentum of a particle given its energy, pseudorapidity, and z position of the primary vertex."""
    pt = E*np.sin(2*np.arctan(np.exp(-1*ZShift(eta,zpv))))
    return pt