from src.tools.models.methylation import OneDimensionalNonCollaborative as OD, OneDimensionalCollaborative as OD_C
from random import random, randint
from src.constants import CONVERGENCE_TOLERANCE
from copy import deepcopy


def test_model_coherence():
    M = 10
    params = {
        'b_0': random(),
        'b_M': random(), 
        'd_0': random(),
        'd_M': random(), 
        'r_uh': random(),
        'r_hm': random(), 
        'r_mh': random(),
        'r_hu': random(), 
        'r_uh_h': 0,
        'r_uh_m': 0,
        'r_hm_h': 0,
        'r_hm_m': 0,
        'r_mh_h': 0,
        'r_mh_u': 0,
        'r_hu_h': 0,
        'r_hu_u': 0,
        'p': 1.0,
    }
    
    params["r_mu"] = params["r_mh"] * params["r_hu"] / (params["r_hu"] + params["r_hm"])
    params["r_um"] = params["r_uh"] * params["r_hm"] / (params["r_hu"] + params["r_hm"])
    coll = OD_C(M)
    non_coll = OD(M)
    coll_ext = coll.calculate_extinction(params)
    non_coll_ext = non_coll.calculate_extinction(params)

    for a, b in zip(coll_ext, non_coll_ext):
        assert(- CONVERGENCE_TOLERANCE < a - b)
        assert(a - b < CONVERGENCE_TOLERANCE)    


def test_model_coherence2():
    M = 10
    base_params = {
        'b_0': 0,
        'b_M': 2, 
        'd_0': 0.5,
        'd_M': 0.5, 
        'r_uh': 0.1,
        'r_hm': 0.1, 
        'r_mh': 0.1,
        'r_hu': 0.1, 
        'r_uh_h': 0.01,
        'r_uh_m': 0.01,
        'r_hm_h': 0.01,
        'r_hm_m': 0.01,
        'r_mh_h': 0.01,
        'r_mh_u': 0.01,
        'r_hu_h': 0.01,
        'r_hu_u': 0.01,
        'p': 1.0,
    }
    
    base_params["r_mu"] = base_params["r_mh"] * base_params["r_hu"] / (base_params["r_hu"] + base_params["r_hm"])
    base_params["r_um"] = base_params["r_uh"] * base_params["r_hm"] / (base_params["r_hu"] + base_params["r_hm"])

    coll = OD_C(M)
    non_coll = OD(M)
    

    
    coll_ext = coll.calculate_extinction(base_params)
    non_coll_ext = non_coll.calculate_extinction(base_params)

    assert(coll_ext[0] > non_coll_ext[0])
    assert(coll_ext[-1] < non_coll_ext[-1])


test_model_coherence()