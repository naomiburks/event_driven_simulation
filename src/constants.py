BASIC_PARAMS = {
    'b_0': 1.0,
    'b_M': 1.0, 
    'd_0': 1.0,
    'd_M': 1.0, 
    'r_um': 1.0,
    'r_mu': 1.0, 
    'M': 10,
    'p': 1.0,
  }

DYING_PARAMS = {
    'b_0': 2.0,
    'b_M': 2.0, 
    'd_0': 3.0,
    'd_M': 1.5, 
    'r_um': 1.0,
    'r_mu': 1.0, 
    'M': 4,
    'p': 1.0,
  }

LIVING_DEATHRATE_PARAMS = {
    'b_0': 1.8,
    'b_M': 1.8, 
    'd_0': 3,
    'd_M': 1, 
    'r_um': 0.2,
    'r_mu': 0.1, 
    'p': 1,
  }

LIVING_BIRTHRATE_PARAMS = {
    'b_0': 0.8,
    'b_M': 2.8, 
    'd_0': 2,
    'd_M': 2, 
    'r_um': 0.2,
    'r_mu': 0.1, 
    'p': 1,
  }

BIRTHRATE_PARAMS_COLL = {
    'b_0': 1.0,
    'b_M': 1.0, 
    'd_0': 1.0,
    'd_M': 1.0, 
    'r_uh': 0.1,
    'r_hm': 0.5, 
    'r_mh': 0.1,
    'r_hu': 0.1, 
    'r_uh_h': 0,
    'r_uh_m': 10,
    'r_hm_h': 10,
    'r_hm_m': 20,
    'r_mh_h': 0,
    'r_mh_u': 1,
    'r_hu_h': 0,
    'r_hu_u': 10,
    'p': 1.0,
  }

BARELY_LIVING_PARAMS = {
    'b_0': 1.8,
    'b_M': 1.8, 
    'd_0': 3,
    'd_M': 1, 
    'r_um': 0.036,
    'r_mu': 0.036, 
    'M': 100,
    'p': 1,
}

BISTABLE_PDMP_PARAMS = {
    'b': 1,
    'r_uh': 0.35,
    'r_hm': 0.5, 
    'r_mh': 0.1,
    'r_hu': 0.1, 
    'r_uh_h': 5.5,
    'r_uh_m': 11,
    'r_hm_h': 10,
    'r_hm_m': 20,
    'r_mh_h': 5,
    'r_mh_u': 10,
    'r_hu_h': 5,
    'r_hu_u': 10,
}

BISTABLE_STRONG_PDMP_PARAMS = {
    'b': 1,
    'r_uh': 1.4,
    'r_hm': 2, 
    'r_mh': 0.2,
    'r_hu': 1.6, 
    'r_uh_h': 22,
    'r_uh_m': 44,
    'r_hm_h': 40,
    'r_hm_m': 80,
    'r_mh_h': 20,
    'r_mh_u': 40,
    'r_hu_h': 20,
    'r_hu_u': 40,
}

def UNMETHYLATED_INITIAL(M : int, cell_count = 100):
    n = [0] * (M + 1)
    n[0] = cell_count
    return n

def HALF_METHYLATED_INITIAL(M: int, cell_count = 100):
    n = [0] * (M + 1)
    n[M // 2] = cell_count
    return n

def ALL_EQUAL(M: int, cell_count = 100):
    n = [cell_count] * (M + 1)
    return n

CONVERGENCE_TOLERANCE = 0.0000000000001
