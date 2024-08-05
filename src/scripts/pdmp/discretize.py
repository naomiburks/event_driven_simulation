from src.tools.models.pdmp import PDMP
from src.constants import BISTABLE_STRONG_PDMP_PARAMS
import numpy as np


def get_hitting_times(hemimethylation_levels, parameters):
    model = PDMP()
    t_mat = []
    for h1 in hemimethylation_levels:
        t_list = []
        for h2 in hemimethylation_levels:
            if h2 <= h1 / 2:
                t_list.append(0)
                continue
            hitting_time = model.get_hitting_time(parameters, h1, h2 * 2)
            t_list.append(hitting_time)
        t_mat.append(t_list)
    return t_mat
    


t_mat = get_hitting_times([i / 5 for i in range(5)], BISTABLE_STRONG_PDMP_PARAMS)
print(np.array(t_mat))




    



