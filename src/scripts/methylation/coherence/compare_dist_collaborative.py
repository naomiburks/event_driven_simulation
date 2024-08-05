from src.constants import BIRTHRATE_PARAMS_COLL
from src.tools.models.methylation import OneDimensionalCollaborative, OneDimensionalNonCollaborative
from src.tools import io
from src.tools import plot
import numpy as np
import copy

Ms = [5, 10, 20, 50, 100, 200, 500]
timepoints = list(range(11))
P = 0.5


def _calculate_initial_condition(P, M):
    condition = [0]*(M + 1)
    condition[int(P * M)] = 10
    return condition

data = []
parameters = copy.deepcopy(BIRTHRATE_PARAMS_COLL)
parameters["b_0"] = 1
parameters["b_M"] = 1
parameters["d_M"] = 0.8
parameters["d_0"] = 0.8
parameters["r_hm"] = 0.1
parameters["r_mh"] = 0.1
parameters["r_uh"] = 0.1
parameters["r_hu"] = 0.1

parameters["r_uh_m"] = 0
parameters["r_mh_u"] = 0

parameters["r_hm_m"] = 0
parameters["r_hu_u"] = 0

#parameters["r_uh_h"] = 0.00
#parameters["r_mh_h"] = 0.00

parameters["r_hm_h"] = 0
parameters["r_hu_h"] = 0

data = {}
for M in Ms: 
    model = OneDimensionalCollaborative(M)
    growth, dist = model.get_deterministic_model().get_long_term_behavior(parameters)
    model_data = {}
    for i in range(M + 1):
        model_data[i / M] = dist[i] * (M + 1)
    data[M] = model_data

plot.plot_dictionary_series(data, scale="linear")
plot.show()
