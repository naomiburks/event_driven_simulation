from src.constants import BIRTHRATE_PARAMS_COLL
from src.tools.models.methylation import OneDimensionalCollaborative, OneDimensionalNonCollaborative
from src.tools import io
import numpy as np
import copy

Ms = [5, 10, 20, 40]
timepoints = list(range(11))
P = 0.5
simulation_count = 1
filename = f"monte_carlo_timepoints_{simulation_count}_coll"




def _calculate_initial_condition(P, M):
    condition = [0]*(M + 1)
    condition[int(P * M)] = 10
    return condition

data = []
parameters = copy.deepcopy(BIRTHRATE_PARAMS_COLL)
parameters["b_0"] = 1
parameters["b_M"] = 1
parameters["d_M"] = 0.6
parameters["d_0"] = 0.6
parameters["r_hm"] = 0.1
parameters["r_uh_m"] = 0.02
parameters["r_hm_m"] = 0.02
parameters["r_uh_h"] = 0.01
parameters["r_hm_h"] = 0.01
parameters["r_hu_u"] = 0.02
parameters["r_hu_h"] = 0.01
parameters["r_mh_u"] = 0.02
parameters["r_mh_h"] = 0.01


for M in Ms:
    print(f"starting {simulation_count} simulations for {M}-site model")
    model = OneDimensionalCollaborative(M)
    initial_condition = _calculate_initial_condition(P, M)
    result = model.generate_simulation_data(parameters, initial_condition, timepoints, sample_count = simulation_count)
    data.append(result)

print(data)

io.write_simulation(data, filename)
