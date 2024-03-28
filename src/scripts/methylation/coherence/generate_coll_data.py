from src.constants import BIRTHRATE_PARAMS_COLL
from src.tools.models.methylation import OneDimensionalColl, OneDimensionalNonCollaborativeMethylation
from src.tools import io
import numpy as np
import copy

Ms = [5]
timepoints = list(range(11))
P = 0.8
simulation_count = 1
filename = "monte_carlo_timepoints_10K_coll"




def _calculate_initial_condition(P, M):
    condition = [0]*(M + 1)
    condition[int(P * M)] = 10
    return condition

data = []
parameters = copy.deepcopy(BIRTHRATE_PARAMS_COLL)
parameters["b_0"] = 1.5
parameters["d_M"] = 0.5

for M in Ms:
    print(f"starting {simulation_count} simulations for {M}-site model")
    model = OneDimensionalColl(M)
    initial_condition = _calculate_initial_condition(P, M)
    result = model.generate_simulation_data(parameters, initial_condition, timepoints, sample_count = simulation_count)
    data.append(result)

print(data)

"""io.write_simulation(data, filename)"""
