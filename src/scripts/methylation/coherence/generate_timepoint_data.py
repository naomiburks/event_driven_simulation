from src.constants import LIVING_BIRTHRATE_PARAMS
from src.tools.models.methylation import OneDimensionalNonCollaborativeMethylation
from src.tools import io
import numpy as np


Ms = [5, 10, 20, 40, 80]
timepoints = list(range(11))
P = 0.8
simulation_count = 10000
filename = "monte_carlo_timepoints_10K"




def _calculate_initial_condition(P, M):
    condition = [0]*(M + 1)
    condition[int(P * M)] = 10
    return condition
data = []
parameters = LIVING_BIRTHRATE_PARAMS

for M in Ms:
    print(f"starting {simulation_count} simulations for {M}-site model")
    model = OneDimensionalNonCollaborativeMethylation(M)
    initial_condition = _calculate_initial_condition(P, M)
    result = model.generate_simulation_data(parameters, initial_condition, timepoints, sample_count = simulation_count)
    data.append(result)

io.write_simulation(data, filename)
