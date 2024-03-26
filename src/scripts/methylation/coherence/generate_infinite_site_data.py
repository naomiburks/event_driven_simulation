from src.tools.models.methylation import OneDimensionalNonCollaborativeMethylation
from src.constants import LIVING_BIRTHRATE_PARAMS
from src.tools import io

timepoints = list(range(11))
P = 0.8
simulation_count = 10000
filename = "monte_carlo_limit_timepoints_10K"


initial_condition = {P: 10}

data = []
parameters = LIVING_BIRTHRATE_PARAMS

model = OneDimensionalNonCollaborativeMethylation.get_limit_model()

result = model.generate_simulation_data(parameters, initial_condition, timepoints, sample_count = simulation_count)

io.write_simulation(result, filename)
