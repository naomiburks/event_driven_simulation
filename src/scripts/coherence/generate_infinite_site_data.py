from src.tools.models.methylation import BranchingDiffusion2
from src.constants import LIVING_BIRTHRATE_PARAMS

diff = lambda x, parameters: parameters["r_um"] * (1 - x) - parameters["r_mu"] * x
r_b = lambda x, parameters: parameters["b_0"] * (1 - x) + parameters["b_M"] * x
r_d = lambda x, parameters: parameters["d_0"] * (1 - x) + parameters["d_M"] * x


initial_condition = {1: 10}

times = list(range(2))
model = BranchingDiffusion2(r_b, r_d, diff)

print(r_b(0, LIVING_BIRTHRATE_PARAMS))

result = model.run(LIVING_BIRTHRATE_PARAMS, initial_condition, 1)
print(result)



extinction = model.calculate_extinction(LIVING_BIRTHRATE_PARAMS, point_count = 10)
print(extinction)
