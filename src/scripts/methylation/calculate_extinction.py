from src.tools.models.methylation import OneDimensionalNonCollaborativeMethylation as OD, BranchingDiffusion as BD
from src.constants import LIVING_BIRTHRATE_PARAMS

Ms = [5, 10, 20, 40, 80]
parameters = LIVING_BIRTHRATE_PARAMS
"""
for M in Ms:
    model = OD(M)
    probabilities = model.calculate_extinction(parameters)
    print(probabilities)
    
    val, vec = model.get_deterministic_model().get_long_term_behavior(parameters)
    print(val)
"""


def r_b(x, parameters):
    return parameters["b_0"] * (1 - x) + parameters["b_M"] * x

def r_d(x, parameters):
    return parameters["d_0"] * (1 - x) + parameters["d_M"] * x

def diffusion(x, parameters):
    return parameters["r_um"] * (1 - x) - parameters["r_mu"] * x






model = BD(r_b, r_d, diffusion)


extinction = model.calculate_extinction(parameters, point_count = 7)
print(extinction)


