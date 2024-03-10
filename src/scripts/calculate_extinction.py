from src.tools.models.homogeneous import Birth, Death, Switch, HomogeneousModel
from src.tools.models.methylation import OneDimensionalNonCollaborativeMethylation as OD
from src.constants import LIVING_BIRTHRATE_PARAMS

Ms = [5, 10, 20, 40, 80, 160, 240, 1000, 10000, 100000]
parameters = LIVING_BIRTHRATE_PARAMS

for M in Ms:
    model = OD(M)
    """probabilities = model.calculate_extinction(parameters)
    print(probabilities)"""
    val, vec = model.get_deterministic_model().get_long_term_behavior(parameters)
    print(val)
