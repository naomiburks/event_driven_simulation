from src.constants import LIVING_DEATHRATE_PARAMS
from src.tools.models.methylation import OneDimensionalNonCollaborativeMethylation
from src.tools.plot import plot_extinction_comparison, show
from src.tools.io import save_figure
from copy import copy


infinite_model = OneDimensionalNonCollaborativeMethylation.get_limit_model()
point_count = 1001
base_parameters = copy(LIVING_DEATHRATE_PARAMS)
del base_parameters["r_mu"]
del base_parameters["r_um"]
del base_parameters["M"]
base_parameters["p"] = 1

r_mus = [0.01, 0.1, 0.3, 1, 10]
ratio = 2
results = []
for r_mu in r_mus:
    parameters = copy(base_parameters)
    parameters["r_mu"] = r_mu
    parameters["r_um"] = ratio * r_mu
    extinction_probabilities = infinite_model.calculate_extinction(
        parameters, point_count=point_count)
    extinction_dict = {}
    for i, probability in enumerate(extinction_probabilities):
        extinction_dict[i / (point_count - 1)] = probability
    result = {
    "data": extinction_dict,
    "label": f"r_mu = {r_mu}"
    }
    results.append(result)


info = "M = infinity\n"
for name, param in base_parameters.items():
    info = info + f"{name} = {param}" + "\n"
info = info + f"r_um / r_mu = {ratio}"

results_dict = {
    "results": results,
    "info": info
}
fig, ax = plot_extinction_comparison(results_dict)
show()
save_figure(fig, "infinite_extinction_by_methylation")

