from src.constants import LIVING_DEATHRATE_PARAMS
from src.tools.models.methylation import OneDimensionalNonCollaborative
from src.tools.plot import plot_extinction_comparison, show
from src.tools.io import save_figure
from copy import copy

finite_models = []
Ms = [3, 10, 30, 100]
for M in Ms:
    finite_models.append(OneDimensionalNonCollaborative(M))

infinite_model = finite_models[0].get_limit_model()
point_count = 1001
base_parameters = copy(LIVING_DEATHRATE_PARAMS)
base_parameters["p"] = 1

results = []
for M, model in zip(Ms, finite_models):
    parameters = copy(base_parameters)
    extinction_probabilities = model.calculate_extinction(parameters)
    extinction_dict = {}
    for i, probability in enumerate(extinction_probabilities):
        extinction_dict[i / M] = probability
    result = {
    "data": extinction_dict,
    "label": f"M = {M}"
    }
    results.append(result)

infinite_extinction_probabilities = infinite_model.calculate_extinction(parameters, point_count=point_count)
extinction_dict = {}
for i, probability in enumerate(infinite_extinction_probabilities):
    extinction_dict[i / (point_count + 1)] = probability
result = {
"data": extinction_dict,
"label": "M = infinite"
}
results.append(result)


info = ""
for name, param in base_parameters.items():
    info = info + f"{name} = {param}" + "\n"

results_dict = {
    "results": results,
    "info": info
}

#print(results_dict)

fig, ax = plot_extinction_comparison(results_dict)
show()
save_figure(fig, "extinction_by_site_count")
