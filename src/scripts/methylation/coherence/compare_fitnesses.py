from src.tools import io
from src.tools.models.methylation import OneDimensionalNonCollaborative as ODModel
import matplotlib.pyplot as plt
from src.constants import LIVING_BIRTHRATE_PARAMS
import copy
import os


plt.rcParams['svg.fonttype'] = 'none'


basic_params = {}
for key in ["b_0", "b_M", "d_0", "d_M"]:
    basic_params[key] = 1
for key in ["r_mu", "r_um"]:
    basic_params[key] = 0.1
basic_params["p"] = 1


updated_parameters = [{}, {"b_0": 0, "b_M": 2}, {"b_0": 0, "b_M": 10}, {"d_0": 2.2, "d_M": 0}, {"d_0": 9.5, "d_M": 0}]
names = ["no fitness", "weak birth fitness (0-2)", "strong birth fitness (0-10)", "weak death fitness (2.2-0)", "strong death fitness (9.5-0)"]

timepoint = 10
M = 40
initial_condition = [10 * int(i == M * 4 // 5) for i in range(M + 1)]


parameter_list = []
for info in updated_parameters:
    params = copy.deepcopy(basic_params)
    for key, val in info.items():
        params[key] = val
    parameter_list.append(params)

xs = list(range(M + 1))
y_lists = []
for parameters in parameter_list:
    model = ODModel(M).get_deterministic_model()
    res = model.run(parameters, initial_condition, timepoint)
    total = sum(res)
    normalized_res = [(M + 1) * datum / total for datum in res]
    y_lists.append(normalized_res)




fig, ax = plt.subplots()


for ys, name in zip(y_lists, names):
    ax.plot(xs, ys, label=name)
"""
ax.bar([i - 0.2 for i in range(M + 1)], simulated_average, width=0.4, align="center", label = "average of 10k simulations")
ax.bar([i + 0.2 for i in range(M + 1)], calculated_average, width=0.4, align="center", label = "calculated average")
"""
ax.set_xlabel("Site Count")
ax.set_ylabel(f"Relative Frequence at time t = {timepoint}")
"""
ax.text(0, 0, f"parameters: {parameters}")

"""
ax.legend()

plt.show()

io.save_figure(fig, f"{os.path.basename(__file__)[0:-3]}")
