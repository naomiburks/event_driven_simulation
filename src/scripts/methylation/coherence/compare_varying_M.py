from src.tools import io
from src.tools.models.methylation import OneDimensionalNonCollaborative
import matplotlib.pyplot as plt
import os
from src.tools import plot

plt.rcParams['svg.fonttype'] = 'none'

filename = "10k_preprocessed"

# loads data
data = io.read_simulation(filename)
# Extracts basic information
Ms = list(data.keys())
sample_count = len(data[Ms[0]]["data"])
timepoints = list(data[Ms[0]]["data"][0].keys())
parameters = data[Ms[0]]["parameters"]

timepoint = timepoints[-1]





simulated_averages = {M: {i / M: 0 for i in range(M + 1)} for M in Ms}

for M in Ms:
    for sample_data in data[M]["data"]:
        for j, datum in enumerate(sample_data[timepoint]):
            simulated_averages[M][j / M] = simulated_averages[M][j / M] + datum

for M in Ms: 
    M_total = sum(simulated_averages[M].values())
    for j in range(M + 1):
        simulated_averages[M][j / M] = simulated_averages[M][j / M] / M_total * M





initial_states = {M: data[M]["data"][0][0] for M in Ms}

average_models = {M: OneDimensionalNonCollaborative(M).get_deterministic_model() for M in Ms}

average_data = {M : average_models[M].generate_simulation_data(parameters, initial_states[M], timepoints)["data"][0] for M in Ms}

calculated_averages = {M: {i / M: average_data[M][timepoint][i] for i in range(M + 1)} for M in Ms}

for M in Ms: 
    M_total = sum(calculated_averages[M].values())
    for j in range(M + 1):
        calculated_averages[M][j / M] = calculated_averages[M][j / M] / M_total * M

for i, data in enumerate([simulated_averages, calculated_averages]):
    fig, ax = plt.subplots()
    if i == 0:
        name = "simulated"
    else:
        name = "calculated"
    for M in Ms:
        series = data[M]
        xs = list(series.keys())
        ys = [series[x] for x in xs]
        ax.plot(xs, ys, label = f"{M} sites", linewidth=2)
    
    ax.set_xlabel("Fraction of sites Methylated")
    ax.set_ylabel(f"Relative Likelihood a living cell will be at this methylation level at time t=10")
    ax.set_ybound(0, 8)
    ax.set_xbound(0, 1)
    ax.text(0, 0, f"parameters: {parameters}")


    ax.legend()
    plt.show()

    io.save_figure(fig, f"{os.path.basename(__file__)[0:-3]}_{name}")





#print(calculated_averages)

#plot.plot_dictionary_series(calculated_averages)

#plot.show()
