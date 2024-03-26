from src.tools import io
from src.tools.models.methylation import OneDimensionalNonCollaborativeMethylation
import matplotlib.pyplot as plt
import os


plt.rcParams['svg.fonttype'] = 'none'

filename = "10k_preprocessed"

# loads data
data = io.read_simulation(filename)
# Extracts basic information
Ms = list(data.keys())
M = 80
sample_count = len(data[M]["data"])
timepoints = list(data[M]["data"][0].keys())
parameters = data[M]["parameters"]
print(parameters)

simulated_average = [0 for _ in timepoints]

for i, time in enumerate(timepoints):
    for sample_data in data[M]["data"]:
        simulated_average[i] = simulated_average[i] + sum(sample_data[time])

for i, x in enumerate(simulated_average):
    simulated_average[i] = x / sample_count

initial_state = data[M]["data"][0][0]

average_model = OneDimensionalNonCollaborativeMethylation(M).get_deterministic_model()

average_data = average_model.generate_simulation_data(parameters, initial_state, timepoints)["data"][0]

calculated_average = [sum(average_data[time]) for time in timepoints]


fig, ax = plt.subplots()

ax.bar([i - 0.2 for i in timepoints], simulated_average, width=0.4, align="center", label = "average of 10k simulations")
ax.bar([i + 0.2 for i in timepoints], calculated_average, width=0.4, align="center", label = "calculated average")

ax.set_xlabel("Timepoint")
ax.set_ylabel(f"Total Cell Count Across All Methylation Levels")

ax.text(0, 0, f"parameters: {parameters}")


ax.legend()
plt.show()

io.save_figure(fig, f"{os.path.basename(__file__)[0:-3]}")