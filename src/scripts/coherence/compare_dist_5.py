from src.tools import io
from src.tools.models.methylation import OneDimensionalNonCollaborativeMethylation
import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'

filename = "10k_preprocessed"

# loads data
data = io.read_simulation(filename)
# Extracts basic information
Ms = list(data.keys())
sample_count = len(data[5]["data"])
timepoints = list(data[5]["data"][0].keys())
parameters = data[5]["parameters"]
print(parameters)
M = Ms[0]
timepoint = timepoints[-1]


simulated_average = [0 for _ in range(M + 1)]

for sample_data in data[M]["data"]:
    timepoint_data = sample_data[timepoint]
    for i, x in enumerate(timepoint_data):
        simulated_average[i] = simulated_average[i] + x
for i, x in enumerate(simulated_average):
    simulated_average[i] = x / sample_count

initial_state = data[5]["data"][0][0]

average_model = OneDimensionalNonCollaborativeMethylation(5).get_deterministic_model()

calculated_average = average_model.generate_simulation_data(parameters, initial_state, [timepoint])["data"][0][timepoint]


fig, ax = plt.subplots()

ax.bar([i - 0.2 for i in range(6)], simulated_average, width=0.4, align="center", label = "average of 10k simulations")
ax.bar([i + 0.2 for i in range(6)], calculated_average, width=0.4, align="center", label = "calculated average")

ax.set_xlabel("Site Count")
ax.set_ylabel(f"Cell Count at Time t={10}")

ax.text(0, 0, f"parameters: {parameters}")


ax.legend()
plt.show()

io.save_figure(fig, "compare_dist_5")
