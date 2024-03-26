import itertools
from math import log
from copy import deepcopy
from src.tools import io
from src.tools.plot import show, plot_dictionary_series

filename = "monte_carlo_timepoints_10k"

# loads data
data = io.read_simulation(filename)
# Extracts basic information
Ms = [len(data[i]["data"][0][0]) - 1 for i in range(len(data))]
sample_counts = len(data[0]["data"])
timepoints = list(data[0]["data"][0].keys())


last_time = timepoints[-1]


def calculate_total_populations(data, timepoint:float, population_index:int):
    population_sizes = []
    population_data = data[population_index]
    for sample in population_data["data"]:
        population_sizes.append(sum(sample[timepoint]))
    return population_sizes

def calculate_average_populations(total_population_sizes):
    return sum(total_population_sizes) / len(total_population_sizes)

def filter_out_extinctions(data):
    new_data = deepcopy(data)
    for i, _ in enumerate(Ms):
        M_data = new_data[i]["data"]
        extinctions = []
        for j, sample in enumerate(M_data):
            if sample[last_time] == 0:
                extinctions.append(j)
        for j in extinctions:
            M_data.pop(j)
    return new_data

filtered_data = filter_out_extinctions(data)


averages = {}

for (i, M), t in itertools.product(enumerate(Ms), timepoints):
    if M not in averages:
        averages[M] = {}
    average = calculate_average_populations(calculate_total_populations(data, t, i))
    averages[M][t] = average

growth = {}



extinction_counts = {}

for i, M in enumerate(Ms):
    sizes = calculate_total_populations(data, last_time, i)
    extinction_count = 0
    for size in sizes:
        if size == 0:
            extinction_count += 1
    extinction_counts[M] = extinction_count

print(averages)
print(extinction_counts)

growth_rates = {}
for M, data in averages.items():
    growth = {a: log(data[a] / data[b]) for a, b in zip(timepoints[1:], timepoints[:-1])}
    growth_rates[M] = growth
fig, ax = plot_dictionary_series(averages, scale="log")
show()
