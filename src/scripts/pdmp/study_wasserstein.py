from src.tools.models.pdmp import PDMP
from src.tools.models.rng import generate_poisson
from math import log as log
import matplotlib.pyplot as plt
import random
import numpy as np
from src.constants import BISTABLE_PDMP_PARAMS
folder_path = "output/simulations/study_wasserstein"
#random.seed(0)
import csv


r_hm = 0.1
r_hu = 0.1
r_hm_h = 10
r_hu_u = 2
"""
params = {
    'b': 1,
    'r_uh': r_hm * 2,
    'r_hm': r_hm, 
    'r_mh': r_hu * 2,
    'r_hu': r_hu, 
    'r_uh_h': r_hm_h * 2,
    'r_uh_m': r_hm_h * 4,
    'r_hm_h': r_hm_h,
    'r_hm_m': r_hm_h * 2,
    'r_mh_h': r_hu_u * 2,
    'r_mh_u': r_hu_u * 4,
    'r_hu_h': r_hu_u,
    'r_hu_u': r_hu_u * 2,
}
"""

params = BISTABLE_PDMP_PARAMS


times = [i / 100 for i in range(101)]
replicates = 1000
initial_conditions = [[1, 0, 0], [0, 0, 1]]
model = PDMP()

  
def get_wasserstein_upper_bound(parameters, time):    
    return 8 / (2 ** 0.5 + 2 * time * (parameters["r_hm_h"] - parameters["r_hu_h"])**0.5 * parameters["r_hm"]**0.5)**2



total_time_differences = {time: 0 for time in times}


for _ in range(replicates):
    jump_times = generate_poisson(params["b"], times[-1])
    data = []
    for initial_condition in initial_conditions:
        data.append(model.generate_timepoint_data(params, initial_condition, times, splitting_times=jump_times))
    distances = {}
    for time in times:
        difference = [data[0][time][i] - data[1][time][i] for i in range(3)]
        total_time_differences[time] = total_time_differences[time] + difference[0] - difference[2]
       
  
  
average_time_differences = {}
for time in times:
    average_time_differences[time] = total_time_differences[time] / replicates
bounds = {time: get_wasserstein_upper_bound(params, time) for time in times}
print(average_time_differences)
print(bounds)
multiplicative_tightness = {time: log(bounds[time]) - log(average_time_differences[time]) for time in times}
print(multiplicative_tightness)


series_data = {"bound tightness": multiplicative_tightness}
scale = "linear"
fig, ax = plt.subplots(figsize=(13, 8))

for name, series in series_data.items():
    xs = list(series.keys())
    ys = [series[x] for x in xs]
    if scale == "linear":
        ax.plot(xs, ys, label=name, linewidth=2)
    elif scale == "log":
        ax.semilogy(xs, ys, label=name, linewidth=2)
ax.set_ylabel("log(bound) - log(actual)")
ax.set_xlabel("time")
ax.legend()
plt.show()


