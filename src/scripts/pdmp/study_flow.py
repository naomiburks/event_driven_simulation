from src.tools.models.pdmp import PDMP
from src.tools.plot import plot_dictionary_series
from src import constants
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import logging

logging.basicConfig(filename="output/scripts_log.txt", level=logging.INFO)
logger = logging.getLogger(__name__)
labels = ["m", "h", "u"]

random.seed(0)
x = 2

params = {
    'b': 1,
    'r_uh': 0.1,
    'r_hm': 0.1, 
    'r_mh': 0.1,
    'r_hu': 0.1, 
    'r_uh_h': 0,
    'r_uh_m': x,
    'r_hm_h': 5.5 * x,
    'r_hm_m': 5.5 * x,
    'r_mh_h': 0,
    'r_mh_u': x,
    'r_hu_h': 0,
    'r_hu_u': x,
}

model = PDMP()

initial_us = [i / 100 for i in range(101)]
duration = 0.4
num_of_timepoints = 8



times = [i * duration / (num_of_timepoints) for i in range(num_of_timepoints + 1)]
results = []
for u in initial_us:
    initial_state = [0, 1 - u, u]
    result = model.generate_timepoint_data(params, initial_state, times, splitting_times=[])
    results.append(result)


unmethylated_times = [i / 100 for i in range(1000)]
u_result = model.generate_timepoint_data(params, [0, 0, 1], unmethylated_times, splitting_times=[])


fig, ax = plt.subplots()

for result in results:
    xs = [result[time][0] for time in times]
    ys = [result[time][2] for time in times]
    if result[0][2] == 1:
        color = (1, 0, 0)
    else: 
        color = (0, 1 - result[0][2], result[0][2])
    """ax.plot(xs, ys, color = color)"""
for time in times:
    print(time)
    xs = [result[time][0] for result in results]
    ys = [result[time][2] for result in results]
    ax.plot(xs, ys, label=f"t = {time}")


xs = [u_result[time][0] for time in unmethylated_times]
ys = [u_result[time][2] for time in unmethylated_times]
ax.plot(xs, ys, label = f"unmethylated trajectory (t < 10)")
plt.legend()
plt.show()
