from src.tools.models.pdmp import PDMP
from src.tools.models import rng
import random
import matplotlib.pyplot as plt
import numpy as np
import os
from src.constants import BISTABLE_STRONG_PDMP_PARAMS
import json

output_path = f"output/simulations/{os.path.basename(__file__)[:-3]}.json"

if input("Resimulate data? y/n: ") == "y":
    
    """random.seed(0)"""

    parameters = BISTABLE_STRONG_PDMP_PARAMS

    


    model = PDMP()
    duration = 10
    frequency = 100
    
    x = 0
    initial_state = [0, 1, 0]
    results = []

    jump_times = [i + 0.01    for i in rng.generate_poisson(parameters["b"], duration)]
    jump_times.append(0.01)
    jump_times = sorted(jump_times)
    
    waiting_times = []
    for j1, j2 in zip(jump_times[1:], jump_times[:-1]):
        waiting_times.append(j1 - j2)
    
    print(waiting_times)
    """jump_times = [j + 0.1 for j in jump_times]
    jump_times.extend([0, 0.1])
    jump_times = sorted(jump_times)"""
    """jump_times = rng.generate_poisson(parameters["b"], duration)"""
    """
    duration = 2000
    jump_times2 = rng.generate_poisson(parameters["b"], duration / 2)
    jump_times_ext = [t + 1025 for t in jump_times]
    jump_times.extend(jump_times_ext)
    """
    times = [i / frequency for i in range(frequency * duration)]
    result = model.generate_timepoint_data(parameters, initial_state, times, splitting_times=jump_times)
    results.append(result)
    json_object = json.dumps(result)

    with open(output_path, "w") as outfile:
        outfile.write(json_object)
else:
    with open(output_path) as infile:
        json_object = infile.read()
    str_key_result = json.loads(json_object)
    result = {}
    for key, val in str_key_result.items():
        result[float(key)] = val
    
fig, ax = plt.subplots()
xs = sorted(list(result.keys()))
ys = [result[x] for x in xs]
ax.plot(xs, ys, label=["m", "h", "u"])
ax.legend()
ax.set_ybound(0, 1)
ax.set_xbound(xs[0], xs[-1])
ax.set_xlabel("Time")
ax.set_ylabel("Fraction of Sites")
plt.show()

