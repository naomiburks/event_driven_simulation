from src.tools.models.pdmp import PDMP
from src.tools.plot import plot_dictionary_series, show
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import logging
from src.constants import BISTABLE_PDMP_PARAMS, BISTABLE_STRONG_PDMP_PARAMS
random.seed(0)

parameters = BISTABLE_STRONG_PDMP_PARAMS

"""
parameters = BISTABLE_PDMP_PARAMS
"""

model = PDMP()
times = [i for i in range(1000)]
splitting_times = [i for i in range(200)]
x = 0



thresholds = [i / 1000 for i in range(999)]

recovery = {}
for threshold in thresholds:
    recovery_time = model.get_hitting_time(parameters, threshold, threshold * 2, stepsize = 0.001)
    recovery[threshold] = recovery_time
    print(f"recovery time for {threshold}: {recovery_time}")


fig, ax = plt.subplots()
x = list(recovery.keys())
y = [recovery[key] for key in x]
ax.plot(x, y)
ax.set_xbound(0, 1)
ax.set_ybound(0, 25)
ax.set_xlabel("Fraction of Hemimethylated Sites Post-split")
ax.set_ylabel("Pre-split Methylation Recovery Time")
plt.show()

#result = model.generate_timepoint_data(parameters, [0, x, 1 - x], times, splitting_times=[])
#result = model.generate_timepoint_data(parameters, [0, 0, 1], times)
#result = model.generate_timepoint_data(parameters, [1, 0, 0], times)

fig, ax = plot_dictionary_series({"PDMP" : result})
show()

#result = model.sample_wasserstein(parameters, times, 100)
#print(result)