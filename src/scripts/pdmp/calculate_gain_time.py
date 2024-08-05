from src.tools.models.pdmp import PDMP
from src.tools.plot import plot_dictionary_series, show
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import logging
from src.constants import BISTABLE_PDMP_PARAMS, BISTABLE_STRONG_PDMP_PARAMS

parameters = BISTABLE_STRONG_PDMP_PARAMS


model = PDMP()
times = [i for i in range(1000)]


starts = [i / 1000 for i in range(1000)]
end = 1.8

hitting = {}
for start in starts:
    hitting_time = model.get_hitting_time(parameters, start, end, stepsize = 0.001)
    hitting[start] = hitting_time
    print(f"hitting time for {start}: {hitting_time}")


fig, ax = plt.subplots()
x = list(hitting.keys())
y = [hitting[key] for key in x]
ax.plot(x, y)
ax.set_xbound(0, 1)
ax.set_ybound(0, 40)
ax.set_xlabel("Fraction of Hemimethylated Sites Post-split")
ax.set_ylabel(f"Time to hit {end} average methylated cytosines per site")
plt.show()

