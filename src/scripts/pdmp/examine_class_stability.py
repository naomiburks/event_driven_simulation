from src.tools.models.pdmp import PDMP
from src.tools.plot import plot_dictionary_series
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

initial_state = [0.4, 0.2, 0.4]
times = [i / 10 for i in range(50001)]
result = model.generate_timepoint_data(params, initial_state, times)

fig, ax = plt.subplots()

xs = list(result.keys())
for i in range(2 + 1):
    ys = []
    for x in xs:
        ys.append(result[x][i])
    ax.plot(xs, ys, label=labels[i])
plt.legend()
plt.show()

