from src.tools.models.pdmp import PDMP
from src.tools.plot import show
import random
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit

from src.constants import BISTABLE_STRONG_PDMP_PARAMS
import json

output_path = f"output/simulations/{os.path.basename(__file__)[:-3]}.json"


if input("Resimulate data? WARNING: may take hours. y/n: ") == "y":
    
    random.seed(0)

    parameters = BISTABLE_STRONG_PDMP_PARAMS

    


    model = PDMP()
    times = [i / 10 for i in range(10000)]
    x = 0
    samples = 1000
    result = model.sample_wasserstein(parameters, times, samples)
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



def curve(x, a, b):
    return a * np.exp(b * x)


xs = sorted(list(result.keys()))


ys = [result[x] for x in xs]

ys2 = []
smoothing = 1
for i in range(len(xs) - smoothing * 2):
    ys2.append(sum(ys[i:i+smoothing * 2 + 1]) / (smoothing * 2 + 1))

ys = ys2
xs = xs[smoothing:-smoothing]


fit, cov = curve_fit(curve, xs, ys, p0=[1.8, -0.01], bounds=([0, -1], [2, 0]))

print(fit)
yfits = [curve(x, fit[0], fit[1]) for x in xs]


fig, ax = plt.subplots()

ax.plot(xs, ys, c="red", label="simulated wasserstein")
ax.plot(xs, yfits, "k", label="exponential fit")
ax.legend()

ax.set_yscale("log")
ax.set_ybound(1, 2)
ax.set_xbound(0, 1000)
ax.set_xlabel("Time")
ax.set_ylabel("Simulated Wasserstein Distance")
show()
