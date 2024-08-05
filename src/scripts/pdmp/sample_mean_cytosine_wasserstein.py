from src.tools.models.pdmp import PDMP
from src.tools.plot import show
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from src.constants import BISTABLE_STRONG_PDMP_PARAMS
import json
import os 

output_path = f"output/simulations/{os.path.basename(__file__)[:-3]}.json"


converted_parameters = PDMP.convert_to_cytosine_mean(BISTABLE_STRONG_PDMP_PARAMS)

print(converted_parameters)

times = [t / 1000 for t in range(1000)]
sample_count = 1000
model = PDMP()
if input("resimulate data? y/n: ") == "y":
    result = model.sample_wasserstein(converted_parameters, times, sample_count)
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


xs = sorted(list(result.keys()))
ys = [result[x] for x in xs]

fig, ax = plt.subplots()

ax.plot(xs, ys, c="red", label="simulated wasserstein")
ax.legend()
ax.set_yscale("log")
ax.set_xbound(0, 1)
ax.set_xlabel("Time")
ax.set_ylabel("Wasserstein Distance")
show()
