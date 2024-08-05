from src.tools.models.pdmp import PDMP
from src.tools.plot import plot_dictionary_series
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import logging


logging.basicConfig(filename="output/scripts_log.txt", level=logging.INFO)
logger = logging.getLogger(__name__)


x_max = 1000
timepoints = 1000
model = PDMP()
initial_state = [0, 0, 1]
timepoints = [i * x_max /(timepoints - 1)  for i in range(timepoints)]
replicates = 10

def analyze_bistability(data, threshold, tolerance):
    regions = []
    in_region = False
    for key, val in data.items():
        m = val[0]
        if not in_region and m > threshold:
            start_of_region = key
            in_region = True
            end_of_region = key + tolerance
        if in_region and m > threshold:
            end_of_region = key + tolerance
        if in_region and key > end_of_region:
            if end_of_region - start_of_region > 6 * tolerance:
                regions.append([start_of_region, end_of_region - tolerance])
            in_region = False
    if in_region and x_max - start_of_region > 5 * tolerance:
        regions.append([start_of_region, x_max])
    return regions



for i in range(replicates):
    print(f"starting replicate {i}")
    params = {
        'b': 1,
        'r_uh': random.random(),
        'r_hm': random.random(), 
        'r_mh': random.random(),
        'r_hu': random.random(), 
        'r_uh_h': random.random() * 10,
        'r_uh_m': random.random() * 10,
        'r_hm_h': random.random() * 10,
        'r_hm_m': random.random() * 10,
        'r_mh_h': random.random() * 10,
        'r_mh_u': random.random() * 10,
        'r_hu_h': random.random() * 10,
        'r_hu_u': random.random() * 10,
    }
    logger.info(f"set params: {params}") 
    result = model.generate_simulation_data(params, initial_state, timepoints)
    logger.info(f"result: {result}")
    stability = analyze_bistability(result["data"][0], 0.6, 2)
    logger.info(f"stability: {stability}")



#for bounds in stability:
#    rect = patches.Rectangle((bounds[0], 0), bounds[1] - bounds[0], 1, facecolor="k", alpha=0.2, zorder=2)
#    ax.add_patch(rect)


#switching_regions = len(stability) * 2
#if stability[0][0] == 0:
#    switching_regions -= 1
#if stability[-1][1] < x_max:
#    switching_regions += 1

#print(sum(bounds[1] - bounds[0] for bounds in stability) / x_max) 
#print(switching_regions)

#a = sum(bounds[1] - bounds[0] for bounds in stability) / x_max
#b = switching_regions - 1

#r_mu = a * b / x_max
#r_um = (1 - a) * b / x_max

#print(f"estimated transition rates: {r_mu}, {r_um}")
