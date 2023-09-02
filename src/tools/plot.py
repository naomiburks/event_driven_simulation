"""
Contains functionality to transform simulation results into plots.
"""


import matplotlib.pyplot as plt

def plot_timepoint_data(result):
    """
    Data is input in the following form:
    {
        "parameters": dict,
        "model": str,
        "data": {
            timepoint : list
        }
    }
    """

    fig, ax = plt.subplots()
    xs = list(result["data"].keys())
    xs.sort()
    site_count = len(result["data"][0]) - 1
    for i in range(site_count + 1):
        ys = []
        for x in xs: 
            ys.append(result["data"][x][i])
        ax.plot(xs, ys, label=f"{i}")
    plt.legend()
    plt.show()


