"""
Contains functionality to transform simulation results into plots.
"""
import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'


def plot_timepoint_data(result):
    """
    Data is input in the following form:
    {
        "parameters": dict,
        "model": str,
        "data": {
            timepoint : list
        }
        "labels": list
    }


    """

    fig, ax = plt.subplots(figsize = (13, 8))


    # plots graph
    xs = list(result["data"].keys())
    xs.sort()
    site_count = len(result["data"][0]) - 1
    for i in range(site_count + 1):
        ys = []
        for x in xs:
            ys.append(result["data"][x][i])
        if "label" in result:
            label = result["label"][i]
        else:
            label = i
        ax.plot(xs, ys, label=f"{label}")

    # set up x-axis

    min_time = xs[0]
    max_time = xs[-1]

    min_y = 0
    max_y = 1.3 * max(result["data"][0])

    ax.set_xlim(min_time, max_time)
    ax.set_ylim(min_y, max_y)
    ax.set_xticks([min_time, max_time])
    ax.set_xlabel("Time")
    ax.set_ylabel("Population Counts")
    # make parameter table


    info = ""
    table = [["parameter"], ["value"]]
    # table = [["parameter", "value"]]
    for key, val in result["parameters"].items():
        info = f"{info} {key} = {val}\n"
        table[0].append(key)
        table[1].append(str(val))
        # table.append([key, str(val)])
    #table = ax.table(table, loc='top', cellLoc='center')
    #table.auto_set_font_size(False)
    #table.set_fontsize(12)

    plt.text(max_time / 10, 0, f"{info}")
    plt.legend()
    return fig, ax


def plot_extinction_comparison(results_dict):
    """
    Data is input in the following form:
    {
        "results":   
            [
                {
                    "parameters": dict,
                    "model": str,
                    "data": {
                        cell_type : extinction_rate
                    },
                    "label": str,
                }
            ],
        "info": str

    """
    fig, ax = plt.subplots(figsize=(13, 8))

    ax.set_ylim(0, 1)
    ax.set_xlim(0, 1)
    fontsize=18
    ax.set_ylabel("Extinction probability starting with a single cell", fontsize=fontsize)
    ax.set_xlabel("Fraction of sites methylated for that cell", fontsize=fontsize)
    ax.grid(True)
    for result in results_dict["results"]:
        xs = []
        ys = []
        for x, y in result["data"].items():
            xs.append(x)
            ys.append(y)
        ax.plot(xs, ys, label=result["label"], linewidth=2)
    ax.legend()
    plt.text(0.21, 0.015, f'{results_dict["info"]}')
    return fig, ax
        

def plot_dictionary_series(series_data, scale="linear"):
    """
    Data is input as 
    {series_name: {timepoint: float}}
    """
    fig, ax = plt.subplots(figsize=(13, 8))

    for name, series in series_data.items():
        xs = list(series.keys())
        ys = [series[x] for x in xs]
        if scale == "linear":
            ax.plot(xs, ys, label=name, linewidth=2)
        elif scale == "log":
            ax.semilogy(xs, ys, label=name, linewidth=2)
    ax.legend()
    return fig, ax



def show():
    plt.show()
