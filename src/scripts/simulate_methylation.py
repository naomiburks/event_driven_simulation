from src.tools import methylation
import matplotlib.pyplot as plt


if __name__ == "__main__":
    site_count = 7
    model = methylation.OneDimensionalMethylation(site_count)
    parameters = {
        "r_um": 1,
        "r_mu": 1,
        "b_0": 0,
        "b_M": 0,
        "d_0": 0,
        "d_M": 0,
    }
    timesteps = [i / 100 for i in range(1001)]
    initial_state = [100] * 8
    result = model.generate_simulation_data(parameters, initial_state, timesteps)
    

    fig, ax = plt.subplots()
    xs =  list(result["data"].keys())
    xs.sort()
    for i in range(site_count + 1):
        ys = []
        for x in xs: 
            ys.append(result["data"][x][i])
        ax.plot(xs, ys, label=f"{i}")
    plt.legend()
    plt.show()
    print(result)