from src.tools.models.methylation import SingleSite, NoncollaborativeSingleCell
from src.tools.plot import plot_timepoint_data, show
from src.tools.io import save_figure
parameters = {
    "r_uh": 0.1,
    "r_hu": 0.1,
    "r_hm": 1,
    "r_mh": 0.1,
    "b": 0.04,
}



site = SingleSite()

stable_distribution = site.get_stable_state(parameters)
print(stable_distribution)
u = stable_distribution["u"]
h = stable_distribution["h"]
m = stable_distribution["m"]

stable_post_split = {
    "u": u + h/2,
    "h": h/2 + m,
    "m": 0,
}

print(stable_post_split)
cell = NoncollaborativeSingleCell()
timesteps = [i/100 for i in range(30001)]
data = cell.generate_simulation_data(parameters, [100, 0, 0], timesteps)
print(data)

fig, ax = plot_timepoint_data(data)
show()
save_figure(fig, "single_cell")

