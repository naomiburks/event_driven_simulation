from src.tools.plot import plot_timepoint_data, show
from src.tools.io import save_figure
from src.tools.models.methylation import OneDimensionalNonCollaborative



parameters = {
    "r_um": 0.2,
    "r_mu": 0.1,
    "b_0": 0.18,
    "b_M": 0.18,
    "d_0": 0.3,
    "d_M": 0.1,
}
SITE_COUNT = 5
model = OneDimensionalNonCollaborative(SITE_COUNT)
print(f"extinction: {model.calculate_extinction(parameters)}")
model2 = model.get_deterministic_model()
initial_state = [0] * (SITE_COUNT + 1)
initial_state[0] = 100
timesteps = [i/10 for i in range(1001)]
data = model.generate_simulation_data(parameters, initial_state, timesteps)


fig, ax = plot_timepoint_data(data)
show()

filename = input("file name: ")
save_figure(fig, filename)

