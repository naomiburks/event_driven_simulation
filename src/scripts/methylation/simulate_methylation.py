"""
This script runs a simulation, saves the results to a file (overwriting the past file),
and plots the simulation.
"""

from src.tools.models import methylation
from src.tools import plot

SIMULATION_PATH = "example"

if __name__ == "__main__":
    SITE_COUNT = 7
    model = methylation.OneDimensionalNonCollaborativeMethylation(SITE_COUNT)
    det_model = model.get_deterministic_model()
    parameters = {
        "r_um": 1,
        "r_mu": 1,
        "b_0": 0.2,
        "b_M": 0.3,
        "d_0": 0.3,
        "d_M": 0.2,
    }
    timesteps = [i / 1000 for i in range(1001)]
    initial_state = [100] * 8
    result = model.generate_simulation_data(
        parameters, initial_state, timesteps)
    det_result = det_model.generate_simulation_data(
        parameters, initial_state, timesteps)
    #fig = plot.plot_timepoint_data(result)
    #plot.show()
    fig = plot.plot_timepoint_data(det_result)
    plot.show()