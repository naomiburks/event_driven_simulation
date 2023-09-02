"""
This script runs a simulation, saves the results to a file (overwriting the past file),
and plots the simulation.
"""

from src.tools import methylation
from src.tools import plot
from src.tools import io

SIMULATION_PATH = "example"

if __name__ == "__main__":
    site_count = 7
    model = methylation.OneDimensionalNonCollaborativeMethylation(site_count)
    parameters = {
        "r_um": 1,
        "r_mu": 1,
        "b_0": 0,
        "b_M": 0,
        "d_0": 0,
        "d_M": 0,
    }
    timesteps = [i / 1000 for i in range(1001)]
    initial_state = [100] * 8
    #result = model.generate_simulation_data(parameters, initial_state, timesteps)
    #io.write_simulation(result, SIMULATION_PATH)
    result = io.read_simulation(SIMULATION_PATH)
    plot.plot_timepoint_data(result)
    