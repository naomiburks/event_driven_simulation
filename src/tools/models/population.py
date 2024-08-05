from scipy.linalg import expm, eig

from src.tools.models.model import Model
import matplotlib.pyplot as plt
import numpy as np


class PopulationModel(Model):
    """In a population model, the state is a list of population counts or fractions
    for populations of various types."""

    def __init__(self, num_of_populations: int):
        super().__init__()
        self.num_of_populations = num_of_populations

    def sample_extinction(self, parameters, duration, num_attempts_per_pop):
        """Uses Monte Carlo to estimate extinction probabilities."""
        extinction_rates = []
        for population_index in range(self.num_of_populations):
            num_extinctions = 0
            initial_state = self._standard_basis_vector(
                population_index, self.num_of_populations)
            for _ in range(num_attempts_per_pop):
                final_state = self.run(parameters, initial_state, duration)
                if sum(final_state) == 0:
                    num_extinctions += 1
            extinction_rates.append(num_extinctions / num_attempts_per_pop)
        return extinction_rates

    @staticmethod
    def _standard_basis_vector(index, length):
        vector = [0] * length
        vector[index] = 1
        return vector

    def plot_sample(self, result, labels):
        fig, ax = plt.subplots()
        timepoint_series = result["data"][0]
        xs = list(timepoint_series.keys())
        population_count = len(timepoint_series[xs[0]])
        ys_list = [[timepoint_series[x][pop] for x in xs] for pop in range(population_count)]
        for ys in ys_list:
            ax.plot(xs, ys, linewidth = 1, alpha = 1)
        leg = ax.legend(labels)
        for line in leg.get_lines():
            line.set_linewidth(4)
        return fig, ax

    def plot_stacked(self, result, labels):
        fig, ax = plt.subplots()
        timepoint_series = result["data"][0]
        xs = list(timepoint_series.keys())
        population_count = len(timepoint_series[xs[0]])
        ys_list = [[sum([timepoint_series[x][smaller_pop] for smaller_pop in range(pop + 1)]) for x in xs] for pop in range(population_count)]
        for i, ys in enumerate(ys_list):
            if i == 0:
                ax.fill_between(xs, ys, facecolor="gyr"[i])
            else:
                ax.fill_between(xs, ys, ys_list[i-1], facecolor="gyr"[i])
        leg = ax.legend(labels)
        for line in leg.get_lines():
            line.set_linewidth(4)
        return fig, ax

    def plot_triangle(self, result, labels):
        fig, ax = plt.subplots()
        timepoint_series = result["data"][0]
        times = list(timepoint_series.keys())
        xs = [timepoint_series[t][0] for t in times]
        ys = [timepoint_series[t][2] for t in times]
        print(labels)
        input()
        print(xs)
        input()
        print(ys)
        ax.plot(xs, ys)
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[2])
        return fig, ax
    
class ExponentialPopulationModel(PopulationModel):
    """
    This model is used to simulate populations growing in accordance
    with the exponential of an instananeous generator matrix M.

    The the entry M[i][j] is the rate at which a type-i individual generates type-j individuals.
    """

    def __init__(self, population_count, get_generator_from_parameters):
        super().__init__(population_count)
        self.get_generator_from_parameters = get_generator_from_parameters

    def run(self, parameters, initial_state, duration):
        transition_matrix = self._get_instantaneous_transition_matrix(parameters)
        return list(initial_state @ expm(duration * transition_matrix))

    def get_long_term_behavior(self, parameters):
        """Returns a double (growth rate: float, stable population fractions: list)"""
        transition_matrix = self._get_instantaneous_transition_matrix(parameters)
        vals, vecs = eig(transition_matrix, left=True, right=False)
        max_val = None
        for i, val in enumerate(vals):
            test_val = np.real(val)
            if max_val is None or test_val > max_val:
                max_val = test_val
                max_index = i
       
        vec = vecs[:, max_index]
        return max_val, list(vec / sum(vec))

    def _get_instantaneous_transition_matrix(self, parameters):
        return self.get_generator_from_parameters(parameters)
