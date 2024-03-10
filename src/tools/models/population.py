from scipy.linalg import expm, eig

from src.tools.models.model import Model
import numpy as np

class PopulationModel(Model):
    """In a population model, the state space is a list of population counts
    for populations of various types."""

    def __init__(self, population_count: int):
        super().__init__()
        self.population_count = population_count

    def sample_extinction(self, parameters, duration, num_attempts_per_pop):
        """Uses Monte Carlo to estimate extinction probabilities."""
        extinction_rates = []
        for population_index in range(self.population_count):
            num_extinctions = 0
            initial_state = self._standard_basis_vector(
                population_index, self.population_count)
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


class ExponentialPopulationModel(PopulationModel):
    """
    This model is used to run deterministic model where populations grow in accordance
    with the exponential of a generator matrix.
    """

    def __init__(self, population_count, generator_function):
        super().__init__(population_count)
        self.generator_function = generator_function

    def run(self, parameters, initial_state, duration):
        transition_matrix = self._get_instantaneous_transition_matrix(parameters)
        return initial_state @ expm(duration * transition_matrix)

    def get_long_term_behavior(self, parameters):
        """Returns (growth rate, stable state)"""
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
        return self.generator_function(parameters)


class ContinuumPopulationModel(Model):
    """The state space of a continuum population model is a list of
    individuals, where each individual has a type between 0 and 1.

    Individuals have a deterministic drift, as well as stochastic birth and death rates. 
    All are paremeter-dependent functions (type, params) -> num.
    """

    def __init__(self, drift, birthrate, deathrate):
        super().__init__()
        self.drift = drift
        self.birthrate = birthrate
        self.deathrate = deathrate

