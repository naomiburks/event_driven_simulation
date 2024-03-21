from src.tools.models.event import TimeIndependentEvent, Event, EventModel
from src.tools.models.population import PopulationModel, ExponentialPopulationModel
import numpy as np
from scipy import optimize
from src.constants import CONVERGENCE_TOLERANCE

class IndependentEvent(TimeIndependentEvent):
    """
    These events have rates varying linearly with population size as is common in
    (multitype) branching processes.
    This class is abstract: an instantiation must say what implementing the event
    does to the state.
    """

    def __init__(self, population_index: int, rate_parameter_name: str):
        super().__init__()
        self.population_index = population_index
        self.rate_parameter_name = rate_parameter_name
        self._necessary_indices = [population_index]

    def get_max_rate(self, state, model_parameters):
        """
        This calculates the frequency the event occurs given a state and parameters
        """
        rate_per_individual = self.get_rate_per_individual(model_parameters)
        return rate_per_individual * state[self.population_index]

    def get_rate_per_individual(self, model_parameters):
        """
        By default, the per-individual rate is given as a parameter. 
        Override this function for other ways to calculate rate. 
        For example, the per-individual rate could depend on parameters without 
        being one itself.
        """
        return model_parameters[self.rate_parameter_name]


class Birth(IndependentEvent):
    """Most commonly birth event in a branching process."""

    def implement(self, state):
        state[self.population_index] = state[self.population_index] + 1
        return state


class Death(IndependentEvent):
    """Most commonly death event in a branching process."""

    def implement(self, state):
        state[self.population_index] = state[self.population_index] - 1
        return state


class Switch(IndependentEvent):
    """Most commonly transition event in a multitype branching process."""

    def __init__(self, population_index: int, new_population_index: int, rate_parameter_name: str):
        super().__init__(population_index, rate_parameter_name)
        self.new_population_index = new_population_index
        self._necessary_indices.append(self.new_population_index)

    def implement(self, state):
        state[self.population_index] = state[self.population_index] - 1
        state[self.new_population_index] = state[self.new_population_index] + 1
        return state


class IndependentModel(EventModel, PopulationModel):
    """
    A homogeneous model is a population model with no interaction between individuals. 
    All event rates depend on a single population and scale linearly with its size.

    PopulationModel is listed after EventModel because super().__init__()
    desired to be EventModel's __init__.
    """
    name = "Linear Model"

    def __init__(self, events: list[IndependentEvent]):
        super().__init__(events)
        population_count = 1 + \
            max([max(event._necessary_indices) for event in events])
        self.population_count = population_count

    def get_deterministic_model(self) -> ExponentialPopulationModel:
        """Returns the deterministic version of the model"""
        model = ExponentialPopulationModel(
            self.population_count, self._calculate_generator)
        model.name = f"{self.name} (deterministic)"
        return model

    def calculate_extinction(self, parameters: dict):
        """
        Calculates the extinction probabilities by solving a 
        first-step conditioning self-similarity equation. 
        Solves the equation via a fixed point scipy solver after initial burn-in self-composition.

        To optimize runtime, we start by caching the following relevant information for each event: 
            - the type of individual that induces the event
            - probability of occurance from the individual
            - result of the event on the individual

        The recursive extinction function maps guess -> first-step conditioning guess. 
        The true extinction probability is a fixed point of this function. This fixed point is the
        guaranteed pointwise limit of iterated self-composition on any nontrivial starting guess.
        """

        total_rates = [0] * self.population_count
        for event in self.events:
            rate = event.get_rate_per_individual(parameters)
            population = event.population_index
            total_rates[population] = total_rates[population] + rate
        probabilities = []
        impacts = []
        populations = []
        for event in self.events:
            population = event.population_index
            impact = event.implement(self._standard_basis_vector(
                population, self.population_count))
            impacts.append(impact)
            probabilities.append(event.get_rate_per_individual(
                parameters) / total_rates[population])
            populations.append(population)

        def recursive_extinction_function(initial_guess):
            # record contribution of each event by rates
            new_guess = [0] * self.population_count
            for probability, impact, population in zip(probabilities, impacts, populations):
                contribution = probability
                for population_count, population_guess in zip(impact, initial_guess):
                    contribution *= population_guess ** population_count
                new_guess[population] = new_guess[population] + contribution
            return np.array(new_guess)

        initial_guess = np.array([0.5] * self.population_count)

        
        return optimize.fixed_point(recursive_extinction_function, initial_guess, 
                                    maxiter=50000, method="iteration", xtol = CONVERGENCE_TOLERANCE)

    def _calculate_generator(self, parameters: dict):
        generator = np.zeros(
            shape=(self.population_count, self.population_count))
        for event in self.events:
            rate = event.get_rate_per_individual(parameters)
            population_index = event.population_index
            effect_vector = self._standard_basis_vector(
                population_index, self.population_count)
            event.implement(effect_vector)
            generator[population_index,
                      population_index] = generator[population_index, population_index] - rate
            for i, entry in enumerate(effect_vector):
                generator[population_index,
                          i] = generator[population_index, i] + rate * entry
        return generator
