"""
Models are implemented in a parameter-agnostic way. 
Parameters instead are supplied at the time of running.
"""
from copy import deepcopy
from random import random
import numpy as np
from scipy.linalg import expm


class Event:
    """
    Abstract class. 
    Event consists of:
     - rate function | (state, parameters) -> nonnegative
     - implementation function | state -> state
    """

    def get_rate(self, state, model_parameters):
        """Outputs a nonnegative number: the frequency the event occurs"""

    def implement(self, state):
        """Mutates the state with the event implemented"""


class Model:
    """
    Abstract class to hand event-driven time-independent models.
    Parameters not included in instantiation.
    """
    model_name = "Abstract Model"
    def __init__(self, events: list[Event]):
        """"""
        self.events = events

    def run(self, parameters: dict, initial_state, duration: float, max_num_steps=None):
        current_state = deepcopy(initial_state)
        current_time = 0
        num_steps = 0
        while True:
            num_steps += 1
            if max_num_steps is not None and num_steps > max_num_steps:
                raise RuntimeError(
                    "Maximum number of steps for single simulation exceeded")
            rates = []
            for event in self.events:
                rates.append(event.get_rate(current_state, parameters))
            total_rate = sum(rates)
            if total_rate == 0:
                break
            else:
                waiting_time = - np.log(random()) / total_rate
            current_time += waiting_time
            if current_time > duration:
                break
            event_index = random() * total_rate
            found_event = None
            for event, rate in zip(self.events, rates):
                if rate >= event_index:
                    found_event = event
                    break
                event_index -= rate
            if found_event is None:
                raise RuntimeError("Event was not able to be found!")
            event.implement(current_state)
        return current_state

    def generate_simulation_data(self, parameters: dict, initial_state, timesteps: list,
                                 max_num_steps=None):
        """
        Useful to run simulations. 
        Data is output all-in-one: the type of model, its parameters, and the timepoint data.
        You do NOT get the model's namespace variables!
        """

        simulation_result = {
            "parameters": parameters,
            "model": self.model_name,
            "data": {0: initial_state},
        }
        last_time = 0
        current_state = initial_state
        for time in timesteps:
            duration = time - last_time
            current_state = self.run(
                parameters, current_state, duration, max_num_steps=max_num_steps)
            simulation_result["data"][time] = current_state
            last_time = time
        return simulation_result


class PopulationModel(Model):
    """In a population model, the state space is a list of population counts."""

    def __init__(self, population_count: int, events: list[Event]):
        super().__init__(events)
        self.population_count = population_count

    def sample_extinction(self, parameters, duration, num_attempts_per_pop):
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


class LinearEvent(Event):
    """
    These events have rates varying linearly with population size.
    This class is abstract: an instantiation must say what implementing the event
    does to the state.
    """
    def __init__(self, population_index: int, rate_parameter_name: str):
        super().__init__()
        self.population_index = population_index
        self.rate_parameter_name = rate_parameter_name
        self._necessary_indices = [population_index]

    def get_rate(self, state, model_parameters):
        """
        This calculates the frequency the event occurs given a state and parameters
        """
        rate_per_individual = self.get_rate_per_individual(model_parameters)
        return rate_per_individual * state[self.population_index]

    def get_rate_per_individual(self, model_parameters):
        """
        By default, the per-individual rate is given as a parameter. 
        Override this function for other ways to calculate rate.
        """
        return model_parameters[self.rate_parameter_name]


class Birth(LinearEvent):
    def implement(self, state):
        state[self.population_index] = state[self.population_index] + 1


class Death(LinearEvent):
    def implement(self, state):
        state[self.population_index] = state[self.population_index] - 1


class Transition(LinearEvent):
    def __init__(self, population_index: int, new_population_index: int, rate_parameter_name: str):
        super().__init__(population_index, rate_parameter_name)
        self.new_population_index = new_population_index
        self._necessary_indices.append(self.new_population_index)

    def implement(self, state):
        state[self.population_index] = state[self.population_index] - 1
        state[self.new_population_index] = state[self.new_population_index] + 1


class LinearModel(PopulationModel):
    """
    A linear model is a population model with no collaboration between individuals. 
    This means that all events should be linear. 
    """

    def __init__(self, events: list[LinearEvent]):
        population_count = 1 + \
            max([max(event._necessary_indices) for event in events])
        super().__init__(population_count, events)

    def run_deterministic(self, parameters: dict, initial_state, duration: float):
        generator = self._calculate_generator(parameters)
        end_state = initial_state @ expm(duration * generator)
        return end_state

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
