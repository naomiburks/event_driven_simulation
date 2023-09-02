"""
Models are implemented in a parameter-agnostic way. 
Parameters instead are supplied at the time of running.
"""
# pylint:disable=arguments-differ
from copy import deepcopy
from random import random
import numpy as np
from scipy.linalg import expm
from scipy import optimize


class Model:
    """Abstract class. Contains a state space and function to run for a duration."""
    name = "Abstract Model"
    def run(self, parameters, initial_state, duration: float, **kwargs):
        """This is the only important function in a model. """
        # pylint:disable=unused-argument
        return initial_state

    def generate_simulation_data(self, parameters: dict, initial_state, timesteps: list,
                                 **kwargs):
        """
        Useful to run simulations. 
        Data is output in json-style:
        {
            "model": [model name],
            "parameters": [parameter dictionary],
            "data" [timepoint data dictionary],
        }
        """

        simulation_result = {
            "parameters": parameters,
            "model": self.name,
            "data": {0: initial_state},
        }
        last_time = 0
        current_state = initial_state
        for time in timesteps:
            duration = time - last_time
            current_state = self.run(
                parameters, current_state, duration, **kwargs)
            simulation_result["data"][time] = current_state
            last_time = time
        return simulation_result

class Event:
    """
    Abstract class. 

    To instantiate an event, you must override:
     - rate function | (state, parameters) -> nonnegative number
     - implement function | state -> state
    """

    def get_rate(self, state, model_parameters):
        """Outputs a nonnegative number: the frequency the event occurs"""

    def implement(self, state):
        """Mutates the state with the event implemented and returns the state"""
        return state

class EventModel(Model):
    """
    Abstract class to handle event-driven time-independent models.
    Parameters not included in instantiation.

    To understand a model you should be comfortable with some concepts: 
        - A state space. This is the type of data that the model will manipulate.
        - A set of possible events. These events occur stochastically and modify the state.
        - Parameters. A single model may behave in different ways under different parameters.
        This is generally accomplished via having the parameters affect the event rates. 


    In order to instantiate a model, you must specify its events. 
    The state space is implicitly determined. Attempting to run models on states 
    outside its state space will generally lead to errors. 
    The parameters must be provided at runtime.

    In order to run a model, you must specify:
        - parameters
        - initial state
        - duration
    """

    name = "Abstract Model"

    def __init__(self, events: list[Event]):
        """"""
        self.events = events

    def run(self, parameters: dict, initial_state, duration: float, max_num_steps=None):
        """
        Runs the model.

        Does not mutate any of the arguments.  
        """

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
        transition_matrix = self._get_transition_matrix(parameters)
        return initial_state @ expm(duration * transition_matrix)

    def _get_transition_matrix(self, parameters):
        return self.generator_function(parameters)

class LinearEvent(Event):
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
        For example, the per-individual rate could depend on parameters without 
        being one itself.
        """
        return model_parameters[self.rate_parameter_name]


class Birth(LinearEvent):
    """Most commonly birth event in a branching process."""

    def implement(self, state):
        state[self.population_index] = state[self.population_index] + 1
        return state


class Death(LinearEvent):
    """Most commonly death event in a branching process."""

    def implement(self, state):
        state[self.population_index] = state[self.population_index] - 1
        return state


class Transition(LinearEvent):
    """Most commonly transition event in a multitype branching process."""

    def __init__(self, population_index: int, new_population_index: int, rate_parameter_name: str):
        super().__init__(population_index, rate_parameter_name)
        self.new_population_index = new_population_index
        self._necessary_indices.append(self.new_population_index)

    def implement(self, state):
        state[self.population_index] = state[self.population_index] - 1
        state[self.new_population_index] = state[self.new_population_index] + 1
        return state


class LinearModel(EventModel, PopulationModel):
    """
    A linear model is a population model with no interaction between individuals. 
    All event rates depend on a single population and scale linearly with its size.

    PopulationModel is listed after EventModel because super().__init__()
    desired to be EventModel's __init__.
    """

    def __init__(self, events: list[LinearEvent]):
        super().__init__(events)
        population_count = 1 + \
            max([max(event._necessary_indices) for event in events])
        self.population_count = population_count

    def run_deterministic(self, parameters: dict, initial_state, duration: float):
        """
        This calculates the output if the stochastic discrete events were 
        instead differential equations. For linear models, this is equivalent to 
        finding the mean behaviour. 
        """
        generator = self._calculate_generator(parameters)
        end_state = initial_state @ expm(duration * generator)
        return end_state

    def get_deterministic_model(self):
        """Returns the deterministic version of the model"""
        return ExponentialPopulationModel(self.population_count, self._calculate_generator)

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

        for _ in range(10):
            initial_guess = recursive_extinction_function(initial_guess)

        return optimize.fixed_point(recursive_extinction_function, initial_guess)

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
