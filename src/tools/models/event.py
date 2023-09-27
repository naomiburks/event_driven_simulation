from copy import deepcopy
from random import random

import numpy as np
from scipy.linalg import solve
from src.tools.models.model import Model


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


class ConstantEvent(Event):
    def __init__(self, initial_state, end_state, rate_parameter_name):
        self.initial_state = initial_state
        self.end_state = end_state
        self.rate_parameter_name = rate_parameter_name

    def get_rate(self, state, model_parameters):
        if state == self.initial_state:
            return model_parameters[self.rate_parameter_name]
        return 0

    def implement(self, state):
        return self.end_state


class TimeDependentEvent(Event):
    """
    Time-dependent events have rates that depend on time rather than on the state. 
    It is more complicated to sample from them.
    """


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

    name = "Event-Driven Model"

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


class ConstantEventModel(EventModel):
    def __init__(self, events: list[ConstantEvent]):
        self.state_space = []
        for event in events:
            for state in [event.initial_state, event.end_state]:
                if state not in self.state_space:
                    self.state_space.append(state)

        super().__init__(events)

    def get_stable_state(self, parameters):
        """
        Finds the left eigenvector of the transition matrix whose eigenvalue is 0 with appropriate norm.
        Accomplishes this by solving the equation Ax=b, where b = [1, 0, .., 0] and A is the matrix
        that takes transpose(transition matrix) and replaces the first row with [1, 0, ..., 0].

        Will return a Singular Matrix error if there the space of these eigenvectors is 2+ dimensional,
        which corresponds to having multiple stable states (only possible for reducible markov process). 
        """

        matrix_to_solve = self._get_transition_matrix(parameters).T
        row = np.zeros(len(self.state_space))
        row[0] = 1
        matrix_to_solve[0] = row
        unscaled_stable_vector = solve(matrix_to_solve, row)
        scale = sum(list(unscaled_stable_vector))
        stable_vector = list(unscaled_stable_vector / scale)

        stable_probabilities = {}
        for state, probability in zip(self.state_space, stable_vector):
            stable_probabilities[state] = probability
        return stable_probabilities

    def _get_transition_matrix(self, parameters):
        state_count = len(self.state_space)
        transition_matrix = np.zeros([state_count, state_count])
        for event in self.events:
            initial_index = self._get_list_index(
                self.state_space, event.initial_state)
            end_index = self._get_list_index(self.state_space, event.end_state)
            rate = event.get_rate(event.initial_state, parameters)
            transition_matrix[initial_index,
                              end_index] = transition_matrix[initial_index, end_index] + rate
            transition_matrix[initial_index,
                              initial_index] = transition_matrix[initial_index, initial_index] - rate
        return transition_matrix

    @staticmethod
    def _get_list_index(items, item):
        for i, entry in enumerate(items):
            if entry == item:
                return i
