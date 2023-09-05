from copy import deepcopy
from random import random

import numpy as np

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

class TimeDependentEvent:
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
