"""Validates events working as properly. Tests rate functions and implementations."""

# pylint:disable=missing-function-docstring
from src.tools.models.homogeneous import IndependentEvent, IndependentBirth, IndependentDeath, IndependentSwitch

def test_linear_event_rate():
    event = IndependentEvent(2, lambda x: x['p'])
    parameters = {'q': 1, 'p': 2}
    state = [1, 3, 5]
    assert event.get_max_rate(state, parameters) == 10

def test_birth_implementation():
    birth = IndependentBirth(1, 'b')
    state = [0, 1, 1]
    birth.implement(state)
    assert state == [0, 2, 1]

def test_death_implementation():
    death = IndependentDeath(0, 'd')
    state = [1, 1, 1]
    death.implement(state)
    assert state == [0, 1, 1]

def test_transition_implementation():
    transition = IndependentSwitch(0, 1, lambda x: x['0->1'])
    state = [1, 0]
    parameters = {'0->1': 3}
    assert transition.get_max_rate(state, parameters) == 3
    transition.implement(state)
    assert state == [0, 1]
