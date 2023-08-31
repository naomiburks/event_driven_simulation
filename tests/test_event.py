"""Validates events working as properly. Tests rate functions and implementations."""

# pylint:disable=missing-function-docstring
from src.tools import model

def test_linear_event_rate():
    event = model.LinearEvent(2, 'p')
    parameters = {'q': 1, 'p': 2}
    state = [1, 3, 5]
    assert event.get_rate(state, parameters) == 10

def test_birth_implementation():
    birth = model.Birth(1, 'b')
    state = [0, 1, 1]
    birth.implement(state)
    assert state == [0, 2, 1]

def test_death_implementation():
    death = model.Death(0, 'd')
    state = [1, 1, 1]
    death.implement(state)
    assert state == [0, 1, 1]

def test_transition_implementation():
    transition = model.Transition(0, 1, '0->1')
    state = [1, 0]
    parameters = {'0->1': 3}
    assert transition.get_rate(state, parameters) == 3
    transition.implement(state)
    assert state == [0, 1]
