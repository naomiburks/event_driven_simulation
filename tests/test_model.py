"""
Tests the model functions to check that they are working properly
"""

# pylint:disable=missing-function-docstring,invalid-name

import numpy as np

from src.tools.models.homogeneous import Birth, Death, HomogeneousModel, Transition


def test_deterministic_run_1():
    e1 = Birth(0, "b")
    e2 = Death(0, "d")
    e3 = Transition(0, 1, "0->1")
    model = HomogeneousModel([e1, e2, e3]).get_deterministic_model()
    p = {
        "b": 1,
        "d": 1,
        "0->1": 1,
    }
    initial = [1, 0]
    d = 0
    assert (model.run(p, initial, d) == np.array([1, 0])).all()


def test_extinction_1():
    b = Birth(0, "b")
    d = Death(0, "d")
    b2 = Birth(1, "b2")
    d2 = Death(1, "d2")
    t = Transition(0, 1, "0->1")
    events = [b, d, b2, d2, t]
    parameters = {
        "b": 1,
        "d": 0.4,
        "b2": 0.5,
        "d2": 1,
        "0->1": 0.1,
    }

    M = HomogeneousModel(events)

    probabilities = M.calculate_extinction(parameters)

    # the correct answer is [0.5, 1] up to floating point errors

    assert probabilities[0] > 0.499999999999
    assert probabilities[0] < 0.500000000001
    assert probabilities[1] > 0.999999999999
    assert probabilities[1] < 1.000000000001


def test_extinction_2():
    b = Birth(0, "b")
    d = Death(0, "d")
    events = [b, d]
    parameters = {
        "b": 1,
        "d": 0.4,
    }

    M = HomogeneousModel(events)

    probability = M.calculate_extinction(parameters)[0]

    assert probability < 0.40000000001
    assert probability > 0.39999999999


def test_extinction_3():
    b1 = Birth(0, "b1")
    b2 = Birth(1, "b2")
    d1 = Death(0, "d1")
    d2 = Death(1, "d2")
    t1 = Transition(0, 1, "0->1")
    t2 = Transition(1, 0, "1->0")
    parameters = {
        "b1": 2,
        "b2": 1,
        "d1": 1,
        "d2": 2,
        "0->1": 1,
        "1->0": 1,
    }
    events = [b1, b2, d1, d2, t1, t2]
    M = HomogeneousModel(events)
    probabilities = M.calculate_extinction(parameters)
    true_probabilities = [0.7632670961354369005, 0.8879150644557030511]
    # calculated as solution to system of equations via wolfram alpha

    for p, true_p in zip(probabilities, true_probabilities):
        assert p < true_p + 0.0000001
        assert p > true_p - 0.0000001


def test_extinction_4():
    b1 = Birth(0, "b0")
    b2 = Birth(1, "b1")
    b3 = Birth(2, "b2")
    d1 = Death(0, "d0")
    d2 = Death(1, "d1")
    d3 = Death(2, "d2")
    t1 = Transition(2, 1, "2->1")
    t2 = Transition(1, 0, "1->0")
    parameters = {
        "b0": 3,
        "d0": 1,
        "b1": 2,
        "d1": 2,
        "1->0": 3,
        "b2": 9,
        "d2": 1,
        "2->1": 8,
    }
    events = [b1, b2, b3, d1, d2, d3, t1, t2]
    M = HomogeneousModel(events)
    probabilities = M.calculate_extinction(parameters)
    true_probabilities = [1/3, 1/2, 1/3]  # calculated by hand :(
    for p, true_p in zip(probabilities, true_probabilities):
        assert p < true_p + 0.0000001
        assert p > true_p - 0.0000001
