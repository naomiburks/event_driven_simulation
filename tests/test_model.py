"""
Tests the model functions to check that they are working properly
"""

# pylint:disable=missing-function-docstring,invalid-name

import numpy as np

from src.tools.models.homogeneous import IndependentBirth, IndependentDeath, IndependentModel, IndependentSwitch

from src.constants import CONVERGENCE_TOLERANCE

def test_deterministic_run_1():
    e1 = IndependentBirth(0, lambda x: x["b"])
    e2 = IndependentDeath(0, lambda x: x["d"])
    e3 = IndependentSwitch(0, 1, lambda x: x["0->1"])
    model = IndependentModel([e1, e2, e3]).get_deterministic_model()
    p = {
        "b": 1,
        "d": 1,
        "0->1": 1,
    }
    initial = [1, 0]
    d = 0
    assert (model.run(p, initial, d) == np.array([1, 0])).all()


def test_extinction_1():
    b = IndependentBirth(0, lambda x: x["b"])
    d = IndependentDeath(0, lambda x: x["d"])
    b2 = IndependentBirth(1, lambda x: x["b2"])
    d2 = IndependentDeath(1, lambda x: x["d2"])
    t = IndependentSwitch(0, 1, lambda x: x["0->1"])
    events = [b, d, b2, d2, t]
    parameters = {
        "b": 1,
        "d": 0.4,
        "b2": 0.5,
        "d2": 1,
        "0->1": 0.1,
    }

    M = IndependentModel(events)

    probabilities = M.calculate_extinction(parameters)

    # the correct answer is [0.5, 1]
    assert probabilities[0] > 0.5 - CONVERGENCE_TOLERANCE
    assert probabilities[0] < 0.5 + CONVERGENCE_TOLERANCE
    assert probabilities[1] > 1 - CONVERGENCE_TOLERANCE
    assert probabilities[1] < 1 + CONVERGENCE_TOLERANCE


def test_extinction_2():
    b = IndependentBirth(0, lambda x: x["b"])
    d = IndependentDeath(0, lambda x: x["d"])
    events = [b, d]
    parameters = {
        "b": 1,
        "d": 0.4,
    }

    M = IndependentModel(events)

    probability = M.calculate_extinction(parameters)[0]

    assert probability < 0.4 + CONVERGENCE_TOLERANCE
    assert probability > 0.4 - CONVERGENCE_TOLERANCE


def test_extinction_3():
    b1 = IndependentBirth(0, lambda x: x["b1"])
    b2 = IndependentBirth(1, lambda x: x["b2"])
    d1 = IndependentDeath(0, lambda x: x["d1"])
    d2 = IndependentDeath(1, lambda x: x["d2"])
    t1 = IndependentSwitch(0, 1, lambda x: x["0->1"])
    t2 = IndependentSwitch(1, 0, lambda x: x["1->0"])
    parameters = {
        "b1": 2,
        "b2": 1,
        "d1": 1,
        "d2": 2,
        "0->1": 1,
        "1->0": 1,
    }
    events = [b1, b2, d1, d2, t1, t2]
    M = IndependentModel(events)
    probabilities = M.calculate_extinction(parameters)
    true_probabilities = [0.7632670961354369005, 0.8879150644557030511]
    # calculated as solution to system of equations via wolfram alpha

    for p, true_p in zip(probabilities, true_probabilities):
        assert p < true_p + 0.00001
        assert p > true_p - 0.00001


def test_extinction_4():
    b1 = IndependentBirth(0, lambda x: x["b0"])
    b2 = IndependentBirth(1, lambda x: x["b1"])
    b3 = IndependentBirth(2, lambda x: x["b2"])
    d1 = IndependentDeath(0, lambda x: x["d0"])
    d2 = IndependentDeath(1, lambda x: x["d1"])
    d3 = IndependentDeath(2, lambda x: x["d2"])
    t1 = IndependentSwitch(2, 1, lambda x: x["2->1"])
    t2 = IndependentSwitch(1, 0, lambda x: x["1->0"])
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
    M = IndependentModel(events)
    probabilities = M.calculate_extinction(parameters)
    true_probabilities = [1/3, 1/2, 1/3]  # calculated by hand :(
    for p, true_p in zip(probabilities, true_probabilities):
        assert p < true_p + CONVERGENCE_TOLERANCE
        assert p > true_p - CONVERGENCE_TOLERANCE
