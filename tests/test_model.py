"""
Tests the model functions to check that they are working properly
"""

# pylint:disable=missing-function-docstring,invalid-name

import numpy as np

from src.tools.model import Birth, Death, LinearModel, Transition


def test_deterministic_run_1():
    e1 = Birth(0, "b")
    e2 = Death(0, "d")
    e3 = Transition(0, 1, "0->1")
    model = LinearModel([e1, e2, e3])
    p = {
        "b": 1,
        "d": 1,
        "0->1": 1,
    }
    initial = [1, 0]
    d = 0
    assert (model.run_deterministic(p, initial, d) == np.array([1, 0])).all()

