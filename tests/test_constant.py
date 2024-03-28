from src.tools.models.event import ConstantEvent, ConstantEventModel
import numpy as np
from scipy.linalg import eig

def test_transition_matrix():
    e1 = ConstantEvent("A", "B", lambda x : x["A->B"])
    e2 = ConstantEvent("B", "A", lambda x : x["B->A"])
    parameters = {
        "A->B": 1,
        "B->A": 1,
    }
    model = ConstantEventModel([e1, e2])
    transition_matrix = model._get_transition_matrix(parameters)
    assert (transition_matrix == np.array([[-1, 1],[1, -1]])).all()


def test_stable_state_1():
    e1 = ConstantEvent("A", "B", lambda x : x["A->B"])
    e2 = ConstantEvent("B", "A", lambda x : x["B->A"])
    parameters = {
        "A->B": 2,
        "B->A": 1,
    }
    model = ConstantEventModel([e1, e2])
    stable_state = model.get_stable_state(parameters)

    assert stable_state == {
        "A": 1/3, 
        "B": 2/3,
    }


def test_stable_state_2():
    e1 = ConstantEvent("A", "B", lambda x : x["A->B"])
    e2 = ConstantEvent("B", "C", lambda x : x["B->C"])
    e3 = ConstantEvent("C", "A", lambda x : x["C->A"])
    parameters = {
        "A->B": 2,
        "B->C": 1,
        "C->A": 1,
    }
    model = ConstantEventModel([e1, e2, e3])
    stable_state = model.get_stable_state(parameters)

    assert stable_state == {
        "A": 1/5, 
        "B": 2/5,
        "C": 2/5,
    }