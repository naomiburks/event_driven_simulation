from src.tools.models import methylation

def test_1d_birth_1():
    event = methylation.OneDimensionalBirth(2, 6)
    parameters = {
        "b_0": 0,
        "b_M": 1,
    }
    state = [1, 1, 1, 0, 0, 0, 0]
    rate = event.get_max_rate(state, parameters)
    assert rate == 1/3

def test_1d_birth_2():
    event = methylation.OneDimensionalBirth(2, 6)
    state = [1, 1, 1, 0, 0, 0, 0]
    event.implement(state)
    assert state == [1, 1, 2, 0, 0, 0, 0]

def test_1d_death_1():
    event = methylation.OneDimensionalDeath(2, 6)
    parameters = {
        "d_0": 0,
        "d_M": 1,
    }
    state = [1, 1, 1, 0, 0, 0, 0]
    rate = event.get_max_rate(state, parameters)
    assert rate == 1/3

def test_1d_death_2():
    event = methylation.OneDimensionalDeath(2, 6)
    state = [1, 1, 1, 0, 0, 0, 0]
    event.implement(state)
    assert state == [1, 1, 0, 0, 0, 0, 0]

def test_noncoll_methylation_1():
    event = methylation.OneDimensionalNonCollaborativeMethylation(2, 4)
    state = [1, 0, 2, 0, 0]
    parameters = {"r_um": 4}
    rate = event.get_max_rate(state, parameters)
    assert rate == 16

def test_noncoll_methylation_2():
    event = methylation.OneDimensionalNonCollaborativeMethylation(2, 4)
    state = [1, 0, 2, 0, 0]
    event.implement(state)
    assert state == [1, 0, 1, 1, 0]


def test_noncoll_demethylation_1():
    event = methylation.OneDimensionalNonCollaborativeDemethylation(4, 4)
    state = [1, 0, 2, 0, 1]
    parameters = {"r_mu": 2}
    rate = event.get_max_rate(state, parameters)
    assert rate == 8

def test_noncoll_demethylation_2():
    event = methylation.OneDimensionalNonCollaborativeDemethylation(2, 4)
    state = [1, 0, 2, 0, 0]
    event.implement(state)
    assert state == [1, 1, 1, 0, 0]