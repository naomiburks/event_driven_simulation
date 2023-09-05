from src.tools.models import methylation

def test_interpolated_linear_event_1():
    event = methylation.InterpolatedLinearEvent(1, 3, "r_min", "r_max")
    parameters = {
        "r_min": 0,
        "r_max": 1,
    }
    state = [1, 1, 1]
    rate = event.get_rate(state, parameters)
    assert rate == 0.5

def test_interpolated_linear_event_2():
    event = methylation.InterpolatedLinearEvent(2, 7, "r_min", "r_max")
    parameters = {
        "r_min": 0,
        "r_max": 1,
    }
    state = [1, 1, 1, 0, 0, 0, 0]
    rate = event.get_rate(state, parameters)
    assert rate == 1/3

def test_interpolated_birth_1():
    event = methylation.InterpolatedBirth(2, 7, "r_min", "r_max")
    parameters = {
        "r_min": 0,
        "r_max": 1,
    }
    state = [1, 1, 1, 0, 0, 0, 0]
    rate = event.get_rate(state, parameters)
    assert rate == 1/3

def test_interpolated_birth_2():
    event = methylation.InterpolatedBirth(2, 7, "r_min", "r_max")
    state = [1, 1, 1, 0, 0, 0, 0]
    event.implement(state)
    assert state == [1, 1, 2, 0, 0, 0, 0]

def test_interpolated_death_1():
    event = methylation.InterpolatedDeath(2, 7, "r_min", "r_max")
    parameters = {
        "r_min": 0,
        "r_max": 1,
    }
    state = [1, 1, 1, 0, 0, 0, 0]
    rate = event.get_rate(state, parameters)
    assert rate == 1/3

def test_interpolated_death_2():
    event = methylation.InterpolatedDeath(2, 7, "r_min", "r_max")
    state = [1, 1, 1, 0, 0, 0, 0]
    event.implement(state)
    assert state == [1, 1, 0, 0, 0, 0, 0]
