import random
from math import log
from scipy.optimize import minimize_scalar


def generate_exponential_waiting_time(rate):
    return -log(random.random()) / rate

def generate_poisson(rate, duration):
    t = 0
    times = []
    while True:
        t += generate_exponential_waiting_time(rate)
        if t > duration:
            return times
        times.append(t)

def generate_poisson_nonhom(rate_func, duration, max_rate = None):
    times = []
    if max_rate is None:
        min_result = minimize_scalar(lambda x : - rate_func(x), duration / 2, bounds = (0, duration))
        max_rate = rate_func(min_result.x)
    potential_times = generate_poisson(max_rate, duration)
    for time in potential_times:
        true_rate = rate_func(time)
        if random.random() < true_rate / max_rate:
            times.append(time)
    return times
