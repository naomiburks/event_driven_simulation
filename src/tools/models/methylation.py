"""
This module contains tools for modeling methylation. 
In particular, we instantiate:
    - Methylation and Demethylation events, where the rate varies linearly with the population index
    - Interpolated linear events, where we can establish linearly varying birth and death rates
    to simulate a fitness landscape.
"""

import math
import random

from scipy.integrate import odeint
from scipy.optimize import root

from src.tools.models.event import (ConstantEvent, ConstantEventModel, Event,
                                    EventModel)
from src.tools.models.homogeneous import (Birth, Death, HomogeneousEvent,
                                          HomogeneousModel, Switch)
from src.tools.models.model import Model


class OneDimensionalMethylation(Switch):
    """Describes a methylation event"""

    def __init__(self, site_index, site_count, rate_parameter_name):
        super().__init__(site_index, site_index + 1, rate_parameter_name)
        self.site_count = site_count

    def get_rate_per_individual(self, model_parameters):
        """
        Overrides the function in the base class, 
        since methylation rate depends on number of unmethylated sites.
        """
        rate_parameter = model_parameters[self.rate_parameter_name]
        return (self.site_count - self.population_index) * rate_parameter


class OneDimensionalDemethylation(Switch):
    """Describes a demethylation event"""

    def __init__(self, site_index, rate_parameter_name):
        super().__init__(site_index, site_index - 1, rate_parameter_name)

    def get_rate_per_individual(self, model_parameters):
        """
        Overrides the function in the base class, 
        since demethylation rate depends on number of methylated sites.
        """
        rate_parameter = model_parameters[self.rate_parameter_name]
        return self.population_index * rate_parameter


class InterpolatedLinearEvent(HomogeneousEvent):
    """
    The rate of an interpolated linear event is calculated 
    linearly between the rate for pop 0 and pop max.
    """

    def __init__(self, population_index, population_count,
                 rate_min_parameter_name, rate_max_parameter_name):
        super().__init__(population_index, None)
        self.population_count = population_count
        self.rate_min_parameter_name = rate_min_parameter_name
        self.rate_max_parameter_name = rate_max_parameter_name

    def get_rate_per_individual(self, model_parameters):
        rate_max = model_parameters[self.rate_max_parameter_name]
        rate_min = model_parameters[self.rate_min_parameter_name]

        rate = (rate_max * self.population_index +
                rate_min * (self.population_count - self.population_index - 1)) / \
            (self.population_count - 1)
        return rate


class InterpolatedBirth(InterpolatedLinearEvent, Birth):
    """
    Interpolated birth is both a Birth and an InterpolatedLinearEvent
    """


class InterpolatedDeath(InterpolatedLinearEvent, Death):
    """
    Interpolated death is both a Death and an InterpolatedLinearEvent
    """


class OneDimensionalNonCollaborativeMethylation(HomogeneousModel):
    """
    Defines a model with the following properties:
        - M + 1 types: type for each of 0 sites methylated through M sites methylated
        - Methylation (x -> x+1) and Demthyation (x -> x-1) events are noncollaborative. 
            - Rates proportional to the number of unmethylated (or methylated) sites.
        - Birth and death rates follow vary linearly with methylation level. 
    """
    name = "One Dimensional Noncollaborative Methylation with Linear Fitness"

    def __init__(self, M: int):
        events = []
        for i in range(M + 1):
            if i != M:  # add methylations
                events.append(OneDimensionalMethylation(i, M, "r_um"))
            if i != 0:  # add demethylations
                events.append(OneDimensionalDemethylation(i, "r_mu"))
            # add births and deaths
            events.append(InterpolatedBirth(i, M + 1, "b_0", "b_M"))
            events.append(InterpolatedDeath(i, M + 1, "d_0", "d_M"))
        super().__init__(events)
        self.name = f"{self.name} ({self.population_count - 1} sites)"

    def get_limit_model(self):
        def r_um(x, parameters):
            return parameters["r_um"]

        def r_mu(x, parameters):
            return parameters["r_mu"]

        def r_b(x, parameters):
            return x * parameters["b_M"] + (1 - x) * parameters["b_0"]

        def r_d(x, parameters):
            return x * parameters["d_M"] + (1 - x) * parameters["d_0"]
        model = InfiniteSiteOneDimensionalLimit(r_b, r_d, r_um, r_mu)
        return model


class HalfConstantEvent(ConstantEvent):
    def get_max_rate(self, state, model_parameters):
        return super().get_max_rate(state, model_parameters) / 2


class SingleSite(ConstantEventModel):
    def __init__(self):
        events = []
        events.append(ConstantEvent("u", "h", "r_uh"))
        events.append(ConstantEvent("h", "m", "r_hm"))
        events.append(ConstantEvent("h", "u", "r_hu"))
        events.append(ConstantEvent("m", "h", "r_mh"))
        events.append(ConstantEvent("m", "u", "b"))
        events.append(HalfConstantEvent("h", "u", "b"))
        super().__init__(events)


class MethylationBirthShift(Event):
    def __init__(self, rate_parameter_name):
        self.rate_parameter_name = rate_parameter_name

    def get_max_rate(self, state, parameters):
        return parameters[self.rate_parameter_name]

    def implement(self, state):
        """Push the state to the side!"""
        extra_u = 0
        for _ in range(state[1]):
            extra_u += math.floor(2 * random.random())
        state[0] = state[0] + extra_u
        state[1] = state[1] - extra_u + state[2]
        state[2] = 0


class NoncollaborativeSingleCell(EventModel):
    """The state space is (N_u, N_h, N_m) | N_u + N_h + N_m = M"""
    name = "Single Cell Noncollaborative Model"

    def __init__(self):
        events = []
        events.append(Switch(0, 1, "r_uh"))
        events.append(Switch(1, 0, "r_hu"))
        events.append(Switch(1, 2, "r_hm"))
        events.append(Switch(2, 1, "r_mh"))
        events.append(MethylationBirthShift("b"))
        super().__init__(events)


class InfiniteSiteOneDimensionalLimit(Model):
    def __init__(self, r_b, r_d, r_um, r_mu):
        self._r_b = r_b
        self._r_d = r_d
        self._r_um = r_um
        self._r_mu = r_mu

    def calculate_extinction(self, parameters, point_count=1000):
        def extinction_derivative(y, x):
            b = self._r_b(x, parameters)
            d = self._r_d(x, parameters)
            um = self._r_um(x, parameters)
            mu = self._r_mu(x, parameters)
            if methylation_evolution(x) == 0:
                if y == 1:
                    return 0
                return (d - b) / (um + mu)
            return ((b * y - d) * (y - 1)) / (mu * x - um * (1 - x))

        def methylation_evolution(i):
            return self._r_mu(i, parameters) * i - self._r_um(i, parameters) * (1 - i)

        initial_x = root(methylation_evolution, 0.5).x[0]
        initial_y = min(self._r_d(initial_x, parameters) /
                        self._r_b(initial_x, parameters), 1)
        all_times = [i / (point_count - 1) for i in range(point_count)]
        low_times = [time for time in all_times if time <
                     initial_x] + [initial_x]
        high_times = [initial_x] + \
            [time for time in all_times if time >= initial_x]

        # get probabilities below the stable x
        low_res = odeint(extinction_derivative, initial_y,
                         list(reversed(low_times)))

        # get probabilities above the stable x
        high_res = odeint(extinction_derivative, initial_y, high_times)

        # combine probabilities into one list
        res = list(reversed(low_res[:, 0][1:])) + list(high_res[:, 0][1:])

        return res

    