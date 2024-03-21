"""
This module contains tools for modeling methylation. 
In particular, we instantiate:
    - Methylation and Demethylation events, where the rate varies linearly with the population index
    - Interpolated linear events, where we can establish linearly varying birth and death rates
    to simulate a fitness landscape.
"""

import math
import random
import copy 

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import root, minimize

from src.tools.models.event import (ConstantEvent, ConstantEventModel, Event,
                                    EventModel)
from src.tools.models.homogeneous import (Birth, Death, IndependentEvent,
                                          IndependentModel, Switch)
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


class InterpolatedLinearEvent(IndependentEvent):
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


class OneDimensionalNonCollaborativeMethylation(IndependentModel):
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
        model = BranchingDiffusion(r_b, r_d, r_um, r_mu)
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


class BranchingDiffusion(Model):
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

        methylation_evolution = lambda x : self._methylation_evolution(x, parameters)

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

    def run(self, parameters:dict, initial_state:dict, duration:float):
        
        
        
        r_max_b = self._r_b(minimize(lambda x: - self._r_b(x[0], parameters), np.array([1/2]), bounds = [(0, 1)]).x, parameters)
        r_max_d = self._r_d(minimize(lambda x: - self._r_d(x[0], parameters), np.array([1/2]), bounds = [(0, 1)]).x, parameters)
         
        r_max_cell = r_max_b + r_max_d
        
        methylation_evolution = lambda x, t : self._methylation_evolution(x, parameters)
        
        current_state = copy.deepcopy(initial_state)
        current_time = 0


        while True:
            cell_count = sum(current_state.values())
            max_rate = r_max_cell * cell_count
            if max_rate == 0:
                break
            else:
                waiting_time = - np.log(random.random()) / max_rate
            next_time = current_time + waiting_time


            if next_time > duration:
                times = [current_time, duration]
            else: 
                times = [current_time, next_time]
            new_state = {}
            for cell_type in current_state:
                next_key = odeint(methylation_evolution, cell_type, times)[1][0]
                new_state[next_key] = current_state[cell_type]
            current_state = new_state
            if next_time > duration:
                break
            current_time = next_time
            
            event_rng = random.random() * max_rate
            
            event_cell = None
            for cell_type, cell_count in current_state.items():
                if cell_count * r_max_cell >= event_rng:
                    event_cell = cell_type
                    event_rng = event_rng / cell_count
                    break
                event_rng -= cell_count * r_max_cell
            if event_cell is None:
                raise ValueError("event cell should always exist")
            
            

            for i, rate in enumerate([r_max_b, r_max_d]):
                if event_rng <= rate:
                    if i == 0:
                        if event_rng <= self._r_b(event_cell, parameters):
                            """BIRTH"""
                            current_state[event_cell] = current_state[event_cell] + 1            
                    elif i == 1:
                        if event_rng <= self._r_d(event_cell, parameters):
                            """DEATH"""
                            current_state[event_cell] = current_state[event_cell] - 1
                            if current_state[event_cell] == 0:
                                del current_state[event_cell] 
                    break 
                event_rng -= rate
            
            

        return current_state
    

    def _methylation_evolution(self, x, parameters):
            return - self._r_mu(x, parameters) * x + self._r_um(x, parameters) * (1 - x)
    


class BranchingDiffusion2(Model):
    def __init__(self, r_b, r_d, diffusion):
        """inputs are functions which take in (x:float, parameters:dict) and return rate:float"""
        self._r_b = r_b
        self._r_d = r_d
        self.diffusion = diffusion

    def calculate_extinction(self, parameters, point_count=1000):
        def extinction_derivative(y, x):
            b = self._r_b(x, parameters)
            d = self._r_d(x, parameters)
            diff = self.diffusion(x, parameters)
            
            if diff == 0:
                if y == 1:
                    return 0
                """Need to use L'hospital's rule to get a real answer"""
                return 0.1
            return ((d - b * y) * (1 - y)) / (diff)

        diffusion = lambda x : self.diffusion(x, parameters)

        initial_x = root(diffusion, 0.5).x[0]
        initial_y = min(self._r_d(initial_x, parameters) /
                        self._r_b(initial_x, parameters), 1)
        
        print(initial_x)
        print(initial_y)
        input()
        
        all_times = [i / (point_count - 1) for i in range(point_count)]
        low_times = [time for time in all_times if time <
                     initial_x] + [initial_x]
        high_times = [initial_x] + \
            [time for time in all_times if time >= initial_x]

        print(low_times)
        print(high_times)

        # get probabilities below the stable x
        low_res = odeint(extinction_derivative, initial_y,
                         list(reversed(low_times)))
        print(low_res)
        input()
        # get probabilities above the stable x
        high_res = odeint(extinction_derivative, initial_y, high_times)

        # combine probabilities into one list
        res = list(reversed(low_res[:, 0][1:])) + list(high_res[:, 0][1:])

        return res

    def run(self, parameters:dict, initial_state:dict, duration:float):
        
        
        
        r_max_b = self._r_b(minimize(lambda x: - self._r_b(x[0], parameters), np.array([1/2]), bounds = [(0, 1)]).x[0], parameters)
        r_max_d = self._r_d(minimize(lambda x: - self._r_d(x[0], parameters), np.array([1/2]), bounds = [(0, 1)]).x[0], parameters)
        
        r_max_cell = r_max_b + r_max_d
        
        def diffusion(x, t):
            return self.diffusion(x, parameters)
        
        current_state = copy.deepcopy(initial_state)
        current_time = 0


        while True:
            cell_count = sum(current_state.values())
            max_rate = r_max_cell * cell_count
            if max_rate == 0:
                # nothing else happens if nothing is alive
                break
            else:
                diffusing_time = - np.log(random.random()) / max_rate
            next_life_time = current_time + diffusing_time
            
            if next_life_time > duration:
                # diffuse until the end
                times = [current_time, duration]
            else: 
                # diffuse until the next birth/death
                times = [current_time, next_life_time]
            # solve diffusion
            new_state = {}
            for cell_type in current_state:
                next_key = odeint(diffusion, cell_type, times)[1][0]
                new_state[next_key] = current_state[cell_type]
            current_state = new_state
            
            # end if necessary, update the time otherwise
            if next_life_time > duration:
                break
            current_time = next_life_time
            
            # generate number to determine behavior
            event_rng = random.random() * max_rate
            
            # use number to determine cell
            event_cell = None
            for cell_type, cell_count in current_state.items():
                if cell_count * r_max_cell >= event_rng:
                    event_cell = cell_type
                    event_rng = event_rng / cell_count
                    break
                event_rng -= cell_count * r_max_cell
            
            # use number to determine if/what event occurrs
            for i, rate in enumerate([r_max_b, r_max_d]):
                if event_rng <= rate:
                    if i == 0:
                        if event_rng <= self._r_b(event_cell, parameters):
                            """BIRTH"""
                            current_state[event_cell] = current_state[event_cell] + 1            
                    elif i == 1:
                        if event_rng <= self._r_d(event_cell, parameters):
                            """DEATH"""
                            current_state[event_cell] = current_state[event_cell] - 1
                            if current_state[event_cell] == 0:
                                del current_state[event_cell] 
                    break 
                event_rng -= rate
            
            

        return current_state
    