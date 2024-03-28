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
from src.tools.models.homogeneous import Birth, Death, IndependentModel, Switch
from src.tools.models.model import Model


class OneDimensionalNonCollaborativeMethylation(Switch):
    """Describes a methylation event"""

    def __init__(self, site_index, site_count):    
        self.site_count = site_count
        def get_rate_from_parameters(parameters):
            return (self.site_count - self.population_index) * parameters["r_um"]

        super().__init__(site_index, site_index + 1, get_rate_from_parameters)
        

class OneDimensionalNonCollaborativeDemethylation(Switch):
    """Describes a demethylation event"""

    def __init__(self, site_index, site_count):
        self.site_count = site_count
        def get_rate_from_parameters(parameters):
            return self.population_index * parameters["r_mu"]
        
        super().__init__(site_index, site_index - 1, get_rate_from_parameters)


class OneDimensionalBirth(Birth):
    def __init__(self, site_index, site_count):
        self.site_count = site_count
        def get_rate_from_parameters(params):
            return (params["b_0"] * (self.site_count - self.population_index) + params["b_M"] * self.population_index) / self.site_count
        super().__init__(site_index, get_rate_from_parameters)


class OneDimensionalDeath(Death):
    def __init__(self, site_index, site_count):
        self.site_count = site_count
        def get_rate_from_parameters(params):
            return (params["d_0"] * (self.site_count - self.population_index) + params["d_M"] * self.population_index) / self.site_count
        super().__init__(site_index, get_rate_from_parameters)


class OneDimensionalNonCollaborative(IndependentModel):
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
                events.append(OneDimensionalNonCollaborativeMethylation(i, M))
            if i != 0:  # add demethylations
                events.append(OneDimensionalNonCollaborativeDemethylation(i, M))
            # add births and deaths
            events.append(OneDimensionalBirth(i, M))
            events.append(OneDimensionalDeath(i, M))
        super().__init__(events)
        self.name = f"{self.name} ({self.population_count - 1} sites)"

    @staticmethod
    def get_limit_model():
        def r_um(x, parameters):
            return parameters["r_um"]

        def r_mu(x, parameters):
            return parameters["r_mu"]

        def r_b(x, parameters):
            return x * parameters["b_M"] + (1 - x) * parameters["b_0"]

        def r_d(x, parameters):
            return x * parameters["d_M"] + (1 - x) * parameters["d_0"]
        model = InfiniteSiteOneDimensional(r_b, r_d, r_um, r_mu)
        return model


class OneDimensionalCollaborativeMethylation(Switch):
    """Describes a methylation event"""

    def __init__(self, site_index, site_count):

        self.site_count = site_count
        def get_rate_from_parameters(parameters):
            M = self.site_count
            x = self.population_index
            r_uh = parameters["r_uh"]
            r_uh_m = parameters["r_uh_m"]
            r_hm = parameters["r_hm"]
            r_hm_h = parameters["r_hm_h"]
            r_hm_m = parameters["r_hm_m"]
            r_hu = parameters["r_hu"]
            r_hu_h = parameters["r_hu_h"]
            r_hu_u = parameters["r_hu_u"]
            numerator = (M - x) * (r_uh + r_uh_m * x) * (r_hm + r_hm_h + r_hm_m * x)
            denominator = r_hu + r_hu_h + r_hu_u * (M - x - 1) + r_hm + r_hm_h + r_hm_m * (x - 1)
            return numerator / denominator

        super().__init__(site_index, site_index + 1, get_rate_from_parameters)
        

class OneDimensionalCollaborativeDemethylation(Switch):
    """Describes a demethylation event"""

    def __init__(self, site_index, site_count):
        self.site_count = site_count
        def get_rate_from_parameters(model_parameters):
            M = self.site_count
            x = self.population_index
            r_mh = model_parameters["r_mh"]
            r_mh_u = model_parameters["r_mh_u"]
            r_hm = model_parameters["r_hm"]
            r_hm_h = model_parameters["r_hm_h"]
            r_hm_m = model_parameters["r_hm_m"]
            r_hu = model_parameters["r_hu"]
            r_hu_h = model_parameters["r_hu_h"]
            r_hu_u = model_parameters["r_hu_u"]
            
            numerator = x * (r_mh + r_mh_u * (M - x)) * (r_hu + r_hu_h + r_hu_u * (M - x))
            denominator = r_hu + r_hu_h + r_hu_u * (M - x) + r_hm + r_hm_h + r_hm_m * (x - 1)
            return numerator / denominator


        super().__init__(site_index, site_index - 1, get_rate_from_parameters)
        


class OneDimensionalCollaborative(IndependentModel):
    """
    Defines the collaborative version of the 1D simplified model. Properties:
        - M + 1 types: type for each of 0 sites methylated through M sites methylated
        - Methylation (x -> x+1) and Demthyation (x -> x-1) events are collaborative. 
        - Birth and death rates vary linearly with methylation level. 
    """
    name = "One Dimensional Collaborative Methylation with Linear Fitness"

    def __init__(self, M: int):
        events = []
        for i in range(M + 1):
            if i != M:  # add methylations
                events.append(OneDimensionalCollaborativeMethylation(i, M))
            if i != 0:  # add demethylations
                events.append(OneDimensionalCollaborativeDemethylation(i, M))
            # add births and deaths
            events.append(OneDimensionalBirth(i, M))
            events.append(OneDimensionalDeath(i, M))
        super().__init__(events)
        self.name = f"{self.name} ({self.population_count - 1} sites)"


class SingleSite(ConstantEventModel):
    def __init__(self):
        events = []
        events.append(ConstantEvent("u", "h", lambda x : x["r_uh"]))
        events.append(ConstantEvent("h", "m", lambda x : x["r_hm"]))
        events.append(ConstantEvent("h", "u", lambda x : x["r_hu"]))
        events.append(ConstantEvent("m", "h", lambda x : x["r_mh"]))
        events.append(ConstantEvent("m", "u", lambda x : x["b"]))
        events.append(ConstantEvent("h", "u", lambda x : x["b"] / 2))
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
    def __init__(self, r_b, r_d, diffusion):
        """inputs are functions which take in (x:float, parameters:dict) and return rate:float"""
        """x should be normalized to be between 0 and 1"""
        self._r_b = r_b # birth rate
        self._r_d = r_d # death rate
        self.diffusion = diffusion # change in x

    def calculate_extinction(self, parameters, point_count=1001):
        def extinction_derivative(y, x):
            b = self._r_b(x, parameters)
            d = self._r_d(x, parameters)
            diff = self.diffusion(x, parameters)
            
            if diff == 0:
                if y == 1:
                    return 0
                """Need to use L'hospital's rule to get a real answer"""
                return 0
            return - (d - b * y) * (1 - y) / (diff)

        def diffusion(x):
            return self.diffusion(x, parameters)

        initial_x = root(diffusion, 0.5).x[0]
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
        r_max_b = self._r_b(minimize(lambda x: - self._r_b(x[0], parameters), np.array([1/2]), bounds = [(0, 1)]).x[0], parameters)
        r_max_d = self._r_d(minimize(lambda x: - self._r_d(x[0], parameters), np.array([1/2]), bounds = [(0, 1)]).x[0], parameters)
        
        r_max_cell = r_max_b + r_max_d
        
        def diffusion(x, t) -> float:
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
    



class InfiniteSiteOneDimensional(BranchingDiffusion):
    def __init__(self, r_b, r_d, r_um, r_mu):
        def diffusion(x, parameters):
            return r_um(x, parameters) * (1 - x) - r_mu(x, parameters) * x
        super().__init__(r_b, r_d, diffusion)
   
    