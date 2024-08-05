from copy import deepcopy
from itertools import product
from math import log
from random import random

from scipy.integrate import odeint, ode

from src.tools.models.population import PopulationModel
from src.tools.models.rng import generate_poisson


class PDMP(PopulationModel):
    site_names = {0: "m", 1: "h", 2: "u"}
    def __init__(self): 
        super().__init__(3)

    def run(self, parameters, initial_state, duration):
        state = deepcopy(initial_state)
        last_diffusion = False
        
        """
        def flow(x, t):
            flow = [0, 0, 0]
            for a, b, c in product(range(3), range(3), range(4)):
                if c == 3: 
                    noncollaborative_rate = self._get_parameter(parameters, a, b)
                    flow[a] = flow[a] - x[a] * noncollaborative_rate
                    flow[b] = flow[b] + x[a] * noncollaborative_rate
                else:
                    collaborative_rate = self._get_parameter(parameters, a, b, c)
                    flow[a] = flow[a] - x[a] * x[c] * collaborative_rate
                    flow[b] = flow[b] + x[a] * x[c] * collaborative_rate
            return flow
        """  

        while True:            
            b = parameters["b"]
            waiting_time = - log(random()) / b
            if waiting_time > duration:
                waiting_time = duration
                last_diffusion = True
            flow = self.flow(parameters)
            state = odeint(flow, state, [0, waiting_time], tfirst=True)[1]
            if last_diffusion:
                return state
            state = [0, state[0] + state[1] / 2, state[1] / 2 + state[2]]

    def generate_timepoint_data(self, parameters, initial_state, times, splitting_times = None):
        result = {0 : list(initial_state)}
        time = 0
        state = deepcopy(initial_state)
        if splitting_times is None:
            splitting_times = generate_poisson(parameters["b"], times[-1])
        flow = self._make_flow_func(parameters)
        
        all_times = times + splitting_times
        all_times.sort()

        chunk = []
        for time in all_times:
            chunk.append(time)
            if time in splitting_times:
                data = odeint(flow, state, chunk)
                for t, d in zip(chunk, data):
                    result[t] = list(d)
                presplit = result[time]
                postsplit = [0, presplit[0] + presplit[1] / 2, presplit[1] / 2 + presplit[2]]
                result[time] = list(postsplit)
                state = postsplit
                chunk = [time]
        data = odeint(flow, state, chunk)
        for t, d in zip(chunk, data):
            result[t] = list(d)
        return result

    def sample_simulataneously(self, parameters, initial_states, times):
        splitting_times = generate_poisson(parameters["b"], times[-1])
        data = []
        for initial_state in initial_states:
            data.append(self.generate_timepoint_data(parameters, initial_state, times, splitting_times=splitting_times))
        return data

    def sample_wasserstein(self, parameters, times, sample_count):
        if not self.verify_wasserstein_lemma(parameters):
            res = input("Warning: coupling may not be optimal. Continue? y/n: ")
            if res == "n":
                return
        
        first_initial = [1, 0, 0]
        second_initial = [0, 0, 1]
        initial_states = [first_initial, second_initial]
        total_distances = {time: 0 for time in times}
        for _ in range(sample_count):
            print(_)
            result = self.sample_simulataneously(parameters, initial_states, times)
            first_result = result[0]
            second_result = result[1]
            for time in times:
                first_data = first_result[time]
                second_data = second_result[time]
                distance = first_data[0] - second_data[0] + second_data[2] - first_data[2]
                total_distances[time] = total_distances[time] + distance
        average_distances = {}
        for time in times:
            average_distances[time] = total_distances[time] / sample_count
        return average_distances
        
    def get_hitting_time(self, parameters, start, end, stepsize = 0.01):
        if not self.verify_wasserstein_lemma(parameters):
            raise NotImplementedError
        initial_condition = [0, start, 1 - start]
        r = ode(self.flow(parameters))
        r.set_initial_value(initial_condition)
        while r.successful() and r.y[0] * 2 + r.y[1] < end:
            r.integrate(r.t + stepsize)
        return r.t




    @classmethod
    def flow(cls, parameters) -> list[float]:
        def flow_func(t, x):
            flow_at_point = [0, 0, 0]
            for a, b, c in product(range(3), range(3), range(4)):
                if c == 3: 
                    noncollaborative_rate = cls._get_parameter(parameters, a, b)
                    flow_at_point[a] = flow_at_point[a] - x[a] * noncollaborative_rate
                    flow_at_point[b] = flow_at_point[b] + x[a] * noncollaborative_rate
                else:
                    collaborative_rate = cls._get_parameter(parameters, a, b, c)
                    flow_at_point[a] = flow_at_point[a] - x[a] * x[c] * collaborative_rate
                    flow_at_point[b] = flow_at_point[b] + x[a] * x[c] * collaborative_rate
            return flow_at_point
        return flow_func
    
    @staticmethod
    def verify_wasserstein_lemma(parameters):
        if parameters["r_hu_u"] < parameters["r_uh_h"]:
            return False
        if parameters["r_hm_m"] < parameters["r_mh_h"]:
            return False
        return True

    @staticmethod
    def convert_to_cytosine_mean(parameters):
        p2 = {"b": parameters["b"]}
        ave = (parameters["r_uh"] / 2 + parameters["r_hm"]) / 2
        p2["r_hm"] = ave
        p2["r_uh"] = ave * 2
        ave = (parameters["r_mh"] / 2 + parameters["r_hu"]) / 2
        p2["r_hu"] = ave
        p2["r_mh"] = ave * 2
        ave = (parameters["r_uh_h"] / 2 + parameters["r_uh_m"] / 4 + parameters["r_hm_h"] + parameters["r_hm_m"] / 2) / 4
        p2["r_hm_h"] = ave
        p2["r_hm_m"] = ave * 2
        p2["r_uh_h"] = ave * 2
        p2["r_uh_m"] = ave * 4
        ave = (parameters["r_mh_h"] / 2 + parameters["r_mh_u"] / 4 + parameters["r_hu_h"] + parameters["r_hu_u"] / 2) / 4
        p2["r_hu_h"] = ave
        p2["r_hu_u"] = ave * 2
        p2["r_mh_h"] = ave * 2
        p2["r_mh_u"] = ave * 4
        return p2
    
    @staticmethod
    def convert_to_noncollaborative(parameters):
        p2 = deepcopy(parameters)        
        p2["r_hm_h"] = 0
        p2["r_hm_m"] = 0
        p2["r_uh_h"] = 0
        p2["r_uh_m"] = 0        
        p2["r_hu_h"] = 0
        p2["r_hu_u"] = 0
        p2["r_mh_h"] = 0
        p2["r_mh_u"] = 0
        return p2


    def _make_flow_func(self, parameters): 
        def flow(x, t):
            flow = [0, 0, 0]
            for a, b, c in product(range(3), range(3), range(4)):
                if c == 3: 
                    noncollaborative_rate = self._get_parameter(parameters, a, b)
                    flow[a] = flow[a] - x[a] * noncollaborative_rate
                    flow[b] = flow[b] + x[a] * noncollaborative_rate
                else:
                    collaborative_rate = self._get_parameter(parameters, a, b, c)
                    flow[a] = flow[a] - x[a] * x[c] * collaborative_rate
                    flow[b] = flow[b] + x[a] * x[c] * collaborative_rate
            return flow
        return flow
    
    @classmethod
    def _get_parameter(cls, parameters, a, b, c = None):
            """gets the parameter r_ab_c for a -> b mediated by c"""
            if (a - b) ** 2 != 1:
                return 0
            if a > b and c == 2:
                return 0
            if a < b and c == 0:
                return 0
            name = "r_" + cls.site_names[a] + cls.site_names[b]
            if c is not None:
                name = name + "_" + cls.site_names[c]
            return parameters[name]
