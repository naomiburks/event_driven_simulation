"""
This module contains tools for modeling methylation. 
In particular, we instantiate:
    - Methylation and Demethylation events, where the rate varies linearly with the population index
    - Interpolated linear events, where we can establish linearly varying birth and death rates
    to simulate a fitness landscape.
"""

from src.tools import model


class Methylation(model.Transition):
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


class Demethylation(model.Transition):
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


class InterpolatedLinearEvent(model.LinearEvent):
    """
    The rate of an interpolated linear event is calculated 
    linearly between the rate for pop 0 and pop max.
    """

    def __init__(self, population_index, population_count, \
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


class InterpolatedBirth(InterpolatedLinearEvent, model.Birth):
    """
    Interpolated birth is both a Birth and an InterpolatedLinearEvent
    """


class InterpolatedDeath(InterpolatedLinearEvent, model.Death):
    """
    Interpolated death is both a Death and an InterpolatedLinearEvent
    """


class OneDimensionalNonCollaborativeMethylation(model.LinearModel):
    """
    Defines a model with the following properties:
        - M + 1 types: type for each of 0 sites methylated through M sites methylated
        - Methylation (x -> x+1) and Demthyation (x -> x-1) events are noncollaborative. 
            - Rates proportional to the number of unmethylated (or methylated) sites.
        - Birth and death rates follow vary linearly with methylation level. 
    """
    name = "One Dimensional Noncollaborative Methylation"
    def __init__(self, M: int):
        events = []
        for i in range(M + 1):
            if i != M: # add methylations
                events.append(Methylation(i, M, "r_um"))
            if i != 0: # add demethylations
                events.append(Demethylation(i, "r_mu"))
            # add births and deaths
            events.append(InterpolatedBirth(i, M + 1, "b_0", "b_M"))
            events.append(InterpolatedDeath(i, M + 1, "d_0", "d_M"))
        super().__init__(events)
        self.name = f"{self.name} ({self.population_count - 1} sites)"