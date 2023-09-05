"""
Models are implemented in a parameter-agnostic way. 
Parameters instead are supplied at the time of running.
"""
# pylint:disable=arguments-differ


class Model:
    """Abstract class. Contains a state space and function to run for a duration."""
    name = "Abstract Model"

    def run(self, parameters, initial_state, duration: float, **kwargs):
        """This is the only important function in a model. """
        # pylint:disable=unused-argument
        return initial_state

    def generate_simulation_data(self, parameters: dict, initial_state, timesteps: list,
                                 **kwargs):
        """
        Useful to run simulations. 
        Data is output in json-style:
        {
            "model": [model name],
            "parameters": [parameter dictionary],
            "data" [timepoint data dictionary],
        }
        """

        simulation_result = {
            "parameters": parameters,
            "model": self.name,
            "data": {0: initial_state},
        }
        last_time = 0
        current_state = initial_state
        for time in timesteps:
            duration = time - last_time
            current_state = self.run(
                parameters, current_state, duration, **kwargs)
            simulation_result["data"][time] = current_state
            last_time = time
        return simulation_result
