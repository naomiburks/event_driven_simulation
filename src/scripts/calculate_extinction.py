from src.tools.models.homogeneous import Birth, Death, Transition, HomogeneousModel

b1 = Birth(0, "b1")
b2 = Birth(1, "b2")
d1 = Death(0, "d1")
d2 = Death(1, "d2")
t2 = Transition(1, 0, "1->0")
t1 = Transition(0, 1, "0->1")
parameters = {
    "b1": 2,
    "b2": 1,
    "d1": 1,
    "d2": 2,
    "0->1": 1,
    "1->0": 1,
}
events = [b1, b2, d1, d2, t1, t2]
M = HomogeneousModel(events)
probabilities = M.calculate_extinction(parameters)
print(probabilities)
