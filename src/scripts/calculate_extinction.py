from src.tools import model


b1 = model.Birth(0, "b1")
b2 = model.Birth(1, "b2")
d1 = model.Death(0, "d1")
d2 = model.Death(1, "d2")
t1 = model.Transition(0, 1, "0->1")
t2 = model.Transition(1, 0, "1->0")
parameters = {
    "b1": 2,
    "b2": 1,
    "d1": 1,
    "d2": 2,
    "0->1": 1,
    "1->0": 1,
}
events = [b1, b2, d1, d2, t1, t2]
M = model.LinearModel(events)
probabilities = M.calculate_extinction(parameters)
print(probabilities)
