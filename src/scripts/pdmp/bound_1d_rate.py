from src.tools.models.pdmp import PDMP
from src.constants import BISTABLE_STRONG_PDMP_PARAMS

flow = PDMP.flow(BISTABLE_STRONG_PDMP_PARAMS)

level = 1

min_hemi = 0
max_hemi = min(level, 1)

def get_growth(h):
    state = [1 - h / 2 - level / 2, h, (level - h) / 2]
    change = flow(0, state)
    