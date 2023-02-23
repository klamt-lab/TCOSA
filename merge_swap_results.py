import copy
from helper import json_load, json_write

lower = json_load("cosa/results_aerobic/swap_results_STANDARDCONC_lower.json")
upper = json_load("cosa/results_aerobic/swap_results_STANDARDCONC_upper.json")

combined = {}
for reac_id in upper.keys():
    combined[reac_id] = {}
    for mu in upper[reac_id].keys():
        combined[reac_id][mu] = copy.deepcopy(upper[reac_id][mu])
    for mu in lower[reac_id].keys():
        combined[reac_id][mu] = copy.deepcopy(lower[reac_id][mu])
json_write("cosa/results_aerobic/swap_results_STANDARDCONC.json", combined)
