"""This script generates a comparison of Concentration Variability Analysis results with the measured concentrations from Bennett et al., 2009.

The result is stored in a JSON file named cva_comparison.json.
"""

from helper import json_load, json_write


paper_data = json_load("resources/in_vivo_concentration_data/final_concentration_values_paper.json")
comparison = {}
for target in ("OPTMDF", "OPTSUBMDF"):
    cva_data = json_load(f"cosa/results_aerobic/cva_{target}_STANDARDCONC.json")

    for metabolite_id in paper_data.keys():
        if metabolite_id in ["h2o_c", "h2o_p", "h2o_e", "h_c", "h_p", "h_e"]:
            continue

        if "x_"+metabolite_id not in cva_data.keys():
            print(f"INFO: Metabolite {metabolite_id} missing in CVA")
            continue

        if metabolite_id not in comparison.keys():
            comparison[metabolite_id] = {}
        if target not in comparison[metabolite_id].keys():
            comparison[metabolite_id][target] = {}

        is_fully_inside_everywhere = False
        is_fully_inside_at_least_once = False

        is_partly_inside_everywhere = True
        is_partly_inside_at_least_once = False

        defective_mus = []
        min_differences = []
        for growth_rate in cva_data["x_"+metabolite_id].keys():
            min_paper = paper_data[metabolite_id]["min"]
            min_cva = cva_data["x_"+metabolite_id][growth_rate]["min"]

            max_paper = paper_data[metabolite_id]["max"]
            max_cva = cva_data["x_"+metabolite_id][growth_rate]["max"]


            if (min_cva >= min_paper) and (max_cva <= max_paper):
                is_fully_inside_at_least_once = True
            else:
                is_fully_inside_everywhere = False


            if ((min_cva >= min_paper) and (max_cva <= max_paper)) or ((min_cva <= min_paper) and (max_cva >= max_paper)) or ((min_cva <= min_paper) and (max_cva >= min_paper)) or ((min_paper <= max_paper) and (max_cva >= max_paper)):
                is_partly_inside_at_least_once = True
            else:
                is_partly_inside_everywhere = False
                defective_mus.append(growth_rate)
                min_differences.append(min(abs(min_cva-min_paper), abs(min_cva-max_paper), abs(max_cva-max_paper), abs(max_cva-min_paper)))


        # comparison[metabolite_id][target]["is_fully_inside_everywhere"] = is_fully_inside_everywhere
        # comparison[metabolite_id][target]["is_fully_inside_at_least_once"] = is_fully_inside_at_least_once

        comparison[metabolite_id][target]["is_partly_inside_everywhere"] = is_partly_inside_everywhere
        comparison[metabolite_id][target]["is_partly_inside_at_least_once"] = is_partly_inside_at_least_once

        comparison[metabolite_id][target]["defective_mus"] = defective_mus
        comparison[metabolite_id][target]["min_differences"] = min_differences

        if len(defective_mus) > 1:
            print(f"INFO: {metabolite_id} has defective µ with max min difference of {max(min_differences)} M!")


json_write("./cosa/cva_comparison.json", comparison)
