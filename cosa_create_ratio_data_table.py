from helper import json_load

jsons = [
    "cosa/results_aerobic/ratio_ratio_test_data.json",
    "cosa/results_aerobic_acetate/ratio_ratio_test_data.json",
    "cosa/results_anaerobic/ratio_ratio_test_data.json",
]
folders = [
    "cosa/results_aerobic/",
    "cosa/results_aerobic_acetate/",
    "cosa/results_anaerobic/",
]

header = "µ [1/h];Min [NADH]/[NAD] / [NADPH] / [NADP] ratio ratio;Max [NADH]/[NAD] / [NADPH] / [NADP] ratio\n"
for i in range(len(jsons)):
    json = jsons[i]
    folder = folders[i]

    data = json_load(json)
    for key in data.keys():
        if "_OPTMDF_" in key:
            target = "MDF"
        else:
            target = "SubMDF"
        if "STANDARDCONC" in key:
            concentrations = "standard_concentrations"
        else:
            concentrations = "paper_concentrations"

        body = ""
        for entry_i in range(len(data[key]["plotted_growth_rates"])):
            mu = data[key]["plotted_growth_rates"][entry_i]
            min_ratio = data[key]["min_ratios"][entry_i]
            max_ratio = data[key]["max_ratios"][entry_i]
            body += f"{mu};{min_ratio};{max_ratio}\n"
        with open(folder+"ratio_ratios_table_"+target+"_"+concentrations+".csv", "w", encoding="utf-8") as f:
            f.write(header+body)
