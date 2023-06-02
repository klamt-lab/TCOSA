from math import log, log10
from helper import json_load, json_write


measurementjson = json_load("resources/in_vivo_concentration_data/2009_Bennet_full_data.json")
for target in "OPTMDF", "OPTSUBMDF":
    outjson = {}
    cvajson = json_load(f"cosa/results_aerobic/cva_{target}_STANDARDCONC.json")
    comparisonjson = json_load("cosa/cva_comparison.json")

    defective_ids = []
    for metabolite_id in comparisonjson.keys():
        if "0,818" in comparisonjson[metabolite_id][target]["defective_mus"]:
            defective_ids.append(metabolite_id)

    for defective_id in defective_ids:
        data = cvajson["x_"+defective_id]["0,818"]
        min_cva = data["min"]
        max_cva = data["max"]

        min_measurement = measurementjson[defective_id]["lb"]
        max_measurement = measurementjson[defective_id]["ub"]
        name_measurement = measurementjson[defective_id]["name"]

        min_difference = min(
            abs(min_measurement - min_cva),
            abs(min_measurement - max_cva),
            abs(max_measurement - min_cva),
            abs(max_measurement - max_cva),
        )

        outjson[defective_id] = {
            "name": name_measurement,
            "min_cva": min_cva,
            "max_cva": max_cva,

            "min_cva_ln": log(min_cva),
            "max_cva_ln": log(max_cva),
            "min_cva_log10": log10(min_cva),
            "max_cva_log10": log10(max_cva),

            "min_measurement": min_measurement,
            "max_measurement": max_measurement,
            "min_difference": min_difference,

            "min_measurement_ln": log(min_measurement),
            "max_measurement_ln": log(max_measurement),
            "min_difference_ln":  log(min_difference),
            "min_measurement_log10": log10(min_measurement),
            "max_measurement_log10": log10(max_measurement),
            "min_difference_log10": log10(min_difference),
        }

    json_write("./cosa/results_aerobic/cva_metastatistics_"+target+".json", outjson)
