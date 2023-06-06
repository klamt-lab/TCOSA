from math import log, log10
from helper import json_load, json_write


measurementjson = json_load("resources/in_vivo_concentration_data/2009_Bennet_full_data.json")
measurementjson["nad_tcosa_c"] = {
    "lb": 2.32e-3,
    "ub": 2.8e-3,
    "name": "NAD",
}
measurementjson["nadh_tcosa_c"] = {
    "lb": 5.45e-5,
    "ub": 1.27e-4,
    "name": "NADH",
}
measurementjson["nadp_tcosa_c"] = {
    "lb": 1.4e-7,
    "ub": 3.11e-5,
    "name": "NADP",
}
measurementjson["nadph_tcosa_c"] = {
    "lb": 1.1e-4,
    "ub": 1.34e-4,
    "name": "NADPH",
}

outjson = {}
for target in "OPTMDF", "OPTSUBMDF":
    outjson[target] = {}

    cvajson = json_load(f"cosa/results_aerobic/cva_{target}_STANDARDCONC.json")
    comparisonjson = json_load("cosa/cva_comparison.json")

    defective_ids = []
    for metabolite_id in comparisonjson.keys():
        if "0,818" in comparisonjson[metabolite_id][target]["defective_mus"]:
            defective_ids.append(metabolite_id)

    outjson[target]["max_absolute_min_difference"] = 0.0
    outjson[target]["max_abs_min_log10_difference"] = 0.0
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

        min_log10_difference = min(
            log10(min_measurement) - log10(min_cva),
            log10(min_measurement) - log10(max_cva),
            log10(max_measurement) - log10(min_cva),
            log10(max_measurement) - log10(max_cva),
        )

        outjson[target][defective_id] = {
            "name": name_measurement,
            "min_cva": min_cva,
            "max_cva": max_cva,

            "min_measurement": min_measurement,
            "max_measurement": max_measurement,
            "min_difference": min_difference,

            "min_log10_difference": min_log10_difference,
        }

        outjson[target]["max_absolute_min_difference"] = max(outjson[target]["max_absolute_min_difference"], min_difference)
        outjson[target]["max_abs_min_log10_difference"] = max(outjson[target]["max_abs_min_log10_difference"], min_log10_difference)

json_write("./cosa/results_aerobic/cva_metastatistics.json", outjson)
