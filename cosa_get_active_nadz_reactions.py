from math import exp
from helper import json_zip_load, get_files

files = [
    "cosa/results_anaerobic_expanded_nadz_change_30/runs/OPTMDF_STANDARDCONC_FLEXIBLE.json",
    "cosa/results_anaerobic_expanded_nadz_change_30/runs/OPTSUBMDF_STANDARDCONC_FLEXIBLE.json",
]

for file in files:
    data = json_zip_load(file)

    print("===")
    print("FILE:", file)

    for growth_rate in data.keys():
        print("~")
        print("GROWTH RATE:", growth_rate)
        for reac_id in data[growth_rate]["values"].keys():
            if reac_id.startswith("f_var_") or reac_id.startswith("dG0_") or reac_id.startswith("var_z_") or reac_id.startswith("z_var_"):
                continue
            if "NADZ" not in reac_id:
                continue
            reac_v = data[growth_rate]["values"][reac_id]
            if abs(reac_v) < 1e-7:
                continue
            print(reac_id, reac_v)
        print("|")
        print("[NAD]:", exp(data[growth_rate]["values"]["x_nad_tcosa_c"]))
        print("[NADH]:", exp(data[growth_rate]["values"]["x_nadh_tcosa_c"]))

        print("[NADP]:", exp(data[growth_rate]["values"]["x_nadp_tcosa_c"]))
        print("[NADPH]:", exp(data[growth_rate]["values"]["x_nadph_tcosa_c"]))

        print("[NADZ]:", exp(data[growth_rate]["values"]["x_nadz_tcosa_c"]))
        print("[NADZH]:", exp(data[growth_rate]["values"]["x_nadzh_tcosa_c"]))

        print("|")
        print("[NAD]/[NADH]:", exp(data[growth_rate]["values"]["x_nad_tcosa_c"]) / exp(data[growth_rate]["values"]["x_nadh_tcosa_c"]))
        print("[NADP]/[NADPH]:", exp(data[growth_rate]["values"]["x_nadp_tcosa_c"]) / exp(data[growth_rate]["values"]["x_nadph_tcosa_c"]))
        print("[NADZ]/[NADZH]:", exp(data[growth_rate]["values"]["x_nadz_tcosa_c"]) / exp(data[growth_rate]["values"]["x_nadzh_tcosa_c"]))
