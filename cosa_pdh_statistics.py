import json
from tokenize import cookie_re
from helper import json_load
import cobra
from math import exp

dG0_values = json_load("cosa/dG0_values.json")
cobra_model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")

def print_reacs(reacs):
    for reac in reacs:
        cobra_reac = cobra_model.reactions.get_by_id(reac)
        dG0 = round(dG0_values[reac]["dG0"], 2)
        print(f"* {reac};  ΔG0={dG0} kJ/mol; {cobra_reac.reaction}")

def get_essential_reacs(data):
    active_reacs = []
    for reac_id in data.keys():
        if reac_id.startswith("x_") or reac_id.startswith("f_var_"):
            continue

        minf = data[reac_id]["min"]
        maxf = data[reac_id]["max"]

        if (minf > 0.0) and (maxf > 0.0):
            active_reacs.append(reac_id)
    return set(active_reacs)

def get_bottleneck_reacs(data, essential_reacs, optsubmdf):
    optsubmdf += 1e-6
    bottleneck_reacs = []
    for f_var_id in data.keys():
        if not f_var_id.startswith("f_var_"):
            continue

        if f_var_id.replace("f_var_", "") not in essential_reacs:
            continue

        minf = data[f_var_id]["min"]
        maxf = data[f_var_id]["max"]

        if (maxf < optsubmdf):
            bottleneck_reacs.append(f_var_id.replace("f_var_", ""))
    return set(bottleneck_reacs)

for aerobicity in ("aerobic", "anaerobic"):
    json_path = f"cosa/results_{aerobicity}/variability_PDH_ORIGINAL_NAD_TCOSA_swap_STANDARDCONCS_at_max_OPTSUBMDF.json"
    full_data = json_load(
        json_path
    )
    swapped_data = full_data["swapped"]
    origina_data = full_data["unswapped"]

    if "anaerobic" in json_path:
        print("\n=ANAEROBIC=")
    else:
        print("\n=AEROBIC=")

    try:
        optmdf_in_swapped = swapped_data["var_B"]["min"]
        optmdf_in_origina = origina_data["var_B"]["min"]
    except KeyError:
        pass
    try:
        optsubmdf_in_swapped = swapped_data["var_B2"]["min"]
        optsubmdf_in_origina = origina_data["var_B2"]["min"]
    except KeyError:
        pass

    essential_reacs_in_swapped = get_essential_reacs(swapped_data)
    essential_reacs_in_origina = get_essential_reacs(origina_data)

    if "OPTMDF" in json_path:
        original_mdf = optmdf_in_origina
        swapped_mdf = optmdf_in_swapped
    else:
        original_mdf = optsubmdf_in_origina
        swapped_mdf = optsubmdf_in_swapped
    bottleneck_reacs_in_swapped = get_bottleneck_reacs(swapped_data, essential_reacs_in_swapped, swapped_mdf)
    bottleneck_reacs_in_origina = get_bottleneck_reacs(origina_data, essential_reacs_in_origina, original_mdf)

    essential_in_origina_but_not_in_swapped = essential_reacs_in_origina.difference(essential_reacs_in_swapped)
    essential_in_swapped_but_not_in_origina = essential_reacs_in_swapped.difference(essential_reacs_in_origina)

    for swap_state in ("unswapped", "swapped"):
        if swap_state == "unswapped":
            name = "origina"
        else:
            print("")
            name = "swapped"
        print(f"concentration_ranges_{name} [M]")
        nad_min =   round(exp(full_data[swap_state]["x_nad_tcosa_c"]["min"]), 6)
        nad_max =   round(exp(full_data[swap_state]["x_nad_tcosa_c"]["max"]), 6)
        nadp_min =  round(exp(full_data[swap_state]["x_nadp_tcosa_c"]["min"]), 6)
        nadp_max =  round(exp(full_data[swap_state]["x_nadp_tcosa_c"]["max"]), 6)
        nadh_min =  round(exp(full_data[swap_state]["x_nadh_tcosa_c"]["min"]), 6)
        nadh_max =  round(exp(full_data[swap_state]["x_nadh_tcosa_c"]["max"]), 6)
        nadph_min = round(exp(full_data[swap_state]["x_nadph_tcosa_c"]["min"]), 6)
        nadph_max = round(exp(full_data[swap_state]["x_nadph_tcosa_c"]["max"]), 6)
        print(f"* NAD: {nad_min} to {nad_max}")
        print(f"* NADH: {nadh_min} to {nadh_max}")
        print(f"* NADP: {nadp_min} to {nadp_max}")
        print(f"* NADPH: {nadph_min} to {nadph_max}")

    print("\nessential_in_origina_but_not_in_swapped:")
    print_reacs(essential_in_origina_but_not_in_swapped)
    print("\nessential_in_swapped_but_not_in_origina:")
    print_reacs(essential_in_swapped_but_not_in_origina)


    print("\nbottleneck_reacs_in_swapped")
    print_reacs(bottleneck_reacs_in_swapped)
    print("\nbottleneck_reacs_in_origina")
    print_reacs(bottleneck_reacs_in_origina)
