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

def get_blocked_reacs(data):
    blocked_reacs = []
    for reac_id in data.keys():
        if reac_id.startswith("x_") or reac_id.startswith("f_var_"):
            continue

        minf = data[reac_id]["min"]
        maxf = data[reac_id]["max"]

        if (minf <= 1e-5) and (maxf <= 1e-5):
            blocked_reacs.append(reac_id)
    return set(blocked_reacs)

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

for aerobicity in ("aerobic",): #, "anaerobic"):
    json_path_wildtype = f"cosa/results_aerobic/variability_TEST_0_568_swap_STANDARDCONCS_at_max_OPTSUBMDF_WILDTYPE.json"
    full_wildtype_data = json_load(
        json_path_wildtype
    )
    wildtype_data = full_wildtype_data["unswapped"]

    json_path_flexible = f"cosa/results_aerobic/variability_TEST_0_568_swap_STANDARDCONCS_at_max_OPTSUBMDF_FLEXIBLE.json"
    full_flexible_data = json_load(
        json_path_flexible
    )
    flexible_data = full_flexible_data["unswapped"]

    if "anaerobic" in aerobicity:
        print("\n=ANAEROBIC=")
    else:
        print("\n=AEROBIC=")

    try:
        optmdf_in_wildtype = wildtype_data["var_B"]["min"]
        optmdf_in_flexible = flexible_data["var_B"]["min"]
    except KeyError:
        pass
    try:
        optsubmdf_in_wildtype = wildtype_data["var_B2"]["min"]
        optsubmdf_in_flexible = flexible_data["var_B2"]["min"]
    except KeyError:
        pass

    essential_reacs_in_wildtype = get_essential_reacs(wildtype_data)
    essential_reacs_in_flexible = get_essential_reacs(flexible_data)

    blocked_reacs_in_wildtype = get_blocked_reacs(wildtype_data)
    blocked_reacs_in_flexible = get_blocked_reacs(flexible_data)

    # if "OPTMDF" in json_path:
    original_mdf = optmdf_in_flexible
    wildtype_mdf = optmdf_in_wildtype
    # else:
    #    original_mdf = optsubmdf_in_flexible
    #    wildtype_mdf = optsubmdf_in_wildtype
    bottleneck_reacs_in_wildtype = get_bottleneck_reacs(wildtype_data, essential_reacs_in_wildtype, wildtype_mdf)
    bottleneck_reacs_in_flexible = get_bottleneck_reacs(flexible_data, essential_reacs_in_flexible, original_mdf)

    essential_in_flexible_but_not_in_wildtype = essential_reacs_in_flexible.difference(essential_reacs_in_wildtype)
    essential_in_wildtype_but_not_in_flexible = essential_reacs_in_wildtype.difference(essential_reacs_in_flexible)

    blocked_in_flexible_but_not_in_wildtype = blocked_reacs_in_flexible.difference(blocked_reacs_in_wildtype)
    blocked_in_wildtype_but_not_in_flexible = blocked_reacs_in_wildtype.difference(blocked_reacs_in_flexible)

    # for swap_state in ("unswapped", "swapped"):
    #     if swap_state == "unswapped":
    #         name = "origina"
    #     else:
    #         print("")
    #         name = "swapped"
    #     print(f"concentration_ranges_{name} [M]")
    #     nad_min =   round(exp(full_data[swap_state]["x_nad_tcosa_c"]["min"]), 6)
    #     nad_max =   round(exp(full_data[swap_state]["x_nad_tcosa_c"]["max"]), 6)
    #     nadp_min =  round(exp(full_data[swap_state]["x_nadp_tcosa_c"]["min"]), 6)
    #     nadp_max =  round(exp(full_data[swap_state]["x_nadp_tcosa_c"]["max"]), 6)
    #     nadh_min =  round(exp(full_data[swap_state]["x_nadh_tcosa_c"]["min"]), 6)
    #     nadh_max =  round(exp(full_data[swap_state]["x_nadh_tcosa_c"]["max"]), 6)
    #     nadph_min = round(exp(full_data[swap_state]["x_nadph_tcosa_c"]["min"]), 6)
    #     nadph_max = round(exp(full_data[swap_state]["x_nadph_tcosa_c"]["max"]), 6)
    #     print(f"* NAD: {nad_min} to {nad_max}")
    #     print(f"* NADH: {nadh_min} to {nadh_max}")
    #     print(f"* NADP: {nadp_min} to {nadp_max}")
    #     print(f"* NADPH: {nadph_min} to {nadph_max}")

    print("\nessential_reacs_in_flexible")
    print(essential_reacs_in_flexible)
    print("\essential_reacs_in_wildtype")
    print(essential_reacs_in_wildtype)
    print(len(essential_reacs_in_wildtype))
    print([(x, wildtype_data[x]["min"], wildtype_data[x]["max"]) for x in essential_reacs_in_wildtype])

    # print("\nessential_in_flexible_but_not_in_wildtype:")
    # print_reacs(essential_in_flexible_but_not_in_wildtype)
    # print("\nessential_in_wildtype_but_not_in_flexible:")
    # print_reacs(essential_in_wildtype_but_not_in_flexible)

    print("\nblocked_in_flexible_but_not_in_wildtype:")
    print_reacs(blocked_in_flexible_but_not_in_wildtype)
    # print("\nblocked_in_wildtype_but_not_in_flexible:")
    # print_reacs(blocked_in_wildtype_but_not_in_flexible)

    print("\nbottleneck_reacs_in_wildtype")
    print_reacs(bottleneck_reacs_in_wildtype)
    print("\nbottleneck_reacs_in_flexible")
    print_reacs(bottleneck_reacs_in_flexible)
