from helper import json_load, json_write

majority_counter = {}
for round in range(10):
    reac_list = json_load(f"minimal_tcosa_reactions_round{round}.json")
    num_original = 0
    num_variant = 0
    for reac in reac_list:
        if "_ORIGINAL_" in reac:
            num_original += 1
        elif "_VARIANT_" in reac:
            num_variant += 1


    if num_original > num_variant:
        majority = "_ORIGINAL_"
    elif num_variant > num_original:
        majority = "_VARIANT_"
    else:
        print(f"No majority in round {round}")
        continue

    for reac in reac_list:
        base_id = reac.replace("_ORIGINAL_NADP_TCOSA", "")\
                      .replace("_ORIGINAL_NAD_TCOSA", "")\
                      .replace("_VARIANT_NAD_TCOSA", "")\
                      .replace("_VARIANT_NADP_TCOSA", "")
        if base_id not in majority_counter.keys():
            majority_counter[base_id] = {
                "majority": 0,
                "minority": 0,
            }
        if majority in reac:
            majority_counter[base_id]["majority"] += 1
        else:
            majority_counter[base_id]["minority"] += 1

json_write("./y.json", majority_counter)
