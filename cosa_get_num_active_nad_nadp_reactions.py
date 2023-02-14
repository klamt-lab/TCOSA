from helper import json_zip_load

results = json_zip_load("cosa/results_aerobic/runs/OPTSUBMDF_STANDARDCONC_FLEXIBLE.json")

for growth_rate in results.keys():
    num_counted = 0
    values = results[growth_rate]["values"]
    for reaction_id in values:
        if not reaction_id.endswith("_TCOSA"):
            continue
        if reaction_id.startswith("f_var"):
            continue
        if reaction_id.startswith("z_var"):
            continue
        if reaction_id.startswith("dG0"):
            continue
        flux = results[growth_rate]["values"][reaction_id]
        if flux > 1e-7:
            print(f"Active reaction {reaction_id} with flux {flux}")
            num_counted += 1
    print(f"Num at growth rate {growth_rate}:", num_counted)
    input()
