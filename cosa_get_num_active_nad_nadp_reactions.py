from helper import json_zip_load, get_files

# json_paths = ["./cosa/results_aerobic/runs/"+x for x in get_files("./cosa/results_aerobic/runs")] +\
#              ["./cosa/results_anaerobic/runs/"+x for x in get_files("./cosa/results_anaerobic/runs")]

json_paths = [
    "./cosa/results_anaerobic/runs/OPTMDF_STANDARDCONC_WILDTYPE.json",
    "./cosa/results_anaerobic/runs/OPTSUBMDF_STANDARDCONC_WILDTYPE.json",
]

active_reactions = []
for json_path in json_paths:
    json_path = json_path.replace(".zip", "")
    print(json_path)
    results = json_zip_load(json_path)

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
            if reaction_id.startswith("var_z"):
                continue
            if reaction_id.startswith("dG0"):
                continue
            flux = results[growth_rate]["values"][reaction_id]
            if flux > 1e-6:
                if reaction_id not in active_reactions:
                    active_reactions.append(
                        reaction_id
                            .replace("_ORIGINAL_NADP_TCOSA", "")
                            .replace("_ORIGINAL_NAD_TCOSA", "")
                            .replace("_VARIANT_NAD_TCOSA", "")
                            .replace("_VARIANT_NADP_TCOSA", "")
                            .replace("_FWD", "")
                            .replace("_REV", "")
                    )
print(active_reactions)
active_reactions = list(set(active_reactions))
print(f"Num  active reactions:", len(active_reactions))
