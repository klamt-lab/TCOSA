import cobra
from helper import json_zip_load, json_load, json_write, get_files
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids

dG0_values = json_load("cosa/dG0_values.json")
cobra_model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")

def print_reacs(reacs):
    for reac in reacs:
        cobra_reac = cobra_model.reactions.get_by_id(reac)
        try:
            dG0 = round(dG0_values[reac]["dG0"], 2)
        except KeyError:
            dG0 = float("NaN")
        print(f"* {reac};  ΔG0={dG0} kJ/mol; {cobra_reac.reaction}")

tcosa_reactions = get_all_tcosa_reaction_ids(cobra_model)


def _valid_id(var_id: str) -> bool:
    if var_id.startswith("z_var_") and ("MDF_A" not in var_id):
        reac_id = var_id.replace("z_var_", "")
        if reac_id not in tcosa_reactions:
            return False
        if reac_id in ["NADK", "NADTRHD", "NADPPPS", "THD2pp"]:
            return False
        return True

reaction_wise_statistics = {}
all_num_originals = []
all_num_variants = []
for folder_number in ["3"]: # ["", "2"]:
    folder = f"cosa/results_aerobic/mult_flux_samplings{folder_number}/"
    for file in get_files(folder):
        result = json_load(folder+file)
        num_originals = 0
        num_variants = 0
        for var_id in result["values"].keys():
            if var_id.startswith("z_var_") and ("MDF_A" not in var_id):
                if not _valid_id(var_id=var_id):
                    continue

                var_value = result["values"][var_id]
                if var_value < 1e-3:
                    continue
                # print("*", var_id, var_value)
                clear_var_id = var_id.replace("_ORIGINAL_NAD_", "")
                clear_var_id = clear_var_id.replace("_ORIGINAL_NADP_", "")
                clear_var_id = clear_var_id.replace("_VARIANT_NADP_", "")
                clear_var_id = clear_var_id.replace("_VARIANT_NAD_", "")
                if clear_var_id not in reaction_wise_statistics.keys():
                    reaction_wise_statistics[clear_var_id] = {
                        "original": 0,
                        "variant": 0,
                    }

                if "_ORIGINAL_NA" in var_id:
                    num_originals += 1
                    reaction_wise_statistics[clear_var_id]["original"] += 1
                else:
                    num_variants += 1
                    reaction_wise_statistics[clear_var_id]["variant"] += 1

        print(f"No.: {file.split('_SAMPLING_')[-1]}, num_originals: {num_originals}, num_variants: {num_variants}")
        all_num_originals.append(num_originals)
        all_num_variants.append(num_variants)

    mean_num_originals = sum(all_num_originals) / len(all_num_originals)
    mean_num_variants = sum(all_num_variants) / len(all_num_variants)
    print(reaction_wise_statistics)
    print(f"MEANS: mean_num_originals: {round(mean_num_originals, 2)}, mean_num_variants: {round(mean_num_variants, 2)}")
    text = ""
    for reac_id in reaction_wise_statistics.keys():
        text += f"{reac_id.replace('TCOSA', '_TCOSA')}; num_original: {reaction_wise_statistics[reac_id]['original']}; num_variant:  {reaction_wise_statistics[reac_id]['variant']}\n"

    with open(f"conditions_{folder_number}.txt", "w") as f:
        f.write(text)


    ## Correlation statistics
    correlation_stats = {}
    counter = 0
    for file in get_files(folder):
        result = json_load(folder+file)
        print(f"File {counter}")
        counter += 1

        for var_id in result["values"].keys():
            if not _valid_id(var_id):
                continue
            var_value = result["values"][var_id]
            if var_value < 1e-3:
                continue

            if var_id not in correlation_stats.keys():
                correlation_stats[var_id] = {
                    "with_original": [],
                    "with_variant": [],
                }
            correlation_stats[var_id]["with_original"].append(0)
            correlation_stats[var_id]["with_variant"].append(0)

            for other_var_id in result["values"].keys():
                if not _valid_id(other_var_id):
                    continue
                other_var_value = result["values"][other_var_id]
                if other_var_value < 1e-3:
                    continue
                if "_ORIGINAL_" in other_var_id:
                    correlation_stats[var_id]["with_original"][-1] += 1
                else:
                    correlation_stats[var_id]["with_variant"][-1] += 1

    correlation_stats_means = {}
    for var_id in correlation_stats.keys():
        # print(var_id)
        # print(correlation_stats[var_id])

        n = len(correlation_stats[var_id]["with_original"]) # + len(correlation_stats[var_id]["with_variant"])
        correlation_stats_means[var_id] = {
            "n": n,
            "mean_with_original": sum(correlation_stats[var_id]["with_original"]) / len(correlation_stats[var_id]["with_original"]),
            "min_with_original": min(correlation_stats[var_id]["with_original"]),
            "max_with_original": max(correlation_stats[var_id]["with_original"]),
            "mean_with_variant": sum(correlation_stats[var_id]["with_variant"]) / len(correlation_stats[var_id]["with_variant"]),
            "min_with_variant": min(correlation_stats[var_id]["with_variant"]),
            "max_with_variant": max(correlation_stats[var_id]["with_variant"]),
        }
    json_write(f"./correlation_stats_means{folder_number}.json", correlation_stats_means)

    num_originals = 0
    num_variants = 0
    major_means = {
        "mean_variant_with_variant": 0.0,
        "mean_variant_with_original": 0.0,
        "mean_original_with_original": 0.0,
        "mean_original_with_variant": 0.0,
    }
    for key in correlation_stats_means.keys():
        if correlation_stats_means[key]["n"] < 10:
            continue

        sum_ = correlation_stats_means[key]["mean_with_variant"] + correlation_stats_means[key]["mean_with_original"]
        if "_VARIANT_" in key:
            num_variants += 1
            major_means["mean_variant_with_variant"] += correlation_stats_means[key]["mean_with_variant"] / sum_
            major_means["mean_variant_with_original"] += correlation_stats_means[key]["mean_with_original"] / sum_
        else:
            num_originals += 1
            major_means["mean_original_with_variant"] += correlation_stats_means[key]["mean_with_variant"] / sum_
            major_means["mean_original_with_original"] += correlation_stats_means[key]["mean_with_original"] / sum_

    major_means["mean_variant_with_variant"] /= num_variants
    major_means["mean_variant_with_original"] /= num_variants
    major_means["mean_original_with_variant"] /= num_originals
    major_means["mean_original_with_original"] /= num_originals

    json_write(f"./c_major_means_cleaned{folder_number}.json", major_means)
