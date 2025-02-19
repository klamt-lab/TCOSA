"""Creates tables with single cofactor swap statistics as used in TCOSA's publication.

These tables are printed to the console.
"""

# IMPORTS #
# External
import cobra
from statistics import mean
# Internal
from helper import json_load

# LOGIC
model = cobra.io.read_sbml_model("resources/iML1515.xml")

make_base = lambda y: [
    x.replace("_ORIGINAL", "").replace("_VARIANT", "").replace("_NADP_TCOSA", "").replace("_NAD_TCOSA", "").replace("_FWD", "").replace("_REV", "") for x in y
]

all_ids = []
for mode in ("GREATER_THAN", "LOWER_THAN"):
    print(f"===={mode}====")
    for concentrations in ("STANDARDCONC", "VIVOCONC",):
        max_optmdf_change = 0
        max_optsubmdf_change = 0
        min_optmdf_change = 0
        min_optsubmdf_change = 0
        print(f"=={concentrations}==")

        if concentrations == "STANDARDCONC":
            aerobicities = ("aerobic", "anaerobic")
        else:
            aerobicities = ("aerobic",)

        for aerobicity in aerobicities:
            print(f"={aerobicity}=")
            data_path = f"./cosa/results_{aerobicity}/swap_results_{concentrations}.json"
            swap_data = json_load(data_path)

            id_optmdf_at_all = []
            id_optsubmdf_at_all =[]
            id_optmdf_at_one = []
            id_optsubmdf_at_one = []

            positive_optmdf_changes = []
            negative_optmdf_changes = []
            positive_optsubmdf_changes = []
            negative_optsubmdf_changes = []
            for reac_id in swap_data.keys():
                if len(swap_data[reac_id].keys()) <= 1:
                    continue

                optmdf_changes = []
                optsubmdf_changes = []
                error = False
                for growth_rate in swap_data[reac_id].keys():
                    if len(swap_data[reac_id][growth_rate].keys()) != 2:
                        error = True
                        break
                    optmdf_changes.append(swap_data[reac_id][growth_rate]["OptMDF"])
                    optsubmdf_changes.append(swap_data[reac_id][growth_rate]["OptSubMDF"])
                if error:
                    continue

                optmdf_at_all = False
                optsubmdf_at_all = False
                optmdf_at_one = False
                optsubmdf_at_one = False
                if mode == "GREATER_THAN":
                    if (max(optmdf_changes) > 0.0) and (min(optmdf_changes) > 0.0) and (0.0 not in optmdf_changes):
                        optmdf_at_all = True
                    elif (max(optmdf_changes) > 0.0):
                        optmdf_at_one = True

                    if (max(optsubmdf_changes) > 0.0) and (min(optsubmdf_changes) > 0.0) and  (0.0 not in optsubmdf_changes):
                        optsubmdf_at_all = True
                    elif (max(optsubmdf_changes) > 0.0):
                        optsubmdf_at_one = True
                elif mode == "LOWER_THAN":
                    if (max(optmdf_changes) < 0.0) and (min(optmdf_changes) < 0.0) and (0.0 not in optmdf_changes):
                        optmdf_at_all = True
                    elif (min(optmdf_changes) < 0.0):
                        optmdf_at_one = True

                    if (max(optsubmdf_changes) < 0.0) and (min(optsubmdf_changes) < 0.0) and (0.0 not in optsubmdf_changes):
                        optsubmdf_at_all = True
                    elif (min(optsubmdf_changes) < 0.0):
                        optsubmdf_at_one = True
                elif mode == "BOTH":
                    if (0.0 not in optmdf_changes):
                        optmdf_at_all = True
                    elif (max(optmdf_changes) != 0.0) or (min(optmdf_changes) != 0.0):
                        optmdf_at_one = True

                    if (0.0 not in optsubmdf_changes):
                        optsubmdf_at_all = True
                    elif (max(optsubmdf_changes) != 0.0) or (min(optsubmdf_changes) != 0.0):
                        optsubmdf_at_one = True
                elif mode == "STRICT_BOTH":
                    if (max(optmdf_changes) < 0.0) and (min(optmdf_changes) > 0.0):
                        optmdf_at_all = True

                    if (max(optsubmdf_changes) < 0.0) and (min(optsubmdf_changes) > 0.0):
                        optsubmdf_at_all = True

                ####
                max_optmdf_change = max(max_optmdf_change, max(optmdf_changes))
                max_optsubmdf_change = max(max_optsubmdf_change, max(optsubmdf_changes))
                min_optmdf_change = min(min_optmdf_change, min(optmdf_changes))
                min_optsubmdf_change = min(min_optsubmdf_change, min(optsubmdf_changes))

                positive_optmdf_changes += [x for x in optmdf_changes if x > 0.0]
                negative_optmdf_changes += [x for x in optmdf_changes if x < 0.0]

                positive_optsubmdf_changes += [x for x in optsubmdf_changes if x > 0.0]
                negative_optsubmdf_changes += [x for x in optsubmdf_changes if x < 0.0]

                if optmdf_at_all:
                    optmdf_at_one = True
                    optmdf_at_all = False
                if optsubmdf_at_all:
                    optsubmdf_at_one = True
                    optsubmdf_at_all = False
                ####

                if optmdf_at_all:
                    id_optmdf_at_all.append(reac_id)
                if optsubmdf_at_all:
                    id_optsubmdf_at_all.append(reac_id)
                if optmdf_at_one:
                    id_optmdf_at_one.append(reac_id)
                if optsubmdf_at_one:
                    id_optsubmdf_at_one.append(reac_id)

            id_optmdf_at_all = list(set(make_base(id_optmdf_at_all)))
            id_optsubmdf_at_all = list(set(make_base(id_optsubmdf_at_all)))
            id_optmdf_at_one = list(set(make_base(id_optmdf_at_one)))
            id_optsubmdf_at_one = list(set(make_base(id_optsubmdf_at_one)))

            # print(f"OptMDF at all (n={len(id_optmdf_at_all)}): {id_optmdf_at_all}")
            # print(f"OptSubMDF at all (n={len(id_optsubmdf_at_all)}): {id_optsubmdf_at_all}")
            print(f"OptMDF at >= 1 (n={len(id_optmdf_at_one)}): {id_optmdf_at_one}")
            print(f"OptSubMDF at >= 1 (n={len(id_optsubmdf_at_one)}): {id_optsubmdf_at_one}")
            print("~")
            print("Mean positive OptMDF change:", mean(positive_optmdf_changes))
            print("Mean negative OptMDF change:", mean(negative_optmdf_changes))
            print("Mean positive OptSubMDF change:", mean(positive_optsubmdf_changes))
            print("Mean negative OptSubMDF change:", mean(negative_optsubmdf_changes))

            all_ids += id_optmdf_at_all + id_optsubmdf_at_all + id_optmdf_at_one + id_optsubmdf_at_one

        print("###########################")
        print("#########END###############")
        print("###########################")

        print(f"Min OptMDF change: {min_optmdf_change} kJ/mol")
        print(f"Max OptMDF change: {max_optmdf_change} kJ/mol")
        print(f"Min OptSubMDF change: {min_optsubmdf_change} kJ/mol")
        print(f"Max OptSubMDF change: {max_optsubmdf_change} kJ/mol")

        all_ids = list(set(all_ids))
        print(f"All IDs (n={len(all_ids)}): {all_ids}")

        # for id_ in all_ids:
        #     reaction = model.reactions.get_by_id(id_)
        #     print(reaction.id, reaction.reaction)
