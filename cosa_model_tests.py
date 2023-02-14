import cobra
import copy
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from helper import json_zip_load


all_base_ids, original_cobra_model, concentration_values_free, concentration_values_paper,\
standardconc_dG0_values, paperconc_dG0_values,\
num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=False, expanded=False)

print("====")

# TEST 1: Check that NADX_ONLY really has only NADX
print("Performing TEST 1: Check that SINGLE_COFACTOR really has only NADX")
only_testmodel = copy.deepcopy(original_cobra_model)
only_testmodel = cosa_get_model_with_nadx_scenario(
    nadx_scenario="SINGLE_COFACTOR",
    cobra_model=only_testmodel,
    randoms_random_base_lists=[],
    randomfixed_random_base_lists=[],
)
only_test_correct = True
for reaction in only_testmodel.reactions:
    reaction: cobra.Reaction = reaction
    if reaction.id.endswith("_NADY"):
        abs_lower_bound = abs(reaction.lower_bound)
        abs_upper_bound = abs(reaction.upper_bound)

        if (abs_lower_bound > 0.0) or (abs_upper_bound > 0.0):
            only_test_correct = False
            print(f"ERROR in {reaction.id}")
if only_test_correct:
    print("TEST 1 SUCCESSFUL")
else:
    print("TEST 1 FAILED! SEE ERROR MESSAGES")
print("====")


# TEST 2: Check that WILDTYPE really has only in vivo reactions
print("Performing TEST 2: Check that WILDTYPE really has only in vivo reactions")
vivo_testmodel = copy.deepcopy(original_cobra_model)
vivo_testmodel = cosa_get_model_with_nadx_scenario(
    nadx_scenario="WILDTYPE",
    cobra_model=vivo_testmodel,
    randoms_random_base_lists=[],
    randomfixed_random_base_lists=[],
)
vivo_test_correct = True
for reaction in vivo_testmodel.reactions:
    reaction: cobra.Reaction = reaction
    if reaction.id.endswith("_NADX") and (reaction.id in nadp_reactions):
        abs_lower_bound = abs(reaction.lower_bound)
        abs_upper_bound = abs(reaction.upper_bound)

        if (abs_lower_bound > 0.0) or (abs_upper_bound > 0.0):
            vivo_test_correct = False
            print(f"ERROR in {reaction.id}")
    if reaction.id.endswith("_NADY") and (reaction.id in nad_reactions):
        abs_lower_bound = abs(reaction.lower_bound)
        abs_upper_bound = abs(reaction.upper_bound)

        if (abs_lower_bound > 0.0) or (abs_upper_bound > 0.0):
            vivo_test_correct = False
            print(f"ERROR in {reaction.id}")
if vivo_test_correct:
    print("TEST 2 SUCCESSFUL")
else:
    print("TEST 2 FAILED! SEE ERROR MESSAGES")
print("====")


# TEST 3: Test
result = json_zip_load("cosa/results_aerobic/runs/OPTMDF_STANDARDCONC_FLEXIBLE.json")["0,868"]
base_ids = list(result["values"].keys())
base_ids = set([x.replace("_NADX", "").replace("_NADY", "") for x in base_ids if (x.endswith("_NADX") or x.endswith("_NADY")) and (not "_var_" in x) and (not "dG0_" in x)])
print(base_ids)

for base_id in base_ids:
    z_id_1 = "z_var_"+base_id+"_NADX"
    z_id_2 = "z_var_"+base_id+"_NADY"
    z_sum = result["values"][z_id_1] + result["values"][z_id_2]
    print(z_sum)
