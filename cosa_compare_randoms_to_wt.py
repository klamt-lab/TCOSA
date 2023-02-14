import cobra
from helper import json_load
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids

cobra_model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")

randoms47 = json_load("cosa/results_aerobic/randoms_rand_lists.json")["47"]

num_same = 0
num_different = 0
for tcosa_reaction_id in get_all_tcosa_reaction_ids(cobra_model):
    if "_VARIANT" in tcosa_reaction_id:
        continue

    if not "_ORIGINAL" in tcosa_reaction_id:
        continue

    if "_ORIGINAL_NADP" in tcosa_reaction_id:
        wt_state = "NADP"
    elif "_ORIGINAL_NAD" in tcosa_reaction_id:
        wt_state = "NAD"
    else:
        print("~", tcosa_reaction_id)
        continue

    base_id = tcosa_reaction_id.replace(f"_ORIGINAL_{wt_state}_TCOSA", "")
    randoms_statenum = randoms47[base_id]

    if randoms_statenum == 0.0:
        randoms_state = "NADP"
    elif randoms_statenum == 1.0:
        randoms_state = "NAD"
    else:
        input("ERROR")

    if randoms_state == wt_state:
        num_same += 1
    else:
        num_different += 1

print("Same:", num_same)
print("Different:", num_different)
