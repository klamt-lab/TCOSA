import json
from tokenize import cookie_re
from helper import json_load
import cobra
from math import exp
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids

dG0_values = json_load("cosa/dG0_values.json")
cobra_model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")

# 1.X.1.X
oxidoreductase_counter = 0
other_counter = 0
for reaction_id in get_all_tcosa_reaction_ids(cobra_model):
    reaction: cobra.Reaction = cobra_model.reactions.get_by_id(reaction_id)

    if "_ORIGINAL_" not in reaction.id:
        if reaction.id not in ["NADK", "NADTRHD", "NADPPPS", "THD2pp"]:
            continue

    if "ec-code" in reaction.annotation.keys():
        ec_codes = reaction.annotation["ec-code"]
        if type(ec_codes) is str:
            ec_codes = [ec_codes]
        is_oxidoreductase = False
        for ec_code in ec_codes:
            ec_numbers = ec_code.split(".")
            if ((ec_numbers[0] == "1") or (ec_numbers[0] == "-")) and ((ec_numbers[2] == "1") or (ec_numbers[2] == "-")):
                is_oxidoreductase = True

        if is_oxidoreductase:
            oxidoreductase_counter += 1
        else:
            print(reaction.id, ec_codes, reaction.reaction)
            other_counter += 1

print("# oxidoreductases: ", oxidoreductase_counter)
print("# other: ", other_counter)
