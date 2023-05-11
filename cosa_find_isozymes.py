import cobra
from helper import json_write

model = cobra.io.read_sbml_model("resources/iML1515.xml")

"""
reactions_per_ec = {}
for reaction in model.reactions:
    reaction: cobra.Reaction = reaction

    ids = [x.id for x in reaction.metabolites]

    if not (("nad_c" in ids) or ("nadh_c" in ids) or ("nadp_c" in ids) or ("nadph_c" in ids)):
        continue

    if "ec-code" not in reaction.annotation.keys():
        continue

    ec_codes = reaction.annotation["ec-code"]

    if type(ec_codes) is str:
        ec_codes = [ec_codes]

    for ec_code in ec_codes:
        if ec_code not in list(reactions_per_ec.keys()):
            reactions_per_ec[ec_code] = []

        reactions_per_ec[ec_code].append(f"{reaction.id}: {reaction.reaction}")

keys = list(reactions_per_ec.keys())
for key in keys:
    if len(reactions_per_ec[key]) == 1:
        del(reactions_per_ec[key])
        continue
"""
def multi_replace(string: str) -> str:
    string = string.replace("nad_c", "nadx").replace("nadp_c", "nadx").replace("nadh_c", "nadxh").replace("nadph_c", "nadxh")
    return string

def make_tuple_from_reaction(reaction: cobra.Reaction):
    ids = [x.id for x in reaction.metabolites]
    if not (("nad_c" in ids) or ("nadh_c" in ids) or ("nadp_c" in ids) or ("nadph_c" in ids)):
        return None

    temp = [
        str(reaction.metabolites[x])+" "+multi_replace(x.id)
        for x in reaction.metabolites
    ]
    temp.sort()
    return (reaction.id, reaction.reaction, tuple(temp))

reactions_per_ec = {}
reaction_tuples = [make_tuple_from_reaction(x) for x in model.reactions]
reaction_tuples = [x for x in reaction_tuples if x is not None]

counter = 0
found_tuples = []
reactions_per_ec = {}
pos_1 = 0
for reaction_tuple in reaction_tuples:
    pos_2 = 0
    for reaction_tuple_2 in reaction_tuples:
        if pos_1 != pos_2:
            if reaction_tuple[2] == reaction_tuple_2[2]:
                reactions_per_ec[counter] = [
                    reaction_tuple[0] + " " + reaction_tuple[1],
                    reaction_tuple_2[0] + " " + reaction_tuple_2[1],
                ]
                counter += 1
        pos_2 += 1
    pos_1 += 1

json_write(
    "./cosa/isozymes.json",
    reactions_per_ec
)

"""
"SSALx: h2o_c + nad_c + sucsal_c --> 2.0 h_c + nadh_c + succ_c"
"SSALy: h2o_c + nadp_c + sucsal_c --> 2.0 h_c + nadph_c + succ_c"


"""