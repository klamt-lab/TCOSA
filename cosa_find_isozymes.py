"""Identifies 'isozyme-like' reactions, i.e., reactions which only differ in the usage of NAD or NADP, all in iML1515.

The found reactions are stored as 'isozymes.json'.
"""

# IMPORTS #
# External
import cobra
# Internal
from helper import json_write


# PUBLIC FUNCTIONS #
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
    return (reaction.id, reaction.reaction, tuple(temp), reaction.gene_reaction_rule)


# LOGIC #
model = cobra.io.read_sbml_model("resources/iML1515.xml")

reactions_per_ec = {}
reaction_tuples = [make_tuple_from_reaction(x) for x in model.reactions]
reaction_tuples = [x for x in reaction_tuples if x is not None]

found_reaction_rules = []
num_same_enzyme_rule = 0
num_different_enzyme_rule = 0
counter = 0
found_tuples = []
reactions_per_ec = {}
found_reactions_up_to_now = []
pos_1 = 0
for reaction_tuple in reaction_tuples:
    pos_2 = 0
    for reaction_tuple_2 in reaction_tuples:
        if pos_1 != pos_2:
            if reaction_tuple[2] == reaction_tuple_2[2]:
                if reaction_tuple[0] in found_reactions_up_to_now:
                    continue
                reactions_per_ec[counter+1] = [
                    reaction_tuple[0] + " " + reaction_tuple[1] + " : " + "".join(reaction_tuple[3]),
                    reaction_tuple_2[0] + " " + reaction_tuple_2[1] + " : " + "".join(reaction_tuple_2[3]),
                ]
                found_reactions_up_to_now.append(reaction_tuple[0])
                found_reactions_up_to_now.append(reaction_tuple_2[0])
                if reaction_tuple[3] == reaction_tuple_2[3]:
                    num_same_enzyme_rule += 1
                else:
                    num_different_enzyme_rule += 1
                found_reaction_rules.append(reaction_tuple[3])
                found_reaction_rules.append(reaction_tuple_2[3])
                counter += 1
        pos_2 += 1
    pos_1 += 1

found_reaction_rules = list(set(found_reaction_rules))

reactions_per_ec["num_same_enzyme_rule"] = num_same_enzyme_rule
reactions_per_ec["num_different_enzyme_rule"] = num_different_enzyme_rule
reactions_per_ec["num_found_different_reaction_rules"] = len(found_reaction_rules)

json_write(
    "./cosa/isozymes.json",
    reactions_per_ec
)
