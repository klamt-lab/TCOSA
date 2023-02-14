import cobra
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from helper import json_load, pickle_load


dG0_values = json_load("cosa/dG0_values.json")
cobra_model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")

total_nad_nadp_reactions = []
nad_reactions_with_fwd_and_rev = []
nadp_reactions_with_fwd_and_rev = []
nad_reversible_reactions = []
nad_reactions = []
nadp_reactions = []
nadp_reversible_reactions = []
nad_and_nadp_reactions = []
nad_total_reactions = []
nadp_total_reactions = []

total_num_reactions = len(cobra_model.reactions)
total_num_metabolites = len(cobra_model.metabolites)
for reaction in cobra_model.reactions:
    metabolite_ids = [x.id for x in reaction.metabolites]
    has_nad = "nad_tcosa_c" in metabolite_ids
    has_nadh = "nadh_tcosa_c" in metabolite_ids
    has_nadp = "nadp_tcosa_c" in metabolite_ids
    has_nadph = "nadph_tcosa_c" in metabolite_ids
    has_any_nad = has_nad or has_nadh
    has_any_nadp = has_nadp or has_nadph

    # Original model
    if has_any_nad and has_any_nadp:
        nad_and_nadp_reactions.append(reaction.id)

    if reaction.id not in get_all_tcosa_reaction_ids(cobra_model):
        continue

    # TCOSA model
    total_nad_nadp_reactions.append(reaction.id)

    if not (reaction.id.endswith("_TCOSA")):
        continue

    if has_any_nad:
        nad_total_reactions.append(reaction.id)
    elif has_any_nadp:
        nadp_total_reactions.append(reaction.id)

    if ("_VARIANT_" in reaction.id):
        continue

    if ("_REV_ORIGINAL_" in reaction.id):
        continue

    if has_any_nad:
        nad_reactions.append(reaction.id)
        if "_FWD_" in reaction.id:
            nad_reversible_reactions.append(reaction.id)
    elif has_any_nadp:
        nadp_reactions.append(reaction.id)
        if "_FWD_" in reaction.id:
            nadp_reversible_reactions.append(reaction.id)

print("==IN TCOSA MODEL==")
print("Total # metabolites:", total_num_metabolites)
print("Total # reactions:", total_num_reactions)
print("Total # TCOSA reactions:", len(total_nad_nadp_reactions))
print("Total # NAD(H) *AND* NADP(H) REACTIONS:", len(nad_and_nadp_reactions))
print("Total # NAD(H) reactions:", len(nad_total_reactions))
print("Total # NADP(H) reactions:", len(nadp_total_reactions))
print("")
print("==ORIGINAL MODEL==")
print("Total # NAD(H) *AND* NADP(H) REACTIONS:", len(nad_and_nadp_reactions))
print("Total # NAD(H) REACTIONS:", len(nad_reactions))
print(" ...OF WHICH REVERSIBLE:", len(nad_reversible_reactions))
print("Total # NADP(H) REACTIONS:", len(nadp_reactions))
print(" ...OF WHICH REVERSIBLE:", len(nadp_reversible_reactions))

print(nad_total_reactions)
print(len(nad_total_reactions))
