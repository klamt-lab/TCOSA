"""Here, starting with the previouly cleaned iML1515 model, the TCOSA-adapted iML1515_TCOSA is generated.

In the course of this process, the model gets all TCOSA reaction duplications.
"""

import cobra
import copy
from helper import ensure_folder_existence, json_write, pickle_write

ensure_folder_existence("./cosa")

# Load model
original_cobra_model: cobra.Model = cobra.io.read_sbml_model("resources/iML1515_irreversible_cleaned.xml")
original_cobra_model.solver = "cplex"

# Deactivate POR5 as in the original OptMDFpathway publication (it has no biological equivalent)
original_cobra_model.reactions.get_by_id("POR5_FWD").lower_bound = 0
original_cobra_model.reactions.get_by_id("POR5_FWD").upper_bound = 0
original_cobra_model.reactions.get_by_id("POR5_REV").lower_bound = 0
original_cobra_model.reactions.get_by_id("POR5_REV").upper_bound = 0

changed_cobra_model = copy.deepcopy(original_cobra_model)

# Original NAD(P)(H) metabolite IDs
nad_ids = ["nad_c", "nadh_c"]
nadp_ids = ["nadp_c", "nadph_c"]

# NAD(P)(H) substitute metabolites
nadx = cobra.Metabolite(id="nad_tcosa_c", name="NAD for TCOSA", compartment="c")
nady = cobra.Metabolite(id="nadp_tcosa_c", name="NADP for TCOSA", compartment="c")
nadz = cobra.Metabolite(id="nadz_tcosa_c", name="3rd NAD(P) for TCOSA", compartment="c")
nadxh = cobra.Metabolite(id="nadh_tcosa_c", name="NADH for TCOSA", compartment="c")
nadyh = cobra.Metabolite(id="nadph_tcosa_c", name="NADPH for TCOSA", compartment="c")
nadzh = cobra.Metabolite(id="nadzh_tcosa_c", name="3rd NAD(P)H for TCOSA", compartment="c")

# Biomass reactions
biomass1 = copy.deepcopy(changed_cobra_model.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M"))
biomass1.add_metabolites({
    changed_cobra_model.metabolites.get_by_id("nad_c"): 0.001831,
    nadx: -0.001831,
    changed_cobra_model.metabolites.get_by_id("nadp_c"): 0.000447,
    nady: -0.000447,
})
biomass2 = copy.deepcopy(changed_cobra_model.reactions.get_by_id("BIOMASS_Ec_iML1515_WT_75p37M"))
biomass2.add_metabolites({
    changed_cobra_model.metabolites.get_by_id("nad_c"): 0.001787,
    nadx: -0.001787,
    changed_cobra_model.metabolites.get_by_id("nadp_c"): 0.000112,
    nady: -0.000112,
    changed_cobra_model.metabolites.get_by_id("nadh_c"): 4.5e-05,
    nadxh: -4.5e-05,
    changed_cobra_model.metabolites.get_by_id("nadph_c"): 0.000335,
    nadyh: -0.000335,
})
changed_cobra_model.remove_reactions([
    changed_cobra_model.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M"),
    changed_cobra_model.reactions.get_by_id("BIOMASS_Ec_iML1515_WT_75p37M")
])
changed_cobra_model.add_reactions([biomass1, biomass2])

# Hard-coded reactions where both NAD and NADP occur
## NADK
nadk = copy.deepcopy(changed_cobra_model.reactions.get_by_id("NADK"))
nadk.add_metabolites({
    changed_cobra_model.metabolites.get_by_id("nad_c"): 1,
    nadx: -1,
    changed_cobra_model.metabolites.get_by_id("nadp_c"): -1,
    nady: 1,
})

## NADTRHD
natrhd = copy.deepcopy(changed_cobra_model.reactions.get_by_id("NADTRHD"))
natrhd.add_metabolites({
    changed_cobra_model.metabolites.get_by_id("nad_c"): 1,
    nadx: -1,
    changed_cobra_model.metabolites.get_by_id("nadph_c"): 1,
    nadyh: -1,
    changed_cobra_model.metabolites.get_by_id("nadh_c"): -1,
    nadxh: 1,
    changed_cobra_model.metabolites.get_by_id("nadp_c"): -1,
    nady: 1,
})
## THD2pp
thd2pp = copy.deepcopy(changed_cobra_model.reactions.get_by_id("THD2pp"))
thd2pp.add_metabolites({
    changed_cobra_model.metabolites.get_by_id("nadh_c"): 1,
    nadxh: -1,
    changed_cobra_model.metabolites.get_by_id("nadp_c"): 1,
    nady: -1,
    changed_cobra_model.metabolites.get_by_id("nad_c"): -1,
    nadx: 1,
    changed_cobra_model.metabolites.get_by_id("nadph_c"): -1,
    nadyh: 1,
})
## NADPPPS
nadppps = copy.deepcopy(changed_cobra_model.reactions.get_by_id("NADPPPS"))
nadppps.add_metabolites({
    changed_cobra_model.metabolites.get_by_id("nadp_c"): 1,
    nady: -1,
    changed_cobra_model.metabolites.get_by_id("nad_c"): -1,
    nadx: 1,
})

single_cofactor_pseudoreaction = cobra.Reaction(
    id="Single_Cofactor_Pseudoreaction",
    lower_bound=0,
    upper_bound=0,
)
single_cofactor_pseudoreaction.add_metabolites({
    nadx: -1,
    nady: 1,
})

changed_cobra_model.remove_reactions([
    changed_cobra_model.reactions.get_by_id("NADK"),
    changed_cobra_model.reactions.get_by_id("NADTRHD"),
    changed_cobra_model.reactions.get_by_id("THD2pp"),
    changed_cobra_model.reactions.get_by_id("NADPPPS"),
])
changed_cobra_model.add_reactions([nadk, natrhd, thd2pp, nadppps, single_cofactor_pseudoreaction])

changed_cobra_model_expanded = copy.deepcopy(changed_cobra_model)

# Create TCOSA reaction duplications
reaction_ids = [x.id for x in original_cobra_model.reactions]
reactions_to_remove = []
dG0s_to_remove = []
for reaction_id in reaction_ids:
    if "BIOMASS_Ec" in reaction_id:
        continue
    if reaction_id in ("NADK", "NADTRHD", "THD2pp", "NADPPPS"):
        continue

    reaction: cobra.Reaction = original_cobra_model.reactions.get_by_id(reaction_id)
    metabolite_ids = [x.id for x in reaction.metabolites.keys()]

    has_nad = False
    has_nadp = False
    for metabolite_id in metabolite_ids:
        if metabolite_id in nad_ids:
            has_nad = True
        if metabolite_id in nadp_ids:
            has_nadp = True


    if has_nad and has_nadp:
        dG0s_to_remove.append(reaction_id)
        print(reaction_id)

    if has_nad or has_nadp:
        current_copy = 1
        for current_nadxs in ([nadx, nadxh], [nady, nadyh], [nadz, nadzh]):
            reaction_copy = copy.deepcopy(reaction)
            if current_copy == 1:
                if has_nad:
                    reaction_copy.id += "_ORIGINAL_NAD_TCOSA"
                else:
                    reaction_copy.id += "_VARIANT_NAD_TCOSA"
            elif current_copy == 2:
                if has_nadp:
                    reaction_copy.id += "_ORIGINAL_NADP_TCOSA"
                else:
                    reaction_copy.id += "_VARIANT_NADP_TCOSA"
            elif current_copy == 3:
                reaction_copy.id += "_NADZ_TCOSA"
            metabolite_dict = reaction.metabolites
            for key in metabolite_dict.keys():
                if (key.id in nad_ids) or (key.id in nadp_ids):
                    stoichiometry = reaction.metabolites[key]
                    if "h" in key.id:
                        current_nadx = current_nadxs[1]
                    else:
                        current_nadx = current_nadxs[0]
                    reaction_copy.add_metabolites({
                        key: -stoichiometry,
                        current_nadx: stoichiometry,
                    })
            if current_copy <= 2:
                changed_cobra_model.add_reactions([reaction_copy])
            changed_cobra_model_expanded.add_reactions([reaction_copy])
            current_copy += 1

        reactions_to_remove.append(reaction_id)

print(changed_cobra_model.reactions.get_by_id("Single_Cofactor_Pseudoreaction"))
changed_cobra_model.remove_reactions(reactions=reactions_to_remove)
changed_cobra_model_expanded.remove_reactions(reactions=reactions_to_remove)
cobra.io.write_sbml_model(changed_cobra_model, "./cosa/iML1515_TCOSA.xml")
pickle_write("./cosa/iML1515_TCOSA.pickle", changed_cobra_model)
cobra.io.write_sbml_model(changed_cobra_model_expanded, "./cosa/iML1515_3TCOSA_expanded.xml")
pickle_write("./cosa/iML1515_3TCOSA_expanded.pickle", changed_cobra_model_expanded)

print("FLEXIBLE")
print(original_cobra_model.objective)
print(original_cobra_model.slim_optimize())
changed_cobra_model.objective = "BIOMASS_Ec_iML1515_core_75p37M"
print(changed_cobra_model.objective)
print(changed_cobra_model.slim_optimize())
