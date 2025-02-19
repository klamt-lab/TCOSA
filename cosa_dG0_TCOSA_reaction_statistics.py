"""Creates TCOSA-reaction related dG0 statistics, i.e., dG0 staistics for all NAD(P)(H)-containing reactions.

This data is especially used for the dG0 sampling in TCOSA's publication, and mentioned elsewhere there.
"""

# IMPORTS #
# External
import cobra
from statistics import mean, median
# Internal
from helper import json_load


# SCRIPT RUN #
model = cobra.io.read_sbml_model("resources/iML1515_irreversible_cleaned.xml")
dG0s = json_load("resources/dG0_iML1515_irreversible_cleaned.json")

nad_to_nadh = []
nadh_to_nad = []
nadp_to_nadph = []
nadph_to_nadp = []

num_with_nad_or_nadp = 0
tcosa_with_dG0 = []
for key in list(dG0s.keys()):
    if key.endswith("_REV"):
        continue
    if key in ["NADK", "NADTRHD", "NADPPPS", "THD2pp"]:
        continue

    dG0 = dG0s[key]["dG0"]
    if dG0s[key]["num_compartments"] > 1:
        continue

    reaction: cobra.Reaction = model.reactions.get_by_id(key)
    educt_met_ids = [x.id for x in reaction.metabolites.keys() if reaction.metabolites[x] < 0.0]
    product_met_ids = [x.id for x in reaction.metabolites.keys() if reaction.metabolites[x] > 0.0]

    nad_is_educt = "nad_c" in educt_met_ids
    nadh_is_educt = "nadh_c" in educt_met_ids
    nadp_is_educt = "nadp_c" in educt_met_ids
    nadph_is_educt = "nadph_c" in educt_met_ids


    nad_is_product = "nad_c" in product_met_ids
    nadh_is_product = "nadh_c" in product_met_ids
    nadp_is_product = "nadp_c" in product_met_ids
    nadph_is_product = "nadph_c" in product_met_ids

    is_nad_nadp = False
    if (nad_is_educt) and (nadh_is_product):
        nad_to_nadh.append(dG0)
        nadh_to_nad.append(-dG0)
        is_nad_nadp = True
    if (nadh_is_educt) and (nad_is_product):
        nadh_to_nad.append(dG0)
        nad_to_nadh.append(-dG0)
        is_nad_nadp = True
    if (nadp_is_educt) and (nadph_is_product):
        nadp_to_nadph.append(dG0)
        nadph_to_nadp.append(-dG0)
        is_nad_nadp = True
    if (nadph_is_educt) and (nadp_is_product):
        nadph_to_nadp.append(dG0)
        nadp_to_nadph.append(-dG0)
        is_nad_nadp = True
    if is_nad_nadp:
        tcosa_with_dG0.append(key.replace("_FWD", ""))
        print(key, dG0)
        num_with_nad_or_nadp += 1

print("DIRECTION MEAN[kJ/mol] MEDIAN[kJ/mol] N")
print("nad_to_nadh",   round(mean(nad_to_nadh), 3),   round(median(nad_to_nadh), 3), len(nad_to_nadh))
print("nadh_to_nad",   round(mean(nadh_to_nad), 3),   round(median(nadh_to_nad), 3), len(nadh_to_nad))
print("nadp_to_nadph", round(mean(nadp_to_nadph), 3), round(median(nadp_to_nadph), 3), len(nadp_to_nadph))
print("nadph_to_nadp", round(mean(nadph_to_nadp), 3), round(median(nadph_to_nadp), 3), len(nadph_to_nadp))
print("num_with_nad_or_nadp", num_with_nad_or_nadp)
print(len(set(tcosa_with_dG0)))
print(tcosa_with_dG0)
