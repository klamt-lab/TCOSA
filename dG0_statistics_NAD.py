import cobra
from statistics import mean, median
from helper import json_load


model = cobra.io.read_sbml_model("resources/iML1515_irreversible_cleaned.xml")
dG0s = json_load("resources/dG0_iML1515_irreversible_cleaned.json")

nad_to_nadh = []
nadh_to_nad = []
nadp_to_nadph = []
nadph_to_nadp = []

for key in list(dG0s.keys()):
    if key.endswith("_REV"):
        continue
    if key in ["NADK", "NADTRHD", "NADPPPS", "THD2pp"]:
        continue
    if dG0s[key]["num_compartments"] > 1:
        continue

    dG0 = dG0s[key]["dG0"]

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

    if (nad_is_educt) and (nadh_is_product):
        nad_to_nadh.append(dG0)
    if (nadh_is_educt) and (nad_is_product):
        nadh_to_nad.append(dG0)
    if (nadp_is_educt) and (nadph_is_product):
        nadp_to_nadph.append(dG0)
    if (nadph_is_educt) and (nadp_is_product):
        nadph_to_nadp.append(dG0)

print("DIRECTION MEAN[kJ/mol] MEDIAN[kJ/mol] N")
print("nad_to_nadh",   round(mean(nad_to_nadh), 3),   round(median(nad_to_nadh), 3), len(nad_to_nadh))
print("nadh_to_nad",   round(mean(nadh_to_nad), 3),   round(median(nadh_to_nad), 3), len(nadh_to_nad))
print("nadp_to_nadph", round(mean(nadp_to_nadph), 3), round(median(nadp_to_nadph), 3), len(nadp_to_nadph))
print("nadph_to_nadp", round(mean(nadph_to_nadp), 3), round(median(nadph_to_nadp), 3), len(nadph_to_nadp))
