import cobra
from helper import json_load, json_write, pickle_load


# Bennett 2009
with open("resources/in_vivo_concentration_data/2009_Bennet_full_raw_data.txt", encoding="utf-8") as f:
    lines = f.readlines()
lines = [x.replace("\n", "") for x in lines]

names = [x.split(";")[0].split(" (")[0].split("  ")[:-1][0].lower() for x in lines]

print(names)

model = cobra.io.read_sbml_model("resources/iML1515_irreversible_cleaned.xml")
name_to_bigg_ids = {}
for name in names:
    for metabolite in model.metabolites:
        lower_metname = metabolite.name.lower()
        if (name == lower_metname) or (f"l-{name}" == lower_metname) or (f"(s)-{name}" == lower_metname) or (f"d-{name}" == lower_metname) or (name.split(" ")[0] == lower_metname.split(" ")[0]):
            name_to_bigg_ids[name] = metabolite.id
            break
    if name not in name_to_bigg_ids.keys():
        name_to_bigg_ids[name] = "N/A"

# Manual setting for found other metabolites
name_to_bigg_ids["2,3-dihydroxybenzoic acid"] = "23dhb_c"  # KEGG ID C00196
name_to_bigg_ids["3-phosphoglycerate"] = "3pg_c"  # KEGG ID C00197
name_to_bigg_ids["6-phosphogluconate"] = "6pgc_c"  # KEGG ID C00345
name_to_bigg_ids["acetylphosphate"] = "actp_c"  # KEGG ID C00227
name_to_bigg_ids["adp-glucose"] = "adpglc_c"  # KEGG ID C00498
name_to_bigg_ids["alpha-ketoglutarate"] = "akg_c"  # KEGG ID C00026
name_to_bigg_ids["carbamylaspartate"] = "cbasp_c"  # KEGG ID C00438
name_to_bigg_ids["deoxyribose-5-p"] = "2dr5p_c"  # KEGG ID C00673
name_to_bigg_ids["fad"] = "fad_c"  # KEGG ID C00016
name_to_bigg_ids["fructose-1,6-bisphosphate"] = "fdp_c"  # MetaNetX MNXM417
name_to_bigg_ids["glucosamine-6-phosphate"] = "gam6p_c"  # MetaNetX MNXM370
name_to_bigg_ids["glycerate"] = "glyc__R"
name_to_bigg_ids["glycerol-3-phosphate"] = "glyc3p_c"
name_to_bigg_ids["n-acetyl-glucosamine-1p"] = "acgam1p_c"
name_to_bigg_ids["n-acetyl-ornithine"] = "acorn_c"
name_to_bigg_ids["propionyl-coa"] = "ppcoa_c"
name_to_bigg_ids["prpp"] = "prpp_c"
name_to_bigg_ids["udp-glucuronate"] = "udpglcur_c"
name_to_bigg_ids["udp-glucose"] = "udpg_c"
# gluconolactone does not seem to exist in iML1515
# hexose-pa is a general term
# pentose-pd is a general term
# NAD(P)(H) are kept free

for name in name_to_bigg_ids.keys():
    if name_to_bigg_ids[name] == "N/A":
        print("MISSING:", name)

json_write("./resources/in_vivo_concentration_data/2009_Bennett_names.json", name_to_bigg_ids)

name_dict = json_load("./resources/in_vivo_concentration_data/2009_Bennett_names.json")

with open("./resources/in_vivo_concentration_data/2009_Bennet_full_raw_data.txt", encoding="utf-8") as f:
    lines = f.readlines()
lines = [x.replace("\n", "") for x in lines]

lbs = []
ubs = []
full_data_dict = {}
for line in lines:
    linesplit = line.split("  ")
    name = linesplit[0]
    bigg_id = name_dict[name.lower()]

    bigg_id = (bigg_id+"\b").replace("_e\b", "_c")
    bigg_id = (bigg_id+"\b").replace("_p\b", "_c")
    bigg_id = bigg_id.replace("\b", "")

    concentration_raw = linesplit[1]
    mean_conc = float(concentration_raw.split(" ")[0])
    lb_conc = float(concentration_raw.split(" ")[1].replace("(", ""))
    ub_conc = float(concentration_raw.split(" ")[3].replace(")", ""))
    full_data_dict[bigg_id] = {
        "lb": lb_conc,
        "ub": ub_conc,
        "mean": mean_conc,
        "name": name,
    }
    lbs.append(lb_conc)
    ubs.append(ub_conc)
    if ub_conc > 0.02:
        print("OVER", full_data_dict[bigg_id])
    if lb_conc < 1e-6:
        print("UNDER", full_data_dict[bigg_id])

print("MIN", min(lbs))
print("MAX", max(ubs))

# Glutathione and glutathione disulfide special treatment
# Gluthathione disulfide seems to be mixed together with oxidized/reduced glutathione in iML1515
# Therefore, only an upper bound for all affected metabolites is set at the sum of the measured
# upper bound of gluthathione and glutathione sulfate concentrations.
# Since the distribution of oxidized and reduced glutathione remains unclear, the default lower
# bound is set.
ub_glutathione = 1.79E-2
ub_glutathione_disulfide = 2.90E-3
sum_ubs_glutathiones = ub_glutathione + ub_glutathione_disulfide
# -> For oxidized glutathione (gthox_c)
full_data_dict["gthox_c"] = {
    "lb": 1e-6,
    "ub": sum_ubs_glutathiones,
    "mean": (1e-6+sum_ubs_glutathiones)/2,
    "name": name,
}
# -> For reduced glutathione (gthrd_c)
full_data_dict["gthrd_c"] = {
    "lb": 1e-6,
    "ub": sum_ubs_glutathiones,
    "mean": (1e-6+sum_ubs_glutathiones)/2,
    "name": name,
}

json_write("./resources/in_vivo_concentration_data/2009_Bennet_full_data.json", full_data_dict)
