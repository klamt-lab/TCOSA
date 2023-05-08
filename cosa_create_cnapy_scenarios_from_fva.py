import cobra
from helper import json_load, json_write, ensure_folder_existence, json_zip_write, pickle_load, get_files

model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")
reaction_ids = [x.id for x in model.reactions]

files = get_files("./cosa")

for anaerobic in (True, False):
    aerobicity_string = "anaerobic" if anaerobic else "aerobic"
    basepath  = f"./cosa/"
    scenpath = "./cosa_cnapy_scenarios/"
    ensure_folder_existence(scenpath)

    variability_files = [
        basepath+x for x in files
        if x.startswith("variability__"+aerobicity_string) and ("WILDTYPE" in x) and ("STANDARDCONCS" in x) and ("STOICHIO" not in x)
    ]

    zip_file_index = 0
    for variability_file in variability_files:
        json_data = json_load(variability_file)
        growth_rate = variability_file.split("_TEST_")[1].split("_STANDARD")[0]
        for target in ("OptMDF", "OptSubMDF"):
            scen_data = {}
            for json_datapoint in json_data[target]:
                reaction_id = json_datapoint[0]
                base_id = reaction_id.replace("_ORIGINAL_NADP_TCOSA", "")\
                    .replace("_ORIGINAL_NAD_TCOSA", "")\
                    .replace("_VARIANT_NAD_TCOSA", "")\
                    .replace("_VARIANT_NADP_TCOSA", "")\
                    .replace("_FWD", "")\
                    .replace("_REV", "")
                if base_id not in scen_data:
                    scen_data[base_id] = [0.0, 0.0]
                abs_flux_low = abs(json_datapoint[1])
                abs_flux_high = abs(json_datapoint[2])

                if "_REV" in reaction_id:
                    scen_data[base_id][0] -= abs_flux_high
                    scen_data[base_id][1] -= abs_flux_low
                elif "_FWD" in reaction_id:
                    scen_data[base_id][0] += abs_flux_low
                    scen_data[base_id][1] += abs_flux_high
                else:
                    scen_data[base_id][0] += abs_flux_low
                    scen_data[base_id][1] += abs_flux_high
            json_write(
                scenpath+growth_rate+"_STANDARDCONCS_"+f"_{aerobicity_string}_{target}.scen",
                scen_data
            )
        zip_file_index += 1
