import cobra
from helper import json_zip_load, json_write, ensure_folder_existence, json_zip_write, pickle_load

model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")
reaction_ids = [x.id for x in model.reactions]

for anaerobic in (True, False):
    aerobicity_string = "anaerobic" if anaerobic else "aerobic"
    basepath  = f"cosa/results_{aerobicity_string}/runs/"
    scenpath = "./cosa_cnapy_scenarios/"
    ensure_folder_existence(scenpath)

    if anaerobic:
        growth_rates = ("0,371",)
    else:
        growth_rates = ("0,868",)

    zip_files = [
        "OPTMDF_STANDARDCONC_WILDTYPE.json",
        # "OPTMDF_STANDARDCONC_FLEXIBLE.json",
        # "OPTMDF_STANDARDCONC_SINGLE_COFACTOR.json",
        # "OPTMDF_STANDARDCONC_RANDOMFIXED_0.json",
        # "OPTMDF_STANDARDCONC_RANDOMS_0.json",
        # "OPTMDF_VIVOCONC_WILDTYPE.json",
        # "OPTMDF_VIVOCONC_FLEXIBLE.json",
        # "OPTMDF_VIVOCONC_SINGLE_COFACTOR.json",
        # "OPTMDF_VIVOCONC_RANDOMFIXED_0.json",
        # "OPTMDF_VIVOCONC_RANDOMS_0.json",

        "OPTSUBMDF_STANDARDCONC_WILDTYPE.json",
        # "OPTSUBMDF_STANDARDCONC_FLEXIBLE.json",
        # "OPTSUBMDF_STANDARDCONC_SINGLE_COFACTOR.json",
        # "OPTSUBMDF_STANDARDCONC_RANDOMFIXED_0.json",
        # "OPTSUBMDF_STANDARDCONC_RANDOMS_0.json",
        # "OPTSUBMDF_VIVOCONC_WILDTYPE.json",
        # "OPTSUBMDF_VIVOCONC_FLEXIBLE.json",
        # "OPTSUBMDF_VIVOCONC_SINGLE_COFACTOR.json",
        # "OPTSUBMDF_VIVOCONC_RANDOMFIXED_0.json",
        # "OPTSUBMDF_VIVOCONC_RANDOMS_0.json",
    ]
    zip_paths = [basepath+x for x in zip_files]

    zip_file_index = 0
    for zip_path in zip_paths:
        zip_data = json_zip_load(zip_path)
        for growth_rate in growth_rates:
            scen_data = {}
            for reaction_id in reaction_ids:
                base_id = reaction_id.replace("_ORIGINAL_NADP_TCOSA", "")\
                    .replace("_ORIGINAL_NAD_TCOSA", "")\
                    .replace("_VARIANT_NAD_TCOSA", "")\
                    .replace("_VARIANT_NADP_TCOSA", "")\
                    .replace("_FWD", "")\
                    .replace("_REV", "")
                if base_id not in scen_data:
                    scen_data[base_id] = [0.0, 0.0]
                abs_flux = abs(zip_data[growth_rate]["values"][reaction_id])
                scen_data[base_id][0] += abs_flux
                scen_data[base_id][1] += abs_flux
            json_write(
                scenpath+zip_files[zip_file_index].replace(".json", "")+"_"+growth_rate.replace(",", "_")+f"_{aerobicity_string}.scen",
                scen_data
            )
        zip_file_index += 1
