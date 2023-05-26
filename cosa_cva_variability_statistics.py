from helper import json_load

top = 15
used_mu = "0,818"

for aerobicity in ("aerobic",):
    for target in ("OPTMDF", "OPTSUBMDF"):
        txtpath = f"cva_vars_{aerobicity}_{target}_{used_mu}_STANDARDCONC.txt"
        jsondata = json_load(f"./cosa/results_{aerobicity}/cva_{target}_STANDARDCONC.json")

        num_not_max_var = 0
        var_to_mets = {}
        mets_to_minmax = {}
        for x_var_id in jsondata.keys():
            met_id = x_var_id.replace("x_var_", "")
            min_conc = round(jsondata[x_var_id][used_mu]["min"], 7)
            max_conc = round(jsondata[x_var_id][used_mu]["max"], 7)

            mets_to_minmax[met_id] = (min_conc, max_conc)

            variability = max_conc - min_conc
            if variability not in var_to_mets.keys():
                var_to_mets[variability] = []

            if variability < 0.02-1e-6:
                num_not_max_var += 1

            var_to_mets[variability].append(met_id)

        sorted_vars = list(var_to_mets.keys())
        sorted_vars.sort()

        txt = ""
        txt += f"Number metabolites with concentration range smaller than [1e-6; 0.02]: {num_not_max_var}\n\n"

        counter = 1
        txt += "TOP OF THE LOW VARIABILITIES:\n"
        for sorted_var in sorted_vars:
            met_ids = var_to_mets[sorted_var]
            txt += f"#{counter}: {sorted_var}\n"
            for met_id in met_ids:
                txt += f" >{met_id}: [{mets_to_minmax[met_id][0]}; {mets_to_minmax[met_id][1]}]\n"

            counter += 1
            if counter > top:
                break

        with open(txtpath, "w", encoding="utf-8") as f:
            f.write(txt)

