import cobra
from helper import json_zip_load

searched_metabolite_id = "nadp_tcosa_c"
model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")
data = json_zip_load("cosa/results_aerobic/runs/OPTMDF_STANDARDCONC_FLEXIBLE.json")

metabolite = model.metabolites.get_by_id(searched_metabolite_id)
production = []
consumption = []
reaction_ids = [x.id for x in model.reactions]
for growth_rate in list(data.keys())[:1]:
    growth_data = data[growth_rate]["values"]
    for key, value in growth_data.items():
        if key not in reaction_ids:
            continue
        if value <= 1e-4:
            continue
        metabolide_ids = [x.id for x in model.reactions.get_by_id(key).metabolites]
        if sum([searched_metabolite_id in x for x in metabolide_ids]) == 0:
            continue
        if model.reactions.get_by_id(key).metabolites[metabolite] > 0.0:
            production.append((key, value))
        else:
            consumption.append((key, value))

output_str = ""
output_str += "PRODUCTION\n"
for data in production:
    output_str += f"{data[0]} {data[1]}\n"
output_str += "\nCONSUMPTION\n"
for data in consumption:
    output_str += f"{data[0]} {data[1]}\n"
with open("./scenario_viewer.txt", "w", encoding="utf-8") as f:
    f.write(output_str)
