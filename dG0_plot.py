from helper import json_load, save_xy_point_plot
from typing import Any, Dict
import math

dgdata: Dict[str, Any] = json_load("resources/dG0_iML1515_irreversible_cleaned.json")

dG0s = []
uncertainties = []
for reaction_id, data in dgdata.items():
    dG0 = data["dG0"]
    uncertainty = abs(data["uncertainty"])
    compartments = data["num_compartments"]

    # if (reaction_id in ["SHCHD2", "IG3PS", "AIRC3_FWD", "AIRC3_REV", "MCTP1App", "ATPPRT", "DHPPDA2", "MECDPS", "KDOCT2", "MALCOAMT"]):
    if uncertainty > 10:
        print(f"{reaction_id} | {dG0} ± {uncertainty} kJ/mol | {compartments} compartment(s)")

    dG0s.append(dG0)
    uncertainties.append(uncertainty)

save_xy_point_plot(
    path="./uncertainty_plot.png",
    xs=dG0s,
    ys=uncertainties,
    title="ΔG'° vs. uncertainty",
    xlabel="ΔG'° [kJ/mol]",
    ylabel="uncertainty [kJ/mol]"
)
