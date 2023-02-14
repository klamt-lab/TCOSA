import json
from tokenize import cookie_re
from helper import json_load
import cobra
from math import exp

dG0_values = json_load("cosa/dG0_values.json")
cobra_model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")

def print_reacs(reacs):
    for reac in reacs:
        cobra_reac = cobra_model.reactions.get_by_id(reac)
        dG0 = round(dG0_values[reac]["dG0"], 2)
        print(f"{reac}; {cobra_reac.name}; {cobra_reac.reaction}; {dG0}")

print("Bottlenecks: standardconcs")
reacs = [
    "KDOCT2",
    "MECDPS",
    "DHPPDA2",
    "ATPPRT",
    "IG3PS",
    "MCTP1App",
    "MALCOAMT",
    "AIRC3_REV",
    "AIRC3_FWD",
    "SHCHD2_ORIGINAL_NAD_TCOSA",
    "SHCHD2_VARIANT_NADP_TCOSA",
]
print_reacs(reacs)
print("\n~~~~\n~~~~\n")
print("Bottlenecks: paperconcs, aerobic")
reacs = [
    "KDOCT2",
    "MECDPS",
    "DHPPDA2",
    "ATPPRT",
    "IG3PS",
    "MCTP1App",
    "MALCOAMT",
    "AIRC3_REV",
    "AIRC3_FWD",
    "SHCHD2_ORIGINAL_NAD_TCOSA",
    "SHCHD2_VARIANT_NADP_TCOSA",
    "ASPTA_REV",
    "ASPTA_FWD",
    "ASAD_REV_VARIANT_NAD_TCOSA",
    "ASAD_REV_ORIGINAL_NADP_TCOSA",
    "ASAD_FWD_VARIANT_NAD_TCOSA",
    "ASAD_FWD_ORIGINAL_NADP_TCOSA",
]
print_reacs(reacs)
print("\n~~~~\n~~~~\n")
print("Bottlenecks: paperconcs, anaerobic")
reacs = [
    "KDOCT2",
    "MECDPS",
    "DHPPDA2",
    "ATPPRT",
    "IG3PS",
    "MCTP1App",
    "MALCOAMT",
    "AIRC3_REV",
    "AIRC3_FWD",
    "SHCHD2_ORIGINAL_NAD_TCOSA",
    "SHCHD2_VARIANT_NADP_TCOSA",
    "ASPK_FWD",
    "ASPK_REV",
    "ACACT2r_FWD",
    "ACACT2r_REV",
    "ACACT4r_FWD",
    "ACACT4r_REV",
    "ACACT6r_FWD",
    "ACACT6r_REV",
    "ACACT1r_FWD",
    "ACACT1r_REV",
]
print_reacs(reacs)
