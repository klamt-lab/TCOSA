"""This script shows, when running, some statistics about the ΔG'° values collected for iML1515."""


# IMPORT SECTION #
# External imports
import cobra
from statistics import mean, median
# Internal imports
from helper import json_load


# ACTUAL LOGIC SECTION #
model = cobra.io.read_sbml_model("resources/iML1515_irreversible_cleaned.xml")
dG0s = json_load("resources/dG0_iML1515_irreversible_cleaned.json")

with open("./deleted_transporters.txt", "r", encoding="utf-8") as f:
    deleted = f.readlines()
deleted = [x.split(" | ")[0].replace("_ORIGINAL_NADP_TCOSA", "")\
                    .replace("_ORIGINAL_NAD_TCOSA", "")\
                    .replace("_VARIANT_NAD_TCOSA", "")\
                    .replace("_VARIANT_NADP_TCOSA", "")\
                    .replace("_FWD", "")\
                    .replace("_REV", "") for x in deleted]

num_with_dG0 = 0
for key in list(dG0s.keys()):
    if key.endswith("_REV"):
        continue
    key = key.replace("_ORIGINAL_NADP_TCOSA", "")\
            .replace("_ORIGINAL_NAD_TCOSA", "")\
            .replace("_VARIANT_NAD_TCOSA", "")\
            .replace("_VARIANT_NADP_TCOSA", "")\
            .replace("_FWD", "")\
            .replace("_REV", "")
    if key in deleted:
        continue
    num_with_dG0 += 1

print("num_with_dG0", num_with_dG0)
