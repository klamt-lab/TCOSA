"""Contains the script for performing the single cofactor swap analysis."""

# IMPORTS #
# External
import cobra
import copy
import os
import numpy
import shutil
from typing import List
# Internal
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from cosa_get_suffix import cosa_get_suffix
from fba import get_fba_base_problem, perform_fba_flux_maximization
from helper import ensure_folder_existence, json_load, json_write, json_zip_load, json_zip_write
from optmdfpathway import STANDARD_R, STANDARD_T, get_optmdfpathway_base_problem
from optimization import perform_variable_maximization
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from cosa_add_promiscuity_constraints import cosa_add_promiscuity_constraints


# PUBLIC FUNCTIONS #
def cosa_single_swap_analysis(anaerobic: bool, c_source: str = "glucose"):
    """Performs the single redox cofactor swap analysis.

    Args:
        anaerobic (bool): Is it anaerobic (True)?
        c_source (str, optional): Either 'glucose' or 'acetate'. Defaults to "glucose".
    """
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=False, c_source=c_source)

    suffix = cosa_get_suffix(anaerobic, expanded=False, c_source=c_source)

    nadx_scenario = "WILDTYPE"
    old_cobra_model = copy.deepcopy(cobra_model)
    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    if (c_source != "glucose") or (anaerobic):
        concentration_scenarios = ("STANDARDCONC",)
    else:
        concentration_scenarios = ("STANDARDCONC", "VIVOCONC",)

    for concentration_scenario in concentration_scenarios:
        swap_json_path = f"./cosa/results{suffix}/swap_results_{concentration_scenario}.json"

        if concentration_scenario == "STANDARDCONC":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
            used_concentration_values = concentration_values_free
        elif concentration_scenario == "VIVOCONC":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper

        print("~~~")
        print(nadx_scenario)

        optmdf_json_path = f"./cosa/results{suffix}/runs/OPTMDF_{concentration_scenario}_{nadx_scenario}.json"
        optsubmdf_json_path = f"./cosa/results{suffix}/runs/OPTSUBMDF_{concentration_scenario}_{nadx_scenario}.json"

        optmdf_json = json_zip_load(optmdf_json_path)
        optsubmdf_json = json_zip_load(optsubmdf_json_path)
        growth_rates = optmdf_json.keys()

        cobra_model = copy.deepcopy(old_cobra_model)
        cobra_model = cosa_get_model_with_nadx_scenario(
            nadx_scenario=nadx_scenario,
            cobra_model=cobra_model,
        )

        print(">Get base OptMDFpathway MILP...")
        optmdfpathway_base_problem = get_optmdfpathway_base_problem(
            cobra_model=cobra_model,
            dG0_values=dG0_values,
            metabolite_concentration_values=used_concentration_values,
            ratio_constraint_data=ratio_constraint_data,
            R=STANDARD_R,
            T=STANDARD_T,
            extra_constraints=[],
            sub_network_ids=get_all_tcosa_reaction_ids(cobra_model),
        )
        print(">Get model variables dictionary")
        optmdfpathway_base_variables = optmdfpathway_base_problem.variablesDict()

        print(">Set no promiscuity constraint")
        optmdfpathway_base_problem = cosa_add_promiscuity_constraints(
            optmdfpathway_base_problem=optmdfpathway_base_problem,
            optmdfpathway_base_variables=optmdfpathway_base_variables,
            cobra_model=cobra_model,
            dG0_values=dG0_values,
        )
        try:
            swap_results = json_load(swap_json_path)
        except FileNotFoundError:
            swap_results = {}
        counter = 0
        for reaction in cobra_model.reactions:
            if (reaction.id not in dG0_values.keys()):
                continue
            if ((reaction.id != "ICDHyr_FWD_ORIGINAL_NADP_TCOSA") and (reaction.id != "PDH_ORIGINAL_NAD_TCOSA")) or (c_source != "acetate"):
                if (reaction.id in swap_results.keys()):
                    continue
                multiplier = 1.0
            else:
                if (reaction.id == "ICDHyr_FWD_ORIGINAL_NADP_TCOSA"):
                    multiplier = .9875 # Numeric error at highest growth rate D:
                else:
                    multiplier = .999


            if reaction.id.endswith("_ORIGINAL_NAD_TCOSA"):
                other_id = reaction.id.replace("_ORIGINAL_NAD_TCOSA", "_VARIANT_NADP_TCOSA")
            elif reaction.id.endswith("_VARIANT_NAD_TCOSA"):
                continue
            elif reaction.id.endswith("_ORIGINAL_NADP_TCOSA"):
                other_id = reaction.id.replace("_ORIGINAL_NADP_TCOSA", "_VARIANT_NAD_TCOSA")
            elif reaction.id.endswith("_VARIANT_NADP_TCOSA"):
                continue
            else:
                continue

            original_real_ub = reaction.upper_bound
            original_other_ub = cobra_model.reactions.get_by_id(other_id).upper_bound

            if original_other_ub > 0.0:
                input(f"ERROR WITH UB in reaction {reaction.id}")

            swap_results[reaction.id] = {}

            for growth_rate_str in growth_rates:
                growth_rate = float(growth_rate_str.replace(",", "."))

                if (not anaerobic) and (growth_rate > 0.4):
                    continue
                swap_results[reaction.id][growth_rate_str] = {}

                print("Set growth to", growth_rate)
                optmdfpathway_base_variables[biomass_reaction_id].bounds(
                    growth_rate*multiplier,
                    1e12
                )

                # Do swap
                optmdfpathway_base_variables[reaction.id].bounds(
                    0.0,
                    0.0,
                )
                optmdfpathway_base_variables[other_id].bounds(
                    0.0,
                    reaction.upper_bound,
                )
                # end of Do swap

                print(">OPTMDF calculations")
                optmdfpathway_result = perform_variable_maximization(
                    optmdfpathway_base_problem,
                    "var_B"
                )
                print(optmdfpathway_result["status"])
                if optmdfpathway_result["status"] != "Optimal":
                    swap_results[reaction.id][growth_rate_str]["OptMDF"] = float("NaN")
                    swap_results[reaction.id][growth_rate_str]["OptSubMDF"] = float("NaN")
                    continue
                print("Growth", optmdfpathway_result["values"][biomass_reaction_id])
                swapped_optmdf = optmdfpathway_result["values"]["var_B"]
                original_optmdf = optmdf_json[growth_rate_str]["values"]["var_B"]
                print("Swapped var_B:", swapped_optmdf, "kJ/mol")
                print("Original var_B: ", original_optmdf, "kJ/mol")
                optmdf_difference = swapped_optmdf - original_optmdf


                print(">SubMDF calculations")
                optmdfpathway_base_variables["var_B"].bounds(MIN_OPTMDF, 1e6)
                optsubmdfpathway_result = perform_variable_maximization(
                    optmdfpathway_base_problem,
                    "var_B2"
                )
                print(optsubmdfpathway_result["status"])
                if optsubmdfpathway_result["status"] != "Optimal":
                    swap_results[reaction.id][growth_rate_str]["OptSubMDF"] = float("NaN")
                    continue
                swapped_optsubmdf = optsubmdfpathway_result["values"]["var_B2"]
                original_optsubmdf = optsubmdf_json[growth_rate_str]["values"]["var_B2"]
                print("Swapped var_B2:", swapped_optsubmdf, "kJ/mol")
                print("Original var_B2: ", original_optsubmdf, "kJ/mol")
                optsubmdf_difference = swapped_optsubmdf - original_optsubmdf

                swap_results[reaction.id][growth_rate_str]["OptMDF"] = round(optmdf_difference, 6)
                swap_results[reaction.id][growth_rate_str]["OptSubMDF"] = round(optsubmdf_difference, 6)
                # shutil.copyfile(f"./cosa/results{suffix}/swap_results.json", f"./cosa/results{suffix}/swap_results2.json")

                if counter >= 10:
                    json_write(swap_json_path, swap_results)
                    counter = 0
                counter += 1

                # Undo swap
                optmdfpathway_base_variables[reaction.id].bounds(
                    0.0,
                    original_real_ub,
                )
                optmdfpathway_base_variables[other_id].bounds(
                    0.0,
                    0.0,
                )
                # end of Undo swap
        json_write(swap_json_path, swap_results)
