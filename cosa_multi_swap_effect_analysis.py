import cobra
import copy
import os
import numpy
import shutil
from typing import List
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from cosa_get_suffix import cosa_get_suffix
from fba import get_fba_base_problem, perform_fba_flux_maximization
from helper import ensure_folder_existence, json_load, json_write, json_zip_load, json_zip_write
from optmdfpathway import STANDARD_R, STANDARD_T, get_optmdfpathway_base_problem, get_thermodynamic_bottlenecks
from optimization import perform_variable_maximization
from cosa_create_table import create_cosa_tables
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from cosa_random_sampling_figures import create_cosa_figures, create_total_cosa_figure
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from cosa_add_promiscuity_constraints import cosa_add_promiscuity_constraints


def cosa_multi_swap_analysis(anaerobic: bool, selected_reactions: List[str]):
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=False)

    suffix = cosa_get_suffix(anaerobic, expanded=False)

    nadx_scenario = "WILDTYPE"
    old_cobra_model = copy.deepcopy(cobra_model)
    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"
    for concentration_scenario in ("STANDARDCONC", "VIVOCONC"):
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

        swap_results = {}
        swap_results[str(selected_reactions)] = {}

        for reaction in selected_reactions:
            if reaction.endswith("_ORIGINAL_NAD_TCOSA"):
                other_id = reaction.replace("_ORIGINAL_NAD_TCOSA", "_VARIANT_NADP_TCOSA")
            elif reaction.endswith("_ORIGINAL_NADP_TCOSA"):
                other_id = reaction.replace("_ORIGINAL_NADP_TCOSA", "_VARIANT_NAD_TCOSA")
            else:
                input(f"ERROR: No ORIGINAL TCOSA ID {reaction}")

            # Do swap
            optmdfpathway_base_variables[reaction].bounds(
                0.0,
                0.0,
            )
            optmdfpathway_base_variables[other_id].bounds(
                0.0,
                1000.0,
            )
            # end of Do swap

        for growth_rate_str in growth_rates:
            swap_results[str(selected_reactions)][growth_rate_str] = {}
            growth_rate = optmdf_json[growth_rate_str]["values"][biomass_reaction_id] - 1e-6

            print("Set growth to", growth_rate)
            optmdfpathway_base_variables[biomass_reaction_id].bounds(
                growth_rate,
                1e12
            )

            print(">OPTMDF calculations")
            optmdfpathway_result = perform_variable_maximization(
                optmdfpathway_base_problem,
                "var_B"
            )
            print(optmdfpathway_result["status"])
            if optmdfpathway_result["status"] != "Optimal":
                break
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
                break
            swapped_optsubmdf = optsubmdfpathway_result["values"]["var_B2"]
            original_optsubmdf = optsubmdf_json[growth_rate_str]["values"]["var_B2"]
            print("Swapped var_B2:", swapped_optsubmdf, "kJ/mol")
            print("Original var_B2: ", original_optsubmdf, "kJ/mol")
            optsubmdf_difference = swapped_optsubmdf - original_optsubmdf

            swap_results[str(selected_reactions)][growth_rate_str]["OptMDF"] = round(optmdf_difference, 6)
            swap_results[str(selected_reactions)][growth_rate_str]["OptSubMDF"] = round(optsubmdf_difference, 6)

            json_write(f"./cosa/results{suffix}/multi_swap_results_{concentration_scenario}.json", swap_results)

cosa_multi_swap_analysis(anaerobic=False, selected_reactions=["GND_ORIGINAL_NADP_TCOSA", "G6PDH2r_FWD_ORIGINAL_NADP_TCOSA", "G6PDH2r_REV_ORIGINAL_NADP_TCOSA"])
cosa_multi_swap_analysis(anaerobic=True, selected_reactions=["GND_ORIGINAL_NADP_TCOSA", "G6PDH2r_FWD_ORIGINAL_NADP_TCOSA", "G6PDH2r_REV_ORIGINAL_NADP_TCOSA"])
