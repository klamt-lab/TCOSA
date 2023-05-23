"""Contains the script for performing the single activation effect analysis."""

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
from optimization import perform_variable_maximization, perform_variable_minimization
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from cosa_add_promiscuity_constraints import cosa_add_promiscuity_constraints


# PUBLIC FUNCTIONS #
def cosa_single_activation_effect_analysis(anaerobic: bool, c_source: str = "glucose"):
    """Performs the single reaction activation effect analysis.

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
        concentration_scenarios = ("STANDARDCONC",) #"VIVOCONC",)

    for concentration_scenario in concentration_scenarios:
        activation_json_path = f"./cosa/results{suffix}/activation_results_{concentration_scenario}.json"

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
            add_optmdf_bottleneck_analysis=False,
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
            activation_results = json_load(activation_json_path)
        except FileNotFoundError:
            activation_results = {}
        counter = 0
        for reaction in cobra_model.reactions:
            if (reaction.id not in dG0_values.keys()):
                continue

            if reaction.id != "PYK":
                continue

            activation_results[reaction.id] = {}

            for growth_rate_str in growth_rates:
                growth_rate = float(growth_rate_str.replace(",", "."))

                activation_results[reaction.id][growth_rate_str] = {}

                print("Set growth to", growth_rate)
                optmdfpathway_base_variables[biomass_reaction_id].bounds(
                    growth_rate,
                    1e12
                )

                # Do activation
                optmdfpathway_base_variables[f"z_var_"+reaction.id].bounds(
                    1.0,
                    1.0,
                )
                # end of activation

                print(">OPTMDF calculations")
                optmdfpathway_result = perform_variable_maximization(
                    optmdfpathway_base_problem,
                    "var_B"
                )
                print(optmdfpathway_result["status"])
                if optmdfpathway_result["status"] != "Optimal":
                    activation_results[reaction.id][growth_rate_str]["OptMDF"] = float("NaN")
                    activation_results[reaction.id][growth_rate_str]["OptSubMDF"] = float("NaN")
                    continue
                print("Growth", optmdfpathway_result["values"][biomass_reaction_id])
                swapped_optmdf = optmdfpathway_result["values"]["var_B"]
                original_optmdf = optmdf_json[growth_rate_str]["values"]["var_B"]
                print("Swapped var_B:", swapped_optmdf, "kJ/mol")
                print("Original var_B: ", original_optmdf, "kJ/mol")
                optmdf_difference = swapped_optmdf - original_optmdf

                """
                if optmdf_difference < 0.0:
                    optmdfpathway_base_problem2 = get_optmdfpathway_base_problem(
                        cobra_model=cobra_model,
                        dG0_values=dG0_values,
                        metabolite_concentration_values=used_concentration_values,
                        ratio_constraint_data=ratio_constraint_data,
                        R=STANDARD_R,
                        T=STANDARD_T,
                        extra_constraints=[],
                        sub_network_ids=[],
                        add_optmdf_bottleneck_analysis=True,
                    )
                    print(">Get model variables dictionary")
                    optmdfpathway_base_variables2 = optmdfpathway_base_problem2.variablesDict()

                    print(">Set no promiscuity constraint")
                    optmdfpathway_base_problem2 = cosa_add_promiscuity_constraints(
                        optmdfpathway_base_problem=optmdfpathway_base_problem2,
                        optmdfpathway_base_variables=optmdfpathway_base_variables2,
                        cobra_model=cobra_model,
                        dG0_values=dG0_values,
                    )

                    optmdfpathway_base_variables2 = optmdfpathway_base_problem2.variablesDict()
                    # Do activation
                    optmdfpathway_base_variables2[f"z_var_"+reaction.id].bounds(
                        1.0,
                        1.0,
                    )
                    # end of activation

                    optmdfpathway_base_variables2[biomass_reaction_id].bounds(
                        growth_rate,
                        1e12
                    )
                    optmdfpathway_base_variables2["var_B"].bounds(
                        original_optmdf,
                        1e12
                    )
                    optmdfpathway_result = perform_variable_minimization(
                        optmdfpathway_base_problem2,
                        "zb_sum_var"
                    )
                    print("Status:", optmdfpathway_result["status"])
                    print(
                        f"Σ of reaction changes to achieve OptMDF of >= {MIN_OPTMDF} kJ/mol (zb_sum):",
                        optmdfpathway_result["values"]["zb_sum_var"],
                        "reaction changes"
                    )
                    print("Reached MDF (lower bound for OptMDF):", optmdfpathway_result["values"]["var_B"], "kJ/mol")

                    print(f"->LIST OF FOUND BOTTLENECK CORRECTIONS FOR {nadx_scenario}:")
                    for key in optmdfpathway_result["values"].keys():
                        dG0_change = optmdfpathway_result["values"][key]
                        if key.startswith("zb_var") and (dG0_change > 1e-3):
                            reaction_id = key.replace('zb_var_', '')
                            text = f"{reaction_id}: {dG0_change} kJ/mol"
                            print(text)
                """

                print(">SubMDF calculations")
                optmdfpathway_base_variables["var_B"].bounds(MIN_OPTMDF, 1e6)
                optsubmdfpathway_result = perform_variable_maximization(
                    optmdfpathway_base_problem,
                    "var_B2"
                )
                print(optsubmdfpathway_result["status"])
                if optsubmdfpathway_result["status"] != "Optimal":
                    activation_results[reaction.id][growth_rate_str]["OptSubMDF"] = float("NaN")
                    continue
                swapped_optsubmdf = optsubmdfpathway_result["values"]["var_B2"]
                original_optsubmdf = optsubmdf_json[growth_rate_str]["values"]["var_B2"]
                print("Swapped var_B2:", swapped_optsubmdf, "kJ/mol")
                print("Original var_B2: ", original_optsubmdf, "kJ/mol")
                optsubmdf_difference = swapped_optsubmdf - original_optsubmdf

                activation_results[reaction.id][growth_rate_str]["OptMDF"] = round(optmdf_difference, 6)
                activation_results[reaction.id][growth_rate_str]["OptSubMDF"] = round(optsubmdf_difference, 6)

                if counter >= 10:
                    json_write(activation_json_path, activation_results)
                    counter = 0
                counter += 1

                # Undo activation
                optmdfpathway_base_variables[f"z_var_"+reaction.id].bounds(
                    0.0,
                    1.0,
                )
                # end of undo activation
        json_write(activation_json_path, activation_results)

cosa_single_activation_effect_analysis(anaerobic=False)
