import cobra
import copy
import os
import numpy
import pulp
from typing import Any, Dict, List
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from cosa_get_suffix import cosa_get_suffix
from fba import get_fba_base_problem, perform_fba_flux_maximization
from helper import ensure_folder_existence, json_write, json_zip_write
from optmdfpathway import STANDARD_R, STANDARD_T
from multoptmdfpathway import get_multoptmdfpathway_base_problem
from optimization import perform_variable_maximization, perform_variable_minimization, solve_current_problem, set_timelimit
from cosa_create_table import create_cosa_tables
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from cosa_random_sampling_figures import create_cosa_figures, create_total_cosa_figure
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from cosa_add_promiscuity_constraints import cosa_add_promiscuity_constraints
import sys

set_timelimit(180)

def cosa_multi_flux_sampling(num_samplings: int):
    expanded = False

    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=False, expanded=expanded)
    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    suffix = cosa_get_suffix(False, expanded)

    ensure_folder_existence("./cosa")
    ensure_folder_existence(f"./cosa/results{suffix}")
    ensure_folder_existence(f"./cosa/results{suffix}/mult_flux_samplings4")

    epsilon = 0.001
    conditions = [
        # [(biomass_reaction_id, 0.868-epsilon, float("inf")), ("var_B", 5.096-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 10.0)],  # OPTMDF_AEROBIC_0_868
        # [(biomass_reaction_id, 0.818-epsilon, float("inf")), ("var_B", 7.303-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 10.0)],  # OPTMDF_AEROBIC_0_818
        # [(biomass_reaction_id, 0.718-epsilon, float("inf")), ("var_B", 7.830-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 10.0)],  # OPTMDF_AEROBIC_0_718
        # [(biomass_reaction_id, 0.668-epsilon, float("inf")), ("var_B", 7.966-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 10.0)],  # OPTMDF_AEROBIC_0_668

        # [(biomass_reaction_id, 0.371-epsilon, float("inf")), ("var_B", 4.948-epsilon, float("inf")), ("EX_o2_e_REV", 0.0, 0.0), ("EX_glc__D_e_REV", 0.0, 20.0)],  # OPTMDF_ANAEROBIC_0_371
        # [(biomass_reaction_id, 0.321-epsilon, float("inf")), ("var_B", 6.998-epsilon, float("inf")), ("EX_o2_e_REV", 0.0, 0.0), ("EX_glc__D_e_REV", 0.0, 20.0)],  # OPTMDF_ANAEROBIC_0_321
        # [(biomass_reaction_id, 0.271-epsilon, float("inf")), ("var_B", 7.830-epsilon, float("inf")), ("EX_o2_e_REV", 0.0, 0.0), ("EX_glc__D_e_REV", 0.0, 20.0)],  # OPTMDF_ANAEROBIC_0_271
        # [(biomass_reaction_id, 0.221-epsilon, float("inf")), ("var_B", 7.920-epsilon, float("inf")), ("EX_o2_e_REV", 0.0, 0.0), ("EX_glc__D_e_REV", 0.0, 20.0)],  # OPTMDF_ANAEROBIC_0_221
        # [(biomass_reaction_id, 0.171-epsilon, float("inf")), ("var_B", 7.966-epsilon, float("inf")), ("EX_o2_e_REV", 0.0, 0.0), ("EX_glc__D_e_REV", 0.0, 20.0)],  # OPTMDF_ANAEROBIC_0_171

        [(biomass_reaction_id, 0.868-epsilon, float("inf")), ("var_B", MIN_OPTMDF, float("inf")), ("var_B2", 8.8760-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 10.0)],  # OPTSUBMDF_AEROBIC_0_868
        # [(biomass_reaction_id, 0.818-epsilon, float("inf")), ("var_B", MIN_OPTMDF, float("inf")), ("var_B2", 14.507-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 10.0)],  # OPTSUBMDF_AEROBIC_0_818
        # [(biomass_reaction_id, 0.718-epsilon, float("inf")), ("var_B", MIN_OPTMDF, float("inf")), ("var_B2", 26.804-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 10.0)],  # OPTSUBMDF_AEROBIC_0_718
        # [(biomass_reaction_id, 0.668-epsilon, float("inf")), ("var_B", MIN_OPTMDF, float("inf")), ("var_B2", 28.740-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 10.0)],  # OPTSUBMDF_AEROBIC_0_668
        [(biomass_reaction_id, 0.568-epsilon, float("inf")), ("var_B", MIN_OPTMDF, float("inf")), ("var_B2", 29.014-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 10.0)],  # OPTSUBMDF_AEROBIC_0_568

        [(biomass_reaction_id, 0.371-epsilon, float("inf")), ("var_B", MIN_OPTMDF, float("inf")), ("var_B2", 10.001-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 20.0)],  # OPTSUBMDF_ANAEROBIC_0_371
        # [(biomass_reaction_id, 0.321-epsilon, float("inf")), ("var_B", MIN_OPTMDF, float("inf")), ("var_B2", 16.195-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 20.0)],  # OPTSUBMDF_ANAEROBIC_0_321
        # [(biomass_reaction_id, 0.271-epsilon, float("inf")), ("var_B", MIN_OPTMDF, float("inf")), ("var_B2", 26.804-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 20.0)],  # OPTSUBMDF_ANAEROBIC_0_271
        [(biomass_reaction_id, 0.071-epsilon, float("inf")), ("var_B", MIN_OPTMDF, float("inf")), ("var_B2", 28.740-epsilon, float("inf")), ("EX_glc__D_e_REV", 0.0, 20.0)],  # OPTSUBMDF_ANAEROBIC_0_071
    ]

    condition_names = [
        "WT_OPTSUBMDF_AEROBIC_0_868",
        "WT_OPTSUBMDF_AEROBIC_0_818",
        "WT_OPTSUBMDF_AEROBIC_0_568",
        "WT_OPTSUBMDF_ANAEROBIC_0_371",
        # "WT_OPTSUBMDF_ANAEROBIC_0_321",
        "WT_OPTSUBMDF_ANAEROBIC_0_071",
        # "OPTMDF_AEROBIC_0_868",
        # "OPTMDF_AEROBIC_0_818",
        # "OPTMDF_AEROBIC_0_718",
        # "OPTMDF_AEROBIC_0_668",
        # "OPTMDF_ANAEROBIC_0_371",
        # "OPTMDF_ANAEROBIC_0_321",
        # "OPTMDF_ANAEROBIC_0_271",
        # "OPTMDF_ANAEROBIC_0_221",
        # "OPTMDF_ANAEROBIC_0_171",
        # "OPTSUBMDF_AEROBIC_0_868",
        # "OPTSUBMDF_AEROBIC_0_818",
        # "OPTSUBMDF_AEROBIC_0_718",
        # "OPTSUBMDF_AEROBIC_0_668",
        # "OPTSUBMDF_AEROBIC_0_568",
        # "OPTSUBMDF_ANAEROBIC_0_371",
        # "OPTSUBMDF_ANAEROBIC_0_321",
        # "OPTSUBMDF_ANAEROBIC_0_271",
        # "OPTSUBMDF_ANAEROBIC_0_071",
    ]
    old_cobra_model = copy.deepcopy(cobra_model)
    concentration_scenario = "STANDARDCONC"
    if concentration_scenario == "STANDARDCONC":
        dG0_values = copy.deepcopy(standardconc_dG0_values)
        used_concentration_values = concentration_values_free
    elif concentration_scenario == "VIVOCONC":
        dG0_values = copy.deepcopy(paperconc_dG0_values)
        used_concentration_values = concentration_values_paper

    """
    print("=PRE-ANALYSIS: SINGLE CONDITION Z SUM MINIMIZATION=")
    minimal_tcosa_reaction_set = []
    for condition_counter in range(len(conditions)):
        current_condition = [conditions[condition_counter]]
        current_condition_name = condition_names[condition_counter]

        print("Condition:", current_condition_name, current_condition)

        cobra_model = copy.deepcopy(old_cobra_model)
        cobra_model = cosa_get_model_with_nadx_scenario(
            nadx_scenario="FLEXIBLE",
            cobra_model=cobra_model,
            randoms_random_base_lists=[],
            randomfixed_random_base_lists=[],
            nad_base_ids=nad_base_ids,
            nadp_base_ids=nadp_base_ids,
        )

        print(">Get base OptMDFpathway MILP...")
        multoptmdf_base_problem = get_multoptmdfpathway_base_problem(
            cobra_model=cobra_model,
            dG0_values=dG0_values,
            metabolite_concentration_values=used_concentration_values,
            ratio_constraint_data=ratio_constraint_data,
            R=STANDARD_R,
            T=STANDARD_T,
            sub_network_ids=get_all_tcosa_reaction_ids(cobra_model),
            conditions=current_condition,
            condition_names=current_condition_name,
            add_promiscuity_constraints=True,
        )
        print(">Get model variables dictionary")
        multoptmdf_base_variables = multoptmdf_base_problem.variablesDict()

        print(">Set no promiscuity constraint")
        multoptmdf_base_problem = cosa_add_promiscuity_constraints(
            optmdfpathway_base_problem=multoptmdf_base_problem,
            optmdfpathway_base_variables=multoptmdf_base_variables,
            cobra_model=cobra_model,
            dG0_values=dG0_values,
        )

        print(">Perform TCOSA z sum minimization")
        optsubmdfpathway_result = perform_variable_minimization(
            multoptmdf_base_problem,
            "global_tcosa_z_var_sum",
        )
        print(optsubmdfpathway_result["status"])
        if optsubmdfpathway_result["status"] != "Optimal":
            print("ERROR")
            has_error = True
            break
        print("Global TCOSA Z sum:", optsubmdfpathway_result["values"]["global_tcosa_z_var_sum"])

        counter = 0
        for tcosa_reaction_id in get_all_tcosa_reaction_ids(cobra_model):
            z_var_id = "z_var_"+tcosa_reaction_id
            if optsubmdfpathway_result["values"][z_var_id] > 1e-3:
                minimal_tcosa_reaction_set.append(z_var_id.replace("z_var_", ""))
                # print("*"+z_var_id+":"+str(optsubmdfpathway_result["values"][z_var_id]))
                # counter += 1
        print(counter)
    minimal_tcosa_reaction_set = list(set(minimal_tcosa_reaction_set))
    json_write("./cosa/results_aerobic/mult_flux_samplings/minimal_tcosa_reaction_set.json", minimal_tcosa_reaction_set)
    input("END")
    """


    cobra_model = copy.deepcopy(old_cobra_model)
    cobra_model = cosa_get_model_with_nadx_scenario(
        nadx_scenario="FLEXIBLE",
        cobra_model=cobra_model,
        randoms_random_base_lists=[],
        randomfixed_random_base_lists=[],
        nad_base_ids=nad_base_ids,
        nadp_base_ids=nadp_base_ids,
    )

    print(">Get base OptMDFpathway MILP...")
    multoptmdf_base_problem = get_multoptmdfpathway_base_problem(
        cobra_model=cobra_model,
        dG0_values=dG0_values,
        metabolite_concentration_values=used_concentration_values,
        ratio_constraint_data=ratio_constraint_data,
        R=STANDARD_R,
        T=STANDARD_T,
        sub_network_ids=get_all_tcosa_reaction_ids(cobra_model),
        conditions=conditions,
        condition_names=condition_names,
        add_promiscuity_constraints=True,
    )
    print(">Get model variables dictionary")
    multoptmdf_base_variables = multoptmdf_base_problem.variablesDict()

    # print(">Set no promiscuity constraint")
    # multoptmdf_base_problem = cosa_add_promiscuity_constraints(
    #     optmdfpathway_base_problem=multoptmdf_base_problem,
    #     optmdfpathway_base_variables=multoptmdf_base_variables,
    #     cobra_model=cobra_model,
    #     dG0_values=dG0_values,
    # )

    for num_sampling in range(num_samplings):
        print("~~~")
        print("Sampling no. ", num_sampling)
        sampling_json_path = f"./cosa/results{suffix}/mult_flux_samplings4/OPTSUBMDF_{concentration_scenario}_FLEXIBLE_SAMPLING_{num_sampling}.json"

        has_error = False
        print(">Perform z sum minimization")

        optsubmdfpathway_result = perform_variable_minimization(
            multoptmdf_base_problem,
            "global_tcosa_z_var_sum",
            # "dG0_OMPDC",
        )

        print(optsubmdfpathway_result["status"])
        if optsubmdfpathway_result["status"] != "Optimal":
            print("ERROR")
            has_error = True
            break
        print("Global Z sum:", optsubmdfpathway_result["values"]["global_z_var_sum"])
        print("Global TCOSA Z sum:", optsubmdfpathway_result["values"]["global_tcosa_z_var_sum"])
        full_results = optsubmdfpathway_result

        print("Set TCOSA integer cut")
        z_sum_expression = 0.0
        z_sum = 0
        for tcosa_reaction in get_all_tcosa_reaction_ids(cobra_model):
            if not tcosa_reaction in dG0_values.keys():
                continue
            z_var_key = f"z_var_{tcosa_reaction}"
            z_value = optsubmdfpathway_result["values"][z_var_key]
            z_sum += z_value
            z_var = multoptmdf_base_variables[z_var_key]
            z_sum_expression += z_value * z_var

            if z_value > 1e-3:
                print(z_var_key, z_value)
                for condition_name in condition_names:
                    condition_z_var_key = f"z_var_{tcosa_reaction}_{condition_name}"
                    condition_z_value = optsubmdfpathway_result["values"][condition_z_var_key]
                    print(" ->", condition_z_var_key, condition_z_value)
        integer_cut_var = pulp.LpVariable(
            name=f"integer_cut_{num_sampling}",
            lowBound=-1e6,
            upBound=z_sum-1,
            cat=pulp.LpContinuous,
        )
        multoptmdf_base_problem += integer_cut_var == z_sum_expression

        if not has_error:
            json_write(
                sampling_json_path,
                full_results,
            )
        print("Global TCOSA Z sum:", optsubmdfpathway_result["values"]["global_tcosa_z_var_sum"])

        # sys.stdout = open("x.txt", "w")
        # print(multoptmdf_base_problem)


cosa_multi_flux_sampling(num_samplings=1000)
