import cobra
import copy
import os
import numpy
import pulp
from typing import List
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from cosa_get_suffix import cosa_get_suffix
from fba import get_fba_base_problem, perform_fba_flux_maximization
from helper import ensure_folder_existence, json_write, json_zip_write
from optmdfpathway import STANDARD_R, STANDARD_T, get_optmdfpathway_base_problem, get_thermodynamic_bottlenecks
from optimization import perform_variable_maximization
from cosa_create_table import create_cosa_tables
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from cosa_random_sampling_figures import create_cosa_figures, create_total_cosa_figure
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from cosa_add_promiscuity_constraints import cosa_add_promiscuity_constraints


def cosa_flux_sampling(anaerobic: bool, growth_key: str, mu: float, optsubmdf: float, num_samplings: int):
    expanded = False

    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=expanded)

    suffix = cosa_get_suffix(anaerobic, expanded)

    ensure_folder_existence("./cosa")
    ensure_folder_existence(f"./cosa/results{suffix}")
    ensure_folder_existence(f"./cosa/results{suffix}/flux_samplings")

    old_cobra_model = copy.deepcopy(cobra_model)
    for concentration_scenario in ("STANDARDCONC",): #"VIVOCONC"):
        if concentration_scenario == "STANDARDCONC":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
            used_concentration_values = concentration_values_free
        elif concentration_scenario == "VIVOCONC":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper


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
            optmdfpathway_base_problem=optmdfpathway_base_variables,
            optmdfpathway_base_variables=optmdfpathway_base_variables,
            cobra_model=cobra_model,
            dG0_values=dG0_values,
        )

        for num_sampling in range(num_samplings):
            print("~~~")
            print(num_sampling)
            sampling_json_path = f"./cosa/results{suffix}/flux_samplings/OPTSUBMDF_{concentration_scenario}_FLEXIBLE_SAMPLING_{num_sampling}.json"

            print(">Perform test FBA...")
            biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"
            print(f" Selected biomass reaction: {biomass_reaction_id}")
            fba_base_problem = get_fba_base_problem(
                cobra_model=cobra_model, extra_constraints=[])
            fba_result = perform_fba_flux_maximization(
                base_problem=fba_base_problem, reaction_id=biomass_reaction_id)
            precise_max_growth = fba_result["values"][biomass_reaction_id]

            print(f" Precise max growth is {precise_max_growth}")

            has_error = False
            full_optsubmdf_results = {}
            print("Set growth to", mu)
            optmdfpathway_base_variables[biomass_reaction_id].bounds(
                mu,
                1e12
            )
            print(">SubMDF calculations")
            optmdfpathway_base_variables["var_B"].bounds(MIN_OPTMDF, 1e6)
            optmdfpathway_base_variables["var_B2"].bounds(optsubmdf, 1e6)
            optsubmdfpathway_result = perform_variable_maximization(
                optmdfpathway_base_problem,
                "var_B2"
            )
            print(optsubmdfpathway_result["status"])
            if optsubmdfpathway_result["status"] != "Optimal":
                has_error = True
                break
            print("SubMDF:", optsubmdfpathway_result["values"]["var_B2"])
            full_optsubmdf_results = optsubmdfpathway_result

            print("Set integer cut")
            z_sum_expression = 0.0
            z_sum = 0
            for tcosa_reaction in get_all_tcosa_reaction_ids(cobra_model):
                if not tcosa_reaction in dG0_values.keys():
                    continue
                z_var_key = f"z_var_{tcosa_reaction}"
                z_value = optsubmdfpathway_result["values"][z_var_key]
                z_sum += z_value
                z_var = optmdfpathway_base_variables[z_var_key]
                z_sum_expression += z_value * z_var
            integer_cut_var = pulp.LpVariable(
                name=f"integer_cut_{num_sampling}",
                lowBound=-1e6,
                upBound=z_sum-1,
                cat=pulp.LpContinuous,
            )
            optmdfpathway_base_problem += integer_cut_var == z_sum_expression

            if not has_error:
                json_write(
                    sampling_json_path,
                    full_optsubmdf_results,
                )


cosa_flux_sampling(anaerobic=False, growth_key="0,568", mu=0.568, optsubmdf=29.0137, num_samplings=100)
