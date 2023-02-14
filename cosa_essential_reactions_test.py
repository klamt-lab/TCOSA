from select import select
import cobra
import copy
import os
import pulp
from typing import List
from fba import get_fba_base_problem, perform_fba_flux_maximization
from helper import ensure_folder_existence, json_write, json_zip_write
from optmdfpathway import STANDARD_R, STANDARD_T, get_optmdfpathway_base_problem, perform_concentration_variability_analysis
from optimization import perform_variable_maximization, perform_variable_minimization
from cosa_create_table import create_cosa_tables
from cosa_load_model_data import (
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,
    standardconc_dG0_values, paperconc_dG0_values, MIN_OPTMDF, nad_reactions, nadp_reactions,
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,
    ratio_constraint_data, nad_base_ids, nadp_base_ids,
)
from cosa_random_sampling_figures import create_cosa_figures

old_cobra_model = copy.deepcopy(cobra_model)
nadx_scenarios = ["WILDTYPE", "SINGLE_COFACTOR", "FLEXIBLE"]
print(nadx_scenarios)

nadx_z_vars = [f"z_var_{x}" for x in nad_reactions+nadp_reactions]
# reaction_ids = [x.id for x in cobra_model.reactions]
# nadx_z_vars = [f"z_var_{x}" for x in reaction_ids]
print("Full number NADX reactions:", len(nadx_z_vars)/2)
for concentration_scenario in ("STANDARDCONC", "VIVOCONC"):
    print(f"====CONCENTRATION SCENARIO: {concentration_scenario}====")
    if concentration_scenario == "STANDARDCONC":
        dG0_values = copy.deepcopy(standardconc_dG0_values)
        used_concentration_values = concentration_values_free
    elif concentration_scenario == "VIVOCONC":
        dG0_values = copy.deepcopy(paperconc_dG0_values)
        used_concentration_values = concentration_values_paper

    for nadx_scenario in nadx_scenarios:
        print("~~~")
        print(nadx_scenario)

        cobra_model = copy.deepcopy(old_cobra_model)
        if nadx_scenario == "WILDTYPE":
            for reaction in cobra_model.reactions:
                reaction: cobra.Reaction = reaction
                key = reaction.id
                if (key in nadp_reactions) and (key.endswith("_NADX")):
                    reaction.lower_bound = 0.0
                    reaction.upper_bound = 0.0
                elif (key in nad_reactions) and (key.endswith("_NADY")):
                    reaction.lower_bound = 0.0
                    reaction.upper_bound = 0.0
        elif nadx_scenario == "SINGLE_COFACTOR":
            for reaction in cobra_model.reactions:
                reaction: cobra.Reaction = reaction
                met_ids = [x.id for x in reaction.metabolites.keys()]
                if ("nady_c" in met_ids) or ("nadyh_c" in met_ids):
                    reaction.lower_bound = 0.0
                    reaction.upper_bound = 0.0

        print(">Get base OptMDFpathway MILP...")
        optmdfpathway_base_problem = get_optmdfpathway_base_problem(
            cobra_model=cobra_model,
            dG0_values=dG0_values,
            metabolite_concentration_values=used_concentration_values,
            ratio_constraint_data=ratio_constraint_data,
            R=STANDARD_R,
            T=STANDARD_T,
            extra_constraints=[],
            sub_network_ids=nad_reactions+nadp_reactions,
        )

        print(">Get model variables dictionary")
        optmdfpathway_base_variables = optmdfpathway_base_problem.variablesDict()

        print(">Add NADX z var sum")
        nadx_z_var_sum_var = pulp.LpVariable(
            name="nadx_z_var_sum",
            cat=pulp.LpContinuous,
        )
        nadx_z_var_sum_constraint = 0.0
        added_z_vars = []
        for reac_id in nadx_z_vars:
            if reac_id in optmdfpathway_base_variables.keys():
                added_z_vars.append(reac_id)
                nadx_z_var_sum_constraint += optmdfpathway_base_variables[reac_id]
        optmdfpathway_base_problem += nadx_z_var_sum_constraint == nadx_z_var_sum_var

        print(">Get model variables dictionary again")
        optmdfpathway_base_variables = optmdfpathway_base_problem.variablesDict()

        print(">Perform test FBA...")
        biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"
        print(f" Selected biomass reaction: {biomass_reaction_id}")
        fba_base_problem = get_fba_base_problem(
            cobra_model=cobra_model,
            extra_constraints=[]
        )
        fba_result = perform_fba_flux_maximization(
            base_problem=fba_base_problem,
            reaction_id=biomass_reaction_id
        )
        precise_max_growth = fba_result["values"][biomass_reaction_id]
        used_max_growth = 0.876  # float(str(precise_max_growth)[:5])

        print(f" Precise max growth is {precise_max_growth}")
        print(f" Used maximal rounded and floored max growth is {used_max_growth}")

        used_growth = used_max_growth
        has_error = False

        print("Set growth to", used_growth)
        optmdfpathway_base_variables[biomass_reaction_id].bounds(
            used_growth,
            1e12
        )

        print(">OPTMDF calculations")
        optmdfpathway_result = perform_variable_maximization(
            optmdfpathway_base_problem,
            "var_B"
        )
        print(optmdfpathway_result["status"])
        print("OPTMDF:", optmdfpathway_result["values"]["var_B"])
        optmdfpathway_base_variables["var_B"].bounds(optmdfpathway_result["values"]["var_B"], 1e6)
        essential_reactions_result = perform_variable_minimization(
            optmdfpathway_base_problem,
            "nadx_z_var_sum"
        )
        print(">Minimal essential reactions:", essential_reactions_result["objective_value"], essential_reactions_result["values"]["nadx_z_var_sum"])

        for added_z_var in added_z_vars:
            if essential_reactions_result["values"][added_z_var] > 1e-6:
                print("*", added_z_var, essential_reactions_result["values"][added_z_var])

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        print(">SubMDF calculations")
        optmdfpathway_base_variables["var_B"].bounds(MIN_OPTMDF, 1e6)
        optsubmdfpathway_result = perform_variable_maximization(
            optmdfpathway_base_problem,
            "var_B2"
        )
        print(optsubmdfpathway_result["status"])
        print("SubMDF:", optsubmdfpathway_result["values"]["var_B2"])
        optmdfpathway_base_variables["var_B2"].bounds(optsubmdfpathway_result["values"]["var_B2"], 1e6)
        essential_reactions_result = perform_variable_minimization(
            optmdfpathway_base_problem,
            "nadx_z_var_sum"
        )
        print(">Minimal essential reactions:", essential_reactions_result["objective_value"], essential_reactions_result["values"]["nadx_z_var_sum"])

        for added_z_var in added_z_vars:
            if essential_reactions_result["values"][added_z_var] > 1e-6:
                print("*", added_z_var, essential_reactions_result["values"][added_z_var])

        optmdfpathway_base_variables["var_B2"].bounds(-1e6, 1e6)
        print("=========================")
