"""Performs the dGf sampling, as done for the TCOSA publication."""

# IMPORTS #
# External
import cobra
import copy
import os
import numpy
from typing import List
# Internal
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from cosa_get_suffix import cosa_get_suffix
from fba import get_fba_base_problem, perform_fba_flux_maximization
from helper import ensure_folder_existence, json_write, json_zip_write
from optmdfpathway import STANDARD_R, STANDARD_T, get_optmdfpathway_base_problem
from optimization import perform_variable_maximization
from cosa_create_table import create_cosa_dG0_sampling_tables
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from cosa_random_sampling_figures import create_cosa_dG0_sampling_figures, create_total_dG0_sampling_figure
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from cosa_add_promiscuity_constraints import cosa_add_promiscuity_constraints
from typing import Dict


# PUBLIC FUNCTIONS #
def cosa_dG0_sampling(anaerobic: bool, expanded: bool, num_samplings: int, step_size: float=0.05, change_range: float=25.0):
    """Performs the TCOSA dGf sampling as given.

    Args:
        anaerobic (bool): Is it anaerobic?
        expanded (bool): Is is a 2-cofactor model (False) or not (True)?
        num_samplings (int): Number of dG0 samplings for each 'base' specificity'.
        step_size (float, optional): The µ step size. Defaults to 0.05.
        change_range (float, optional): The maximal dGf change. Defaults to 25.0.

    Returns:
        _type_: _description_
    """
    suffix = cosa_get_suffix(anaerobic, expanded)
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=expanded)


    ensure_folder_existence("./cosa")
    ensure_folder_existence(f"./cosa/results{suffix}")
    ensure_folder_existence(f"./cosa/results{suffix}/dG0_sampling_range{change_range}")
    ensure_folder_existence(f"./cosa/results{suffix}/dG0_sampling_range{change_range}/runs")

    print("Get randoms random lists")
    rng = numpy.random.default_rng(seed=22)

    def get_random_dGf_change_dict() -> Dict[str, float]:
        random_changes = {
            x.id: float(rng.uniform(-change_range, +change_range))
            for x in cobra_model.metabolites
        }
        return random_changes

    excluded_metabolites = ["h_c", "h_p", "h_e", "h2o_c", "h2o_p", "h2o_e", "o2_c", "o2_p", "o2_e", "co2_c", "co2_p", "co2_e"]

    random_standardconc_dG0_values = {}
    for i in range(num_samplings):
        random_dGf_change_dict = get_random_dGf_change_dict()
        random_standardconc_dG0_values[i] = copy.deepcopy(standardconc_dG0_values)
        for reaction in cobra_model.reactions:
            if reaction.id not in random_standardconc_dG0_values[i].keys():
                continue
            for key, value in reaction.metabolites.items():
                if key.id in excluded_metabolites:
                    continue
                stoichiometry = -value
                random_dGf_change = copy.deepcopy(random_dGf_change_dict[key.id]) * stoichiometry
                random_standardconc_dG0_values[i][reaction.id]["dG0"] = copy.deepcopy(random_dGf_change) + copy.deepcopy(random_standardconc_dG0_values[i][reaction.id]["dG0"])

    random_paperconc_dG0_values = {}
    for i in range(num_samplings):
        random_dGf_change_dict = get_random_dGf_change_dict()
        random_paperconc_dG0_values[i] = copy.deepcopy(paperconc_dG0_values)
        for reaction in cobra_model.reactions:
            if reaction.id not in random_paperconc_dG0_values[i].keys():
                continue
            for key, value in reaction.metabolites.items():
                if key.id in excluded_metabolites:
                    continue
                stoichiometry = -value
                random_dGf_change = random_dGf_change_dict[key.id] * stoichiometry
                random_paperconc_dG0_values[i][reaction.id]["dG0"] += random_dGf_change

    json_zip_write(f"./cosa/results{suffix}/dG0_sampling_range{change_range}/random_standard_dG0_changed_list.json", random_standardconc_dG0_values)
    json_zip_write(f"./cosa/results{suffix}/dG0_sampling_range{change_range}/random_dG0_changes_list.json", random_dGf_change_dict)
    json_zip_write(f"./cosa/results{suffix}/dG0_sampling_range{change_range}/random_paper_dG0_changed_list.json", random_paperconc_dG0_values)

    old_cobra_model = copy.deepcopy(cobra_model)
    nadx_scenarios = ["SINGLE_COFACTOR", "WILDTYPE", "FLEXIBLE"] +\
        [f"RANDOM_SINGLE_COFACTOR_{i}" for i in range(num_samplings)] +\
        [f"RANDOM_WILDTYPE_{j}" for j in range(num_samplings)] +\
        [f"RANDOM_FLEXIBLE_{k}" for k in range(num_samplings)]
    print(nadx_scenarios)
    original_used_growth = used_growth

    for concentration_scenario in ("STANDARDCONC",): # "VIVOCONC",
        if concentration_scenario == "STANDARDCONC":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
            used_concentration_values = concentration_values_free
            random_dG0_dicts_list = copy.deepcopy(random_standardconc_dG0_values)
        elif concentration_scenario == "VIVOCONC":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper
            random_dG0_dicts_list = copy.deepcopy(random_paperconc_dG0_values)
        old_dG0_values = copy.deepcopy(dG0_values)

        for nadx_scenario in nadx_scenarios:
            print("~~~")
            print(nadx_scenario)

            optmdf_json_path = f"./cosa/results{suffix}/dG0_sampling_range{change_range}/runs/OPTMDF_{concentration_scenario}_{nadx_scenario}.json"
            optsubmdf_json_path = f"./cosa/results{suffix}/dG0_sampling_range{change_range}/runs/OPTSUBMDF_{concentration_scenario}_{nadx_scenario}.json"

            if os.path.exists(optmdf_json_path+".zip") and os.path.exists(optsubmdf_json_path+".zip"):
                print("Already calculated!")
                continue

            if "SINGLE_COFACTOR" in nadx_scenario:
                base_nadx_scenario = "SINGLE_COFACTOR"
            elif "WILDTYPE" in nadx_scenario:
                base_nadx_scenario = "WILDTYPE"
            elif "FLEXIBLE" in nadx_scenario:
                base_nadx_scenario = "FLEXIBLE"
            cobra_model = copy.deepcopy(old_cobra_model)
            cobra_model = cosa_get_model_with_nadx_scenario(
                nadx_scenario=base_nadx_scenario,
                cobra_model=cobra_model,
                randoms_random_base_lists=[],
                randomfixed_random_base_lists=[],
            )
            dG0_values = copy.deepcopy(old_dG0_values)
            if nadx_scenario.startswith("RANDOM_"):
                random_id = int(nadx_scenario.split("_")[-1])
                dG0_values = random_dG0_dicts_list[random_id]

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

            print(">Perform test FBA...")
            biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"
            print(f" Selected biomass reaction: {biomass_reaction_id}")
            fba_base_problem = get_fba_base_problem(
                cobra_model=cobra_model, extra_constraints=[])
            fba_result = perform_fba_flux_maximization(
                base_problem=fba_base_problem, reaction_id=biomass_reaction_id)
            precise_max_growth = fba_result["values"][biomass_reaction_id]
            used_max_growth = original_used_growth

            print(f" Precise max growth is {precise_max_growth}")
            print(f" Used maximal rounded and floored max growth is {used_max_growth}")

            used_growth = used_max_growth
            has_error = False
            full_optmdf_results = {}
            full_optsubmdf_results = {}
            is_last_round = False
            while True:
                if used_growth <= 0.05:
                    used_growth = 0.05
                    is_last_round = True

                rounded_used_growth = str(round(used_growth, 3)).replace(".", ",")

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
                if optmdfpathway_result["status"] != "Optimal":
                    has_error = True
                    break
                print("var_B", optmdfpathway_result["values"]["var_B"], "kJ/mol")
                print("Growth", optmdfpathway_result["values"][biomass_reaction_id], "kJ/mol")

                if optmdfpathway_result["values"]["var_B"] < MIN_OPTMDF:
                    # has_error = True
                    # break
                    used_minoptmdf = optmdfpathway_result["values"]["var_B"]
                else:
                    used_minoptmdf = MIN_OPTMDF
                full_optmdf_results[rounded_used_growth] = optmdfpathway_result

                print(">SubMDF calculations")
                optmdfpathway_base_variables["var_B"].bounds(used_minoptmdf, 1e6)
                optsubmdfpathway_result = perform_variable_maximization(
                    optmdfpathway_base_problem,
                    "var_B2"
                )
                print(optsubmdfpathway_result["status"])
                if optsubmdfpathway_result["status"] != "Optimal":
                    has_error = True
                    break
                print("SubMDF:", optsubmdfpathway_result["values"]["var_B2"])
                full_optsubmdf_results[rounded_used_growth] = optsubmdfpathway_result

                used_growth -= step_size
                if is_last_round:
                    break
            if not has_error:
                json_zip_write(
                    optmdf_json_path,
                    full_optmdf_results,
                )
                json_zip_write(
                    optsubmdf_json_path,
                    full_optsubmdf_results,
                )

    create_cosa_dG0_sampling_tables(data_path=f"cosa/results{suffix}/dG0_sampling_range{change_range}/runs", output_path=f"cosa/results{suffix}/dG0_sampling_range{change_range}")
    create_cosa_dG0_sampling_figures(data_path=f"./cosa/results{suffix}/dG0_sampling_range{change_range}/", figures_path=f"./cosa/results{suffix}/dG0_sampling_range{change_range}/figures/", anaerobic=anaerobic, num_samplings=num_samplings)
