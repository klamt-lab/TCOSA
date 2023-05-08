import cobra
import copy
import os
import numpy
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


def cosa_random_sampling(anaerobic: bool, expanded: bool, num_randoms_random: int, num_randomfixed_random: int, c_source: str="glucose", step_size: float=0.05, step_number: float=9):
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=expanded, c_source=c_source)

    suffix = cosa_get_suffix(anaerobic, expanded, c_source)
    ensure_folder_existence("./cosa")
    ensure_folder_existence(f"./cosa/results{suffix}")
    ensure_folder_existence(f"./cosa/results{suffix}/runs")

    print("Get randoms random lists")
    rng = numpy.random.default_rng(seed=11)

    def randoms_random_base_list(all_base_ids: List[str]):
        output = {
            base_id: float(rng.integers(0, 1, endpoint=True))
            for base_id in all_base_ids
        }
        for base_id in output.keys():
            if "_REV_" in (base_id+"_"):
                output[base_id] = output[base_id.replace("_REV", "_FWD")]
        return output

    randoms_random_base_lists = {
        str(i): randoms_random_base_list(all_base_ids) for i in range(num_randoms_random)
    }

    json_write(f"./cosa/results{suffix}/randoms_rand_lists.json", randoms_random_base_lists)

    print("Get randomfixed random lists")
    rng = numpy.random.default_rng(seed=22)

    def randomfixed_random_base_list(nad_base_ids: List[str], nadp_base_ids: List[str]):
        reversible_reactions = []
        for reaction in nad_base_ids+nadp_base_ids:
            if reaction.endswith("_FWD"):
                reversible_reactions.append(reaction.replace("_FWD", ""))
        nad_base_ids_new = []
        for reaction in nad_base_ids:
            if reaction.endswith("_FWD"):
                nad_base_ids_new.append(reaction.replace("_FWD", ""))
            elif reaction.endswith("_REV"):
                continue
            else:
                nad_base_ids_new.append(reaction)
        nadp_base_ids_new = []
        for reaction in nadp_base_ids:
            if reaction.endswith("_FWD"):
                nadp_base_ids_new.append(reaction.replace("_FWD", ""))
            elif reaction.endswith("_REV"):
                continue
            else:
                nadp_base_ids_new.append(reaction)

        prelist = [0.0 for _ in range(len(nad_base_ids_new))] + [1.0 for _ in range(len(nadp_base_ids_new))]
        rng.shuffle(prelist)
        all_base_ids = nad_base_ids_new + nadp_base_ids_new
        print(nadp_base_ids_new)
        output = {
            all_base_ids[i]: prelist[i]
            for i in range(len(prelist))
        }
        for reversible_reaction in reversible_reactions:
            value = output[reversible_reaction]
            output[reversible_reaction+"_FWD"] = value
            output[reversible_reaction+"_REV"] = value
            del(output[reversible_reaction])

        return output

    randomfixed_random_base_lists = {
        str(i): randomfixed_random_base_list(nad_base_ids, nadp_base_ids) for i in range(num_randomfixed_random)
    }

    json_write(f"./cosa/results{suffix}/randomfixed_rand_lists.json", randomfixed_random_base_lists)

    print(len(randomfixed_random_base_lists))
    print(len(randomfixed_random_base_lists.keys()))
    old_cobra_model = copy.deepcopy(cobra_model)
    nadx_scenarios = ["SINGLE_COFACTOR", "WILDTYPE", "FLEXIBLE"] +\
        [f"RANDOMFIXED_{i}" for i in range(len(randomfixed_random_base_lists.keys()))] +\
        [f"RANDOMS_{j}" for j in range(len(randoms_random_base_lists.keys()))]
    print(nadx_scenarios)
    original_used_growth = used_growth

    for concentration_scenario in ("STANDARDCONC", "VIVOCONC"):
        if concentration_scenario == "STANDARDCONC":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
            used_concentration_values = concentration_values_free
        elif concentration_scenario == "VIVOCONC":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper

        for nadx_scenario in nadx_scenarios:
            print("~~~")
            print(nadx_scenario)

            optmdf_json_path = f"./cosa/results{suffix}/runs/OPTMDF_{concentration_scenario}_{nadx_scenario}.json"
            optsubmdf_json_path = f"./cosa/results{suffix}/runs/OPTSUBMDF_{concentration_scenario}_{nadx_scenario}.json"

            if os.path.exists(optmdf_json_path+".zip") and os.path.exists(optsubmdf_json_path+".zip"):
                print("Already calculated!")
                continue

            cobra_model = copy.deepcopy(old_cobra_model)
            cobra_model = cosa_get_model_with_nadx_scenario(
                nadx_scenario=nadx_scenario,
                cobra_model=cobra_model,
                randoms_random_base_lists=randoms_random_base_lists,
                randomfixed_random_base_lists=randomfixed_random_base_lists,
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
            is_prelast_round = False
            is_last_round = False
            while True: # (used_growth > 0.05) and (used_growth > original_used_growth - step_size * step_number):
                if is_prelast_round:
                    used_growth = 0.03
                    is_last_round = True
                if (used_growth <= 0.05) and (not is_last_round):
                    used_growth = 0.05
                    is_prelast_round = True
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
                    has_error = True
                    break
                full_optmdf_results[rounded_used_growth] = optmdfpathway_result

                print(">SubMDF calculations")
                optmdfpathway_base_variables["var_B"].bounds(MIN_OPTMDF, 1e6)
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

    create_cosa_tables(data_path=f"cosa/results{suffix}/runs", output_path=f"cosa/results{suffix}")
    create_cosa_figures(data_path=f"./cosa/results{suffix}/", figures_path=f"./cosa/results{suffix}/figures/", anaerobic=anaerobic)


cosa_random_sampling(anaerobic=False, expanded=False, num_randoms_random=50, num_randomfixed_random=50, c_source="acetate")
# cosa_random_sampling(anaerobic=True, expanded=False, num_randoms_random=1, num_randomfixed_random=1, c_source="acetate")
"""
cosa_random_sampling(anaerobic=False, expanded=False, num_randoms_random=500, num_randomfixed_random=500)
cosa_random_sampling(anaerobic=True, expanded=False, num_randoms_random=500, num_randomfixed_random=500)
cosa_random_sampling(anaerobic=False, expanded=True, num_randoms_random=1, num_randomfixed_random=1)
cosa_random_sampling(anaerobic=True, expanded=True, num_randoms_random=1, num_randomfixed_random=1)
create_total_cosa_figure()
"""
