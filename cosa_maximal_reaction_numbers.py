from tkinter.messagebox import NO
import matplotlib.pyplot as plt
import cobra
import copy
import pulp
import math
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from cosa_get_suffix import cosa_get_suffix
from helper import json_load, json_write, json_zip_load
from typing import List
from optmdfpathway import (
    STANDARD_R, STANDARD_T, get_optmdfpathway_base_problem,
    add_differential_reactions_constraints, get_z_variable_status,
    get_thermodynamic_bottlenecks
)
from optimization import perform_variable_minimization, perform_variable_maximization
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from typing import Dict
from helper import ensure_folder_existence


def cosa_maximal_reaction_numbers(anaerobic: bool, expanded: bool, growth_epsilon: float = 0.01) -> None:
    suffix = cosa_get_suffix(anaerobic, expanded)
    figures_path = f"./cosa/results{suffix}/figures/"
    ensure_folder_existence(figures_path)
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=expanded)

    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    original_cobra_model = copy.deepcopy(cobra_model)
    for concentrations in ("STANDARDCONC", "VIVOCONC"):
        print(f"CONCENTRATIONS: {concentrations}")
        print(f"=CONCENTRATION RANGES: {concentrations}=")
        if concentrations == "STANDARDCONC":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
            used_concentration_values = concentration_values_free
        elif concentrations == "VIVOCONC":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper

        for target in ("OPTMDF", "OPTSUBMDF"):
            print(f"TARGET: {target}")
            print(f"===OPTIMIZATION TARGET: {target}===")

            for distribution in ("WILDTYPE", "FLEXIBLE", "SINGLE_COFACTOR"):
                print(f"DISTRIBUTION: {distribution}")
                print(f"cosa/results{suffix}/runs/{target}_{concentrations}_{distribution}.json")
                jsondata_invivo = json_zip_load(f"cosa/results{suffix}/runs/{target}_{concentrations}_{distribution}.json")
                growth_rates = jsondata_invivo.keys()
                report = ""

                cobra_model = copy.deepcopy(original_cobra_model)
                cobra_model = cosa_get_model_with_nadx_scenario(
                    nadx_scenario=distribution,
                    cobra_model=cobra_model,
                )

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
                optmdfpathway_base_variables: Dict[str, pulp.LpVariable] = optmdfpathway_base_problem.variablesDict()
                z_sum = 0.0
                for reaction in get_all_tcosa_reaction_ids(cobra_model):
                    z_sum += optmdfpathway_base_variables[f"z_var_{reaction}"]
                z_sum_var = pulp.LpVariable(
                    name=f"z_sum_var",
                    lowBound=-float("inf"),
                    upBound=float("inf"),
                    cat=pulp.LpContinuous,
                )
                optmdfpathway_base_problem += z_sum_var == z_sum


                for growth_rate in growth_rates:
                    growth_rate_float = float(growth_rate.replace(",", "."))
                    optmdfpathway_base_variables[biomass_reaction_id].bounds(
                        growth_rate_float-growth_epsilon,
                        1e12,
                    )

                    if target == "OPTMDF":
                        min_target = jsondata_invivo[growth_rate]["values"]["var_B"] - 0.001
                        optmdfpathway_base_variables["var_B"].bounds(
                            min_target,
                            1e12,
                        )
                    elif target == "OPTSUBMDF":
                        min_target = jsondata_invivo[growth_rate]["values"]["var_B2"] - 0.001
                        optmdfpathway_base_variables["var_B"].bounds(
                            MIN_OPTMDF,
                            1e12,
                        )
                        optmdfpathway_base_variables["var_B2"].bounds(
                            min_target,
                            1e12,
                        )

                    print(f" @ µ [1/h] of {growth_rate_float} and min {target} of {min_target} kJ/mol")
                    report += f" @ µ [1/h] of {growth_rate_float} and min {target} of {min_target} kJ/mol\n"

                    print("RUN")
                    minimization_result = perform_variable_maximization(
                        optmdfpathway_base_problem,
                        "z_sum_var",
                    )
                    all_tcosas = get_all_tcosa_reaction_ids(cobra_model)
                    active_reactions = []
                    for var_name in minimization_result["values"].keys():
                        if not var_name.startswith("z_var_"):
                            continue
                        reac_id = var_name.replace("z_var_", "")
                        if reac_id not in all_tcosas:
                            continue
                        z_value = minimization_result["values"][var_name]
                        reaction_string = cobra_model.reactions.get_by_id(reac_id).reaction
                        if z_value > 0.1:
                            report_line = f"* {reac_id} | {round(dG0_values[reac_id]['dG0'], 3)} kJ/mol | {reaction_string}"
                            active_reactions.append(reac_id)
                            # print(report_line)
                            report += report_line + "\n"
                    report += f"Reached MDF: {minimization_result['values']['var_B']}\n"
                    report += f"Reached SubMDF: {minimization_result['values']['var_B2']}\n"
                    report += f"Minimal sum of active reactions: {len(active_reactions)}\n"
                    # print(f"Reached MDF: {minimization_result['values']['var_B']} kJ/mol")
                    # print(f"Reached SubMDF: {minimization_result['values']['var_B2']} kJ/mol")

                    aerobic = "anaerobic" if anaerobic else "aerobic"
                with open(f"./cosa/max_reacs_{aerobic}_{distribution}_{target}_{concentrations}.txt", "w", encoding="utf-8") as f:
                    f.write(report)


cosa_maximal_reaction_numbers(anaerobic=False, expanded=False)
cosa_maximal_reaction_numbers(anaerobic=True, expanded=False)
