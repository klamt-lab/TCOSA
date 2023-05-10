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
)
from optimization import perform_variable_minimization, perform_variable_maximization
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from typing import Dict
from helper import ensure_folder_existence


def cosa_test_switch(anaerobic: bool, expanded: bool, growth_epsilon: float = 0.01, c_source: str="glucose") -> None:
    suffix = cosa_get_suffix(anaerobic, expanded, c_source)
    figures_path = f"./cosa/results{suffix}/figures/"
    ensure_folder_existence(figures_path)
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=expanded, c_source=c_source)

    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    report = ""
    original_cobra_model = copy.deepcopy(cobra_model)
    for concentrations in ("STANDARDCONC", "VIVOCONC"):
        print(f"=CONCENTRATION RANGES: {concentrations}=")
        report += f"=CONCENTRATION RANGES: {concentrations}=\n"
        if concentrations == "STANDARDCONC":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
            used_concentration_values = concentration_values_free
        elif concentrations == "VIVOCONC":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper

        for target in ("OPTMDF", "OPTSUBMDF"):
            print(f"===OPTIMIZATION TARGET: {target}===")
            report += f"===OPTIMIZATION TARGET: {target}===\n"

            jsondata_invivo = json_zip_load(f"cosa/results{suffix}/runs/{target}_{concentrations}_FLEXIBLE.json")
            growth_rates = jsondata_invivo.keys()

            for growth_rate in growth_rates:
                cobra_model = copy.deepcopy(original_cobra_model)
                cobra_model = cosa_get_model_with_nadx_scenario(
                    nadx_scenario="FLEXIBLE",
                    cobra_model=cobra_model,
                )

                for reaction in cobra_model.reactions:
                    if (abs(jsondata_invivo[growth_rate]["values"][reaction.id]) < 1e-9):
                        continue
                    switched_id = copy.deepcopy(reaction.id)
                    if "_ORIGINAL_NAD_" in switched_id:
                        switched_id = switched_id.replace("_ORIGINAL_NAD_", "_VARIANT_NADP_")
                    elif "_ORIGINAL_NADP_" in switched_id:
                        switched_id = switched_id.replace("_ORIGINAL_NADP_", "_VARIANT_NAD_")
                    elif "_VARIANT_NAD_" in switched_id:
                        switched_id = switched_id.replace("_VARIANT_NAD_", "_ORIGINAL_NADP_")
                    elif "_VARIANT_NADP_" in switched_id:
                        switched_id = switched_id.replace("_VARIANT_NADP_", "_ORIGINAL_NAD_")
                    else:
                        continue
                    cobra_model.reactions.get_by_id(switched_id).lower_bound = cobra_model.reactions.get_by_id(reaction.id).lower_bound
                    cobra_model.reactions.get_by_id(switched_id).upper_bound = cobra_model.reactions.get_by_id(reaction.id).upper_bound

                    cobra_model.reactions.get_by_id(reaction.id).upper_bound = 0.0
                    cobra_model.reactions.get_by_id(reaction.id).upper_bound = 0.0

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
                optmdfpathway_base_variables: Dict[str, pulp.LpVariable] = optmdfpathway_base_problem.variablesDict()

                growth_rate_float = float(growth_rate.replace(",", "."))
                optmdfpathway_base_variables[biomass_reaction_id].bounds(
                    growth_rate_float-growth_epsilon,
                    1e12,
                )

                if target == "OPTMDF":
                    # min_target = jsondata_invivo[growth_rate]["values"]["var_B"]
                    # optmdfpathway_base_variables["var_B"].bounds(
                    #     min_target,
                    #     1e12,
                    # )
                    target_var = "var_B"
                elif target == "OPTSUBMDF":
                    min_target = jsondata_invivo[growth_rate]["values"]["var_B2"]
                    optmdfpathway_base_variables["var_B"].bounds(
                        MIN_OPTMDF,
                        1e12,
                    )
                    optmdfpathway_base_variables["var_B2"].bounds(
                        min_target,
                        1e12,
                    )
                    target_var = "var_B2"

                print(f" @ {target} @ µ [1/h] of {growth_rate_float} kJ/mol")
                report += f" @ {target} @ µ [1/h] of {growth_rate_float}  kJ/mol\n"

                maximization_result = perform_variable_maximization(
                    optmdfpathway_base_problem,
                    target_var,
                )

                old_target = jsondata_invivo[growth_rate]["values"][target_var]
                new_target = maximization_result["values"][target_var]

                print("COMPARISON: ", old_target, new_target)

cosa_test_switch(anaerobic=False, expanded=False, c_source="glucose")

# cosa_ratio_ratio_test(anaerobic=False, expanded=False)
# cosa_ratio_ratio_test(anaerobic=True, expanded=False)
# cosa_create_full_ratio_ratio_test_figure_one_panel()
# cosa_create_full_ratio_ratio_test_figure_two_panels()
# cosa_create_full_ratio_ratio_test_figure_four_panels()
