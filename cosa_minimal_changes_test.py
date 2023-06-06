"""Contains the function for the minimal number of reactions from wild-type to optimum analysis."""

# IMPORTS #
# Internal
import copy
import cobra
import pulp
from typing import Dict
import matplotlib.pyplot as plt
# External
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from cosa_get_suffix import cosa_get_suffix
from helper import ensure_folder_existence, json_load, json_write, json_zip_load
from optmdfpathway import (
    STANDARD_R, STANDARD_T, get_optmdfpathway_base_problem,
    add_differential_reactions_constraints, get_z_variable_status,
)
from optimization import perform_variable_minimization, perform_variable_maximization
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)


# PUBLIC FUNCTION #
def cosa_minimal_changes_test(anaerobic: bool, disallowed_changed_reaction: str="", c_source: str="glucose") -> None:
    """Performs the minimal number of reactions from wild-type to optimum analysis.

    Args:
        anaerobic (bool): Is it anaerobic (True)?
        disallowed_changed_reaction (str, optional): A reaction which is not allowed to be cofactor switched (only if given). Defaults to "".
        c_source (str, optional): Either 'glucose' or 'acetate'. Defaults to "glucose".
    """
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=False, c_source=c_source)

    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    suffix = cosa_get_suffix(anaerobic, False, c_source)

    figures_path = f"./cosa/results{suffix}/figures/"
    ensure_folder_existence(figures_path)

    if c_source == "glucose":
        concentrations = ("PAPERCONCS", "STANDARDCONCS")
    else:
        concentrations = ("STANDARDCONCS",)

    report = ""
    for concentrations in concentrations:
        print(f"=CONCENTRATION RANGES: {concentrations}=")
        report += f"=CONCENTRATION RANGES: {concentrations}=\n"

        if concentrations == "STANDARDCONCS":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
            used_concentration_values = concentration_values_free
        elif concentrations == "PAPERCONCS":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper

        num_changes_dict = {}
        for target in ("OptSubMDF", "OptMDF"):
            full_scenario_key = f"{target} target, {'aerobic' if (not anaerobic) else 'anaerobic'} conditions"
            num_changes_dict[full_scenario_key] = []

            print(f"===OPTIMIZATION TARGET: {target}===")
            report += f"===OPTIMIZATION TARGET: {target}===\n"
            if concentrations == "PAPERCONCS":
                jsondata_flexible = json_zip_load(f"cosa/results{suffix}/runs/{target.upper()}_VIVOCONC_FLEXIBLE.json")
            else:
                jsondata_flexible = json_zip_load(f"cosa/results{suffix}/runs/{target.upper()}_STANDARDCONC_FLEXIBLE.json")
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

            error_sum = 0.0
            for reaction in cobra_model.reactions:
                reaction: cobra.Reaction = reaction

                if (not reaction.id.endswith("_TCOSA")):
                    continue

                if f"dG0_{reaction.id}" not in optmdfpathway_base_variables.keys():
                    continue

                if "_ORIGINAL_" in reaction.id:
                    continue

                flux_variable = optmdfpathway_base_variables[reaction.id]
                flux_variable.bounds(
                    reaction.lower_bound,
                    1e12,
                )

                error_z_variable = pulp.LpVariable(
                    name=f"z_{reaction.id}_error",
                    cat=pulp.LpBinary,
                )

                if reaction.id.endswith("_ORIGINAL_NAD_TCOSA"):
                    other_id = reaction.id.replace("_ORIGINAL_NAD_TCOSA", "_VARIANT_NADP_TCOSA")
                elif reaction.id.endswith("_VARIANT_NAD_TCOSA"):
                    other_id = reaction.id.replace("_VARIANT_NAD_TCOSA", "_ORIGINAL_NADP_TCOSA")
                elif reaction.id.endswith("_ORIGINAL_NADP_TCOSA"):
                    other_id = reaction.id.replace("_ORIGINAL_NADP_TCOSA", "_VARIANT_NAD_TCOSA")
                elif reaction.id.endswith("_VARIANT_NADP_TCOSA"):
                    other_id = reaction.id.replace("_VARIANT_NADP_TCOSA", "_ORIGINAL_NAD_TCOSA")
                else:
                    continue

                current_z_variable = optmdfpathway_base_variables[f"z_var_{reaction.id}"]
                other_z_variable = optmdfpathway_base_variables[f"z_var_{other_id}"]

                optmdfpathway_base_problem += flux_variable <= error_z_variable * 1_000
                optmdfpathway_base_problem += current_z_variable + other_z_variable <= 1

                if reaction.id == disallowed_changed_reaction:
                    optmdfpathway_base_problem += error_z_variable <= 1e-5

                error_sum += error_z_variable

            error_sum_var = pulp.LpVariable(
                name=f"reaction_error_sum",
                cat=pulp.LpContinuous,
            )
            optmdfpathway_base_problem += error_sum_var == error_sum

            growth_rates = jsondata_flexible.keys()
            for growth_rate in jsondata_flexible.keys():
                growth_rate_float = float(growth_rate.replace(",", "."))
                optmdfpathway_base_variables[biomass_reaction_id].bounds(
                    growth_rate_float,
                    1e12,
                )

                if target == "OptSubMDF":
                    min_target = jsondata_flexible[growth_rate]["values"]["var_B2"]
                    optmdfpathway_base_variables["var_B"].bounds(
                        MIN_OPTMDF,
                        1e12,
                    )
                    optmdfpathway_base_variables["var_B2"].bounds(
                        min_target,
                        1e12,
                    )
                else:
                    min_target = jsondata_flexible[growth_rate]["values"]["var_B"]
                    optmdfpathway_base_variables["var_B"].bounds(
                        min_target,
                        1e12,
                    )
                print(f"@ µ of {growth_rate_float} 1/h and minimal {target} of {min_target} kJ/mol:")
                report += f"@ µ of {growth_rate_float} 1/h and minimal {target} of {min_target} kJ/mol:\n"
                minimization_result = perform_variable_minimization(
                    optmdfpathway_base_problem,
                    "reaction_error_sum",
                )

                num_changes = minimization_result['objective_value']
                print(f"  Minimal number of NADX/Y changes compared to in vivo distribution is {num_changes}")
                report += f"  Minimal number of NADX/Y changes compared to in vivo distribution is {num_changes}\n"
                print(num_changes_dict[full_scenario_key])
                num_changes_dict[full_scenario_key].append(num_changes)

                if num_changes > 1e-6:
                    print("  Affected reactions:")
                    report += "  Affected reactions:\n"
                    for var_id in minimization_result["values"]:
                        if not var_id.endswith("_error"):
                            continue
                        if minimization_result["values"][var_id] > 1e-6:
                            report += f" *{var_id}\n"
                            print(f" *{var_id}")

            print(num_changes_dict)
            for full_scenario_key in num_changes_dict.keys():
                num_changes = num_changes_dict[full_scenario_key]
                print(full_scenario_key)

                if ("OptMDF" in full_scenario_key) and ("aerobic" in full_scenario_key):
                    color="red"
                elif ("OptMDF" in full_scenario_key) and ("anaerobic" in full_scenario_key):
                    color="salmon"
                elif ("SubMDF" in full_scenario_key) and ("aerobic" in full_scenario_key):
                    color="blue"
                elif ("SubMDF" in full_scenario_key) and ("anaerobic" in full_scenario_key):
                    color="deepskyblue"

                plt.plot(
                    [float(x.replace(",", ".")) for x in growth_rates], # x
                    num_changes_dict[full_scenario_key], # y
                    label=full_scenario_key,
                    linestyle="-",
                    color=color,
                    linewidth=1.0,
                )
            figurename = f"4_{concentrations}.jpg"
            plt.legend(loc="best")
            plt.title("Minimal number of NADX/NADY swaps to reach given target")
            plt.xlabel("Growth rate [1/h]")
            plt.ylabel("Num changes")
            plt.xlim(min(growth_rates), max(growth_rates))
            plt.savefig(f"{figures_path}{figurename}")
            plt.close()
            print("")
            report += "\n"
    if disallowed_changed_reaction != "":
        disallowment_addition = "_no_"+disallowed_changed_reaction+"_change"
    else:
        disallowment_addition = ""
    with open(f"./cosa/results{suffix}/figures/4_report{disallowment_addition}.txt", "w") as f:
        f.write(report)


# LOGIC #
cosa_minimal_changes_test(anaerobic=False)
cosa_minimal_changes_test(anaerobic=False, c_source="acetate")
# cosa_minimal_changes_test(anaerobic=False, disallowed_changed_reaction="PDH_NADY")
# cosa_minimal_changes_test(anaerobic=False, disallowed_changed_reaction="NADH16pp_NADY")
# cosa_minimal_changes_test(anaerobic=False, disallowed_changed_reaction="FLDR2_NADX")
# cosa_minimal_changes_test(anaerobic=False, disallowed_changed_reaction="ICDHyr_FWD_NADX")
cosa_minimal_changes_test(anaerobic=True)
