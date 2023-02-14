import cobra
import copy
import matplotlib.pyplot as plt
import pulp
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from cosa_get_suffix import cosa_get_suffix
from helper import json_zip_load, ensure_folder_existence
from optmdfpathway import (
    STANDARD_R, STANDARD_T, get_optmdfpathway_base_problem,
    add_differential_reactions_constraints, get_z_variable_status,
)
from optimization import perform_variable_minimization, perform_variable_maximization
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from typing import Dict


def cosa_optimal_test(anaerobic: bool) -> None:
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=False)

    suffix = cosa_get_suffix(anaerobic, False)

    figures_path = f"./cosa/results{suffix}/figures/"
    ensure_folder_existence(figures_path)

    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    for concentrations in ("STANDARDCONCS", ):  # "PAPERCONCS",
        print(f"=CONCENTRATION RANGES: {concentrations}=")
        if concentrations == "STANDARDCONCS":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
        elif concentrations == "PAPERCONCS":
            dG0_values = copy.deepcopy(paperconc_dG0_values)

        if concentrations == "PAPERCONCS":
            used_concentration_values = concentration_values_paper
        else:
            used_concentration_values = concentration_values_free
        for target in ("OPTMDF", "OPTSUBMDF"):
            print(f"===OPTIMIZATION TARGET: {target}===")
            ylabel = ""
            title = f'Optimal '
            if target == "OPTMDF":
                title += "OptMDF"
                ylabel += "OptMDF"
            else:
                title += "SubMDF"
                ylabel += "SubMDF"
            title += " values of various theoretically maximal distributions"
            ylabel += "[kJ/mol]"
            if concentrations == "PAPERCONCS":
                jsondata_optimal = json_zip_load(f"cosa/results{suffix}/runs/{target}_VIVOCONC_FLEXIBLE.json")
            else:
                jsondata_optimal = json_zip_load(f"cosa/results{suffix}/runs/{target}_STANDARDCONC_FLEXIBLE.json")

            growth_rates = list(jsondata_optimal.keys())
            for mu in jsondata_optimal.keys():
                print(f"==USED GROWTH RATE DISTRUBUTION: {mu.replace(',', '.')} 1/h==")
                z_status_optimal = get_z_variable_status(jsondata_optimal[mu], "_COSA_")

                z_status_ids = set([x.replace("_NADX", "").replace("_NADY", "") for x in list(z_status_optimal.keys())])
                z_status_refined = {}
                for z_status_id in z_status_ids:
                    one_id = z_status_id+"_NADX"
                    two_id = z_status_id+"_NADY"
                    if z_status_optimal[one_id] + z_status_optimal[two_id] == 0:
                        continue
                    z_status_refined[one_id] = z_status_optimal[one_id]
                    z_status_refined[two_id] = z_status_optimal[two_id]
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

                for z_var_id in z_status_refined.keys():
                    if z_status_refined[z_var_id] == 0:
                        optmdfpathway_base_problem += optmdfpathway_base_variables[f"{z_var_id}"] <= 1e-6

                optimal_values = []
                reached_values = []
                for growth_rate in growth_rates:
                    growth_rate_float = float(growth_rate.replace(",", "."))
                    print(f"@ µ={growth_rate_float} 1/h")

                    optmdfpathway_base_variables[biomass_reaction_id].bounds(
                        growth_rate_float,
                        1e12,
                    )
                    if target == "OPTSUBMDF":
                        optimal_target = jsondata_optimal[growth_rate]["values"]["var_B2"]
                        optmdfpathway_base_variables["var_B"].bounds(
                            MIN_OPTMDF,
                            1e12,
                        )
                        maximization_result = perform_variable_maximization(
                            optmdfpathway_base_problem,
                            "var_B2",
                        )
                    else:
                        optimal_target = jsondata_optimal[growth_rate]["values"]["var_B"]
                        maximization_result = perform_variable_maximization(
                            optmdfpathway_base_problem,
                            "var_B",
                        )
                    print(f" Theoretically optimal {target}:                                  {optimal_target} kJ/mol")
                    print(f" Reached {target} with optimal NADX/Y distribution @ µ={mu} 1/h: {-maximization_result['objective_value']} kJ/mol")
                    reached_values.append(-maximization_result['objective_value'])
                    optimal_values.append(optimal_target)

                print("")
                plt.plot(
                    growth_rates, # x
                    reached_values, # y
                    label=f"Distribution of {mu.replace(',', '.')} 1/h",
                    linestyle="-",
                    # color="blue",
                    linewidth=1.0,
                )
            figurename = f"3_{target}_{concentrations}.jpg"
            plt.plot(
                growth_rates, # x
                optimal_values, # y
                label=f"Theoretical maximum",
                linestyle="--",
                color="black",
                linewidth=1.25,
            )
            plt.legend(loc="best")
            plt.title(title)
            plt.xlabel("Growth rate [1/h]")
            plt.ylabel(ylabel)
            plt.xlim(min(growth_rates), max(growth_rates))
            plt.savefig(f"{figures_path}{figurename}")
            plt.close()

cosa_optimal_test(anaerobic=False)
cosa_optimal_test(anaerobic=True)
