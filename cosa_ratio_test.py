import matplotlib.pyplot as plt
import cobra
import copy
import pulp
import math
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from cosa_get_suffix import cosa_get_suffix
from helper import json_zip_load
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


def cosa_ratio_test(anaerobic: bool, expanded: bool, growth_epsilon: float = 0.01, fixed_concentrations=()) -> None:
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=expanded)

    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    suffix = cosa_get_suffix(anaerobic, expanded)

    figures_path = f"./cosa/results{suffix}/figures/"
    ensure_folder_existence(figures_path)

    ratios = [
        # ("nadh_tcosa_c", "nad_tcosa_c"),
        # ("nadph_tcosa_c", "nadp_tcosa_c"),
        ("nadp_tcosa_c", "nadph_tcosa_c"),
    ]

    if expanded:
        ratios += [
            ("nad_tcosa_c", "nadzh_tcosa_c"),
            ("nad_tcosa_c", "nadz_tcosa_c"),
            ("nadp_tcosa_c", "nadz_tcosa_c"),
            ("nadh_tcosa_c", "nadzh_tcosa_c"),
            ("nadph_tcosa_c", "nadzh_tcosa_c"),
        ]

    report = ""
    original_cobra_model = copy.deepcopy(cobra_model)
    for concentrations in ("STANDARDCONCS", "PAPERCONCS"):
        print(f"=CONCENTRATION RANGES: {concentrations}=")
        report += f"=CONCENTRATION RANGES: {concentrations}=\n"
        if concentrations == "STANDARDCONCS":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
            used_concentration_values = concentration_values_free
            titleaddition = ""
        elif concentrations == "PAPERCONCS":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper
            titleaddition = "\nwith adapted measured concentration ranges"

        for fixed_concentration in fixed_concentrations:
            met_id = fixed_concentration[0]
            lb = fixed_concentration[1]
            ub = fixed_concentration[2]
            used_concentration_values[met_id] = {}
            used_concentration_values[met_id]["min"] = lb
            used_concentration_values[met_id]["max"] = ub

        for target in ("OPTSUBMDF", "OPTMDF"):
            print(f"===OPTIMIZATION TARGET: {target}===")
            report += f"===OPTIMIZATION TARGET: {target}===\n"
            cobra_model = copy.deepcopy(original_cobra_model)
            cobra_model = cosa_get_model_with_nadx_scenario(
                nadx_scenario="WILDTYPE",
                cobra_model=cobra_model,
            )

            jsondata_invivo = json_zip_load(f"cosa/results{suffix}/runs/{target}_VIVOCONC_WILDTYPE.json")

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

            current_ratio = 0
            for ratio in ratios:
                ylabel = ""
                title = f'Ratio range of {ratio[0].replace("_c", "").upper()} to {ratio[1].replace("_c", "").upper()}'
                title += " at growth-rate-dependent\nmaximal "
                if target == "OPTMDF":
                    title += "OptMDF"
                else:
                    title += "OptSubMDF"
                title += f" under {'anaerobic' if anaerobic else 'aerobic'} conditions"
                title += titleaddition
                ylabel = "Ratio"

                current_ratio += 1
                ratio_var = pulp.LpVariable(
                    name=f"ratio_var_{current_ratio}",
                    cat=pulp.LpContinuous,
                )
                optmdfpathway_base_problem +=\
                    optmdfpathway_base_variables["x_"+ratio[0]] - optmdfpathway_base_variables["x_"+ratio[1]] == ratio_var

                print(f">RATIO OF {ratio[0]} to {ratio[1]}")
                report += f">RATIO OF {ratio[0]} to {ratio[1]}\n"
                min_ratios: List[float] = []
                max_ratios: List[float] = []
                growth_rates = jsondata_invivo.keys()
                figurename = f'2B_{ratio[0].replace("_tcosa_c", "").upper()}_to_{ratio[1].replace("_tcosa_c", "").upper()}_{target}_{concentrations}_{fixed_concentrations}.jpg'
                for growth_rate in growth_rates:
                    growth_rate_float = float(growth_rate.replace(",", "."))
                    optmdfpathway_base_variables[biomass_reaction_id].bounds(
                        growth_rate_float-growth_epsilon,
                        1e12,
                    )

                    if target == "OPTMDF":
                        min_target = jsondata_invivo[growth_rate]["values"]["var_B"]
                        optmdfpathway_base_variables["var_B"].bounds(
                            min_target,
                            1e12,
                        )
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

                    print(f" @ µ [1/h] of {growth_rate_float} and min {target} of {min_target} kJ/mol")
                    report += f" @ µ [1/h] of {growth_rate_float} and min {target} of {min_target} kJ/mol\n"


                    minimization_result = perform_variable_minimization(
                        optmdfpathway_base_problem,
                        ratio_var.name,
                    )
                    maximization_result = perform_variable_maximization(
                        optmdfpathway_base_problem,
                        ratio_var.name,
                    )
                    print(f" Min ratio: {math.exp(minimization_result['values'][ratio_var.name])}")
                    report += f" Min ratio: {math.exp(minimization_result['values'][ratio_var.name])}\n"
                    print(f" Max ratio: {math.exp(maximization_result['values'][ratio_var.name])}")
                    report += f" Max ratio: {math.exp(maximization_result['values'][ratio_var.name])}\n"
                    min_ratios.append(math.exp(minimization_result['values'][ratio_var.name]))
                    max_ratios.append(math.exp(maximization_result['values'][ratio_var.name]))

                plotted_growth_rates = [x.replace(",", ".") for x in growth_rates]
                plt.plot(
                    plotted_growth_rates, # x
                    min_ratios, # y
                    "bo",
                    label="Minimal ratio",
                    color="red",
                    linewidth=1.0,
                )
                plt.plot(
                    plotted_growth_rates, # x
                    max_ratios, # y
                    "ro",
                    label="Maximal ratio",
                    color="blue",
                    linewidth=1.0,
                )
                plt.legend(loc="best")
                plt.title(title)
                plt.xlabel("Growth rate [1/h]")
                plt.ylabel(ylabel)
                plt.xlim(min(plotted_growth_rates), max(plotted_growth_rates))
                plt.savefig(f"{figures_path}{figurename}")
                plt.close()
                print("")

    with open(f"./cosa/results{suffix}/figures/2B_report_{fixed_concentrations}.txt", "w", encoding="utf-8") as f:
        f.write(report)


cosa_ratio_test(anaerobic=False, expanded=False)
cosa_ratio_test(anaerobic=True, expanded=False)
# cosa_ratio_test(anaerobic=False, expanded=False, fixed_concentrations=[("nad_tcosa_c", 2.32E-3, 2.80E-3), ("nadh_tcosa_c", 5.45E-5, 1.27E-4)])
# cosa_ratio_test(anaerobic=True, expanded=False, fixed_concentrations=[("nadp_tcosa_c", 1.40E-7, 3.11E-5), ("nadph_tcosa_c", 1.10E-4, 1.34E-4)])
