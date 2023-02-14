import matplotlib.pyplot as plt
import cobra
import copy
import pulp
import math
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from cosa_get_suffix import cosa_get_suffix
from helper import json_zip_load
from typing import List
from optmdfpathway import (
    STANDARD_R, STANDARD_T, get_optmdfpathway_base_problem,
    add_differential_reactions_constraints, get_z_variable_status,
    perform_concentration_variability_analysis
)
from optimization import perform_variable_minimization, perform_variable_maximization
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from typing import Dict
from helper import ensure_folder_existence


def cosa_cva_test(anaerobic: bool, expanded: bool, growth_epsilon: float = 0.01) -> None:
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=expanded)

    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    suffix = cosa_get_suffix(anaerobic, expanded)

    figures_path = f"./cosa/results{suffix}/figures/"
    ensure_folder_existence(figures_path)

    metabolites = [
        "nadx_c", "nadxh_c",
        "nady_c", "nadyh_c",
    ]

    if expanded:
        metabolites += ["nadz_c", "nadzh_c"]

    report = ""
    original_cobra_model = copy.deepcopy(cobra_model)
    for concentrations in ("STANDARDCONCS", "PAPERCONCS"):
        print(f"=CONCENTRATION RANGES: {concentrations}=")
        report += f"=CONCENTRATION RANGES: {concentrations}=\n"

        if concentrations == "STANDARDCONCS":
           dG0_values = copy.deepcopy(standardconc_dG0_values)
           used_concentration_values = concentration_values_free
        elif concentrations == "PAPERCONCS":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper

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
                sub_network_ids=nad_reactions+nadp_reactions,
            )
            optmdfpathway_base_variables: Dict[str, pulp.LpVariable] = optmdfpathway_base_problem.variablesDict()

            for metabolite in metabolites:
                concentration_var_name = f"x_{metabolite}"

                ylabel = ""
                title = f'Min/max concentrations of {metabolite}'
                title += " at growth-rate-dependent maximal "
                if target == "OPTMDF":
                    title += "OptMDF"
                else:
                    title += "SubMDF"
                ylabel = "Concentration [M]"

                print(f">MIN/MAX CONCENTRATIONS OF {metabolite}")
                report += f">MIN/MAX CONCENTRATIONS OF {metabolite}\n"
                min_concentrations: List[float] = []
                max_concentrations: List[float] = []
                growth_rates = jsondata_invivo.keys()
                figurename = f'2A_{metabolite}_{target}_{concentrations}.jpg'
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
                        concentration_var_name,
                    )
                    maximization_result = perform_variable_maximization(
                        optmdfpathway_base_problem,
                        concentration_var_name,
                    )
                    print(f" Min concentration: {math.exp(minimization_result['values'][concentration_var_name])}")
                    report += f" Min concentration: {math.exp(minimization_result['values'][concentration_var_name])}\n"
                    print(f" Max concentration: {math.exp(maximization_result['values'][concentration_var_name])}")
                    report += f" Max concentration: {math.exp(maximization_result['values'][concentration_var_name])}\n"
                    min_concentrations.append(math.exp(minimization_result['values'][concentration_var_name]))
                    max_concentrations.append(math.exp(maximization_result['values'][concentration_var_name]))

                plt.plot(
                    growth_rates, # x
                    min_concentrations, # y
                    label="Minimal concentration",
                    linestyle="-",
                    color="red",
                    linewidth=1.0,
                )
                plt.plot(
                    growth_rates, # x
                    max_concentrations, # y
                    label="Maximal concentration",
                    linestyle="-",
                    color="blue",
                    linewidth=1.0,
                )
                plt.legend(loc="best")
                plt.title(title)
                plt.xlabel("Growth rate [1/h]")
                plt.ylabel(ylabel)
                plt.xlim(min(growth_rates), max(growth_rates))
                plt.savefig(f"{figures_path}{figurename}")
                plt.close()

                print("")

    with open(f"./cosa/results{suffix}/figures/2A_report.txt", "w", encoding="utf-8") as f:
        f.write(report)


cosa_cva_test(anaerobic=False, expanded=False)
cosa_cva_test(anaerobic=True, expanded=False)
# cosa_cva_test(anaerobic=False, expanded=True)
# cosa_cva_test(anaerobic=True, expanded=True)
