"""In this script, the (NAD/NADH) / (NADP/NADPH) ratio of ratios variability function and its figure creation are definedd."""

# IMPORTS #
# External
from tkinter.messagebox import NO
import matplotlib.pyplot as plt
import cobra
import copy
import pulp
import math
# Internal
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


# PUBLIC FUNCTIONS SECTION #
def cosa_ratio_ratio_test(anaerobic: bool, expanded: bool, growth_epsilon: float = 0.01, c_source: str="glucose") -> None:
    """Perform (NAD/NADH) / (NADP/NADPH) ratio of ratios variability analysis under the given settings.

    Args:
        anaerobic (bool): Is is anaerobic (no oxygen)?
        expanded (bool): Is it a 2-cofactor (False) or 3-cofactor (True) model?
        growth_epsilon (float, optional): The ε (numerical value to go below) for the µ. Defaults to 0.01.
        c_source (str, optional): Either 'glucose' or 'acetate'. Defaults to "glucose".
    """
    suffix = cosa_get_suffix(anaerobic, expanded, c_source)
    figures_path = f"./cosa/results{suffix}/figures/"
    ensure_folder_existence(figures_path)
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=expanded, c_source=c_source)

    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    ratio_ratios = [
        (("nadh_tcosa_c", "nad_tcosa_c"), ("nadph_tcosa_c", "nadp_tcosa_c")),
        # (("nadph_tcosa_c", "nadp_tcosa_c"), ("nadh_tcosa_c", "nad_tcosa_c")),
    ]
    report = ""
    ratio_ratio_test_data = {}
    original_cobra_model = copy.deepcopy(cobra_model)

    if (c_source != "glucose") or (anaerobic) or (expanded):
        concentration_scenarios = ("STANDARDCONC",)
    else:
        concentration_scenarios = ("STANDARDCONC", "VIVOCONC",)

    for concentrations in concentration_scenarios:
        print(f"=CONCENTRATION RANGES: {concentrations}=")
        report += f"=CONCENTRATION RANGES: {concentrations}=\n"
        if concentrations == "STANDARDCONC":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
            used_concentration_values = concentration_values_free
            titleaddition = ""
        elif concentrations == "VIVOCONC":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper
            titleaddition = "\nwith adapted measured concentration ranges"

        for target in ("OPTMDF", "OPTSUBMDF"):
            print(f"===OPTIMIZATION TARGET: {target}===")
            report += f"===OPTIMIZATION TARGET: {target}===\n"
            cobra_model = copy.deepcopy(original_cobra_model)
            cobra_model = cosa_get_model_with_nadx_scenario(
                nadx_scenario="WILDTYPE",
                cobra_model=cobra_model,
            )

            jsondata_invivo = json_zip_load(f"cosa/results{suffix}/runs/{target}_{concentrations}_WILDTYPE.json")

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
            for ratio_ratio in ratio_ratios:
                ylabel = ""
                title = f'Ratio ratio range of ({ratio_ratio[0][0].replace("_tcosa_c", "").upper()}:{ratio_ratio[0][1].replace("_tcosa_c", "").upper()}) '\
                        f'to ({ratio_ratio[1][0].replace("_tcosa_c", "").upper()}:{ratio_ratio[1][1].replace("_tcosa_c", "").upper()})'
                title += " at growth-rate-dependent\nmaximal "
                if target == "OPTMDF":
                    title += "OptMDF"
                else:
                    title += "OptSubMDF"
                title += f" under {'anaerobic' if anaerobic else 'aerobic'} conditions\n"
                title += titleaddition
                ylabel = "Ratio"

                current_ratio += 1
                ratio_ratio_var = pulp.LpVariable(
                    name=f"ratio_ratio_var_{current_ratio}",
                    cat=pulp.LpContinuous,
                )
                # Set ratio of ratios expression in linear form, i.e.,
                # ln((a/a')/(b/b')) = ln(a) - ln(a') - ln(b) + ln(b')
                optmdfpathway_base_problem +=\
                    optmdfpathway_base_variables["x_"+ratio_ratio[0][0]] - optmdfpathway_base_variables["x_"+ratio_ratio[0][1]] - optmdfpathway_base_variables["x_"+ratio_ratio[1][0]] + optmdfpathway_base_variables["x_"+ratio_ratio[1][1]] == ratio_ratio_var

                report += title
                print(title)
                min_ratios: List[float] = []
                max_ratios: List[float] = []
                growth_rates = jsondata_invivo.keys()
                figurename = f'2C_{ratio_ratio[0][0].replace("_tcosa_c", "").upper()}_to_{ratio_ratio[0][1].replace("_tcosa_c", "").upper()}___to___'\
                             f'{ratio_ratio[1][0].replace("_tcosa_c", "").upper()}_to_{ratio_ratio[1][1].replace("_tcosa_c", "")}_{target}_{concentrations}.jpg'

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

                    print(growth_rate)

                    print(f" @ µ [1/h] of {growth_rate_float} and min {target} of {min_target} kJ/mol")
                    report += f" @ µ [1/h] of {growth_rate_float} and min {target} of {min_target} kJ/mol\n"


                    minimization_result = perform_variable_minimization(
                        optmdfpathway_base_problem,
                        ratio_ratio_var.name,
                    )
                    maximization_result = perform_variable_maximization(
                        optmdfpathway_base_problem,
                        ratio_ratio_var.name,
                    )
                    print(f" Min ratio ratio: {math.exp(minimization_result['values'][ratio_ratio_var.name])}")
                    report += f" Min ratio ratio: {math.exp(minimization_result['values'][ratio_ratio_var.name])}\n"
                    print(f" Max ratio ratio: {math.exp(maximization_result['values'][ratio_ratio_var.name])}")
                    report += f" Max ratio ratio: {math.exp(maximization_result['values'][ratio_ratio_var.name])}\n"
                    min_ratios.append(math.exp(minimization_result['values'][ratio_ratio_var.name]))
                    max_ratios.append(math.exp(maximization_result['values'][ratio_ratio_var.name]))

                plotted_growth_rates = [x.replace(",", ".") for x in growth_rates]
                ratio_ratio_test_data[figurename] = {
                    "min_ratios": min_ratios,
                    "max_ratios": max_ratios,
                    "plotted_growth_rates": plotted_growth_rates,
                }
                json_write(f"./cosa/results{suffix}/ratio_ratio_test_data.json", ratio_ratio_test_data)

    with open(f"./cosa/results{suffix}/figures/2C_report.txt", "w", encoding="utf-8") as f:
        f.write(report)

    ratio_ratio_test_data = json_load(f"./cosa/results{suffix}/ratio_ratio_test_data.json")

    for figurename in ratio_ratio_test_data.keys():
        plotted_growth_rates = ratio_ratio_test_data[figurename]["plotted_growth_rates"]
        min_ratios = ratio_ratio_test_data[figurename]["min_ratios"]
        max_ratios = ratio_ratio_test_data[figurename]["max_ratios"]
        plt.plot(
            plotted_growth_rates, # x
            min_ratios, # y
            "bo",
            label="Minimal ratio of ratios",
            linewidth=1.0,
        )
        plt.plot(
            plotted_growth_rates, # x
            max_ratios, # y
            "ro",
            label="Maximal ratio of ratios",
            linewidth=1.0,
        )
        plt.legend(loc="best")
        plt.xlabel("Growth rate [1/h]")
        plt.ylabel("OptSubMDF [kJ/mol]" if "OPTSUBMDF" in figurename else "OptMDF [kJ/mol]")
        plt.xlim(min(plotted_growth_rates), max(plotted_growth_rates))
        plt.savefig(f"{figures_path}{figurename}")
        plt.close()


def get_latex_scientific_notation(x: float):
    in_notation = str('{:.1E}'.format(x))
    print(in_notation)
    if "E-" in in_notation:
        base = in_notation.split("E-")[0]
        exponent = in_notation.split("E-")[1]
        return r"$\mathrm{"+f"{int(float(base))} "+r" \cdot "+" 10^"+"{-"+f"{int(exponent)}"+"}"+r"}$"
    else:
        return "0.0"


def cosa_create_full_ratio_ratio_test_figure_four_panels():
    """Creates the metabolite ratio of ratios figure in TCOSA's publication."""
    ratio_ratio_test_data_aerobic = json_load("cosa/results_aerobic/ratio_ratio_test_data.json")
    ratio_ratio_test_data_anaerobic = json_load("cosa/results_anaerobic/ratio_ratio_test_data.json")
    concentrations = ("STANDARDCONC",) #"VIVOCONC")
    for concentration in concentrations:
        figurenames_to_plots = {
            ("Aerob", "OptMDF", "A", f"2C_NADH_to_NAD___to___NADPH_to_nadp_OPTMDF_{concentration}.jpg"): (0, 0),
            ("Aerob", "OptSubMDF", "B", f"2C_NADH_to_NAD___to___NADPH_to_nadp_OPTSUBMDF_{concentration}.jpg"): (0, 1),
            ("Anaerob", "OptMDF", "C", f"2C_NADH_to_NAD___to___NADPH_to_nadp_OPTMDF_{concentration}.jpg"): (1, 0),
            ("Anaerob", "OptSubMDF", "D", f"2C_NADH_to_NAD___to___NADPH_to_nadp_OPTSUBMDF_{concentration}.jpg"): (1, 1),
        }
        first = True

        fig, axs = plt.subplots(nrows=2, ncols=2, dpi=500, figsize=(19, 10)) #sharex=True, figsize=(50, 25), dpi=120, facecolor="white")
        fig.tight_layout(pad=3.75)
        for figurename_tuple in figurenames_to_plots.keys():
            if figurename_tuple[0] == "Aerob":
                ratio_ratio_test_data = ratio_ratio_test_data_aerobic
                is_aerobic = True
            else:
                ratio_ratio_test_data = ratio_ratio_test_data_anaerobic
                is_aerobic = False

            if first:
                min_label = "Minimales Verhältnis"
                max_label = "Maximales Verhältnis"
                first = False
            else:
                min_label = None
                max_label = None

            title = f"{figurename_tuple[2]} {figurename_tuple[0]}, {figurename_tuple[1].replace('Opt', '')}"

            figurename = figurename_tuple[3]
            if is_aerobic:
                plotted_growth_rates = ratio_ratio_test_data[figurename]["plotted_growth_rates"][:-1]
                axs_index = figurenames_to_plots[figurename_tuple]
                min_ratios = ratio_ratio_test_data[figurename]["min_ratios"][:-1]
                max_ratios = ratio_ratio_test_data[figurename]["max_ratios"][:-1]
            else:
                plotted_growth_rates = ratio_ratio_test_data[figurename]["plotted_growth_rates"][:-1]
                axs_index = figurenames_to_plots[figurename_tuple]
                min_ratios = ratio_ratio_test_data[figurename]["min_ratios"][:-1]
                max_ratios = ratio_ratio_test_data[figurename]["max_ratios"][:-1]
            plotted_growth_rates = [float(x) for x in plotted_growth_rates]
            min_ratios = [float(x) for x in min_ratios]
            max_ratios = [float(x) for x in max_ratios]
            axs[axs_index].plot(
                plotted_growth_rates[::-1], # x
                min_ratios[::-1], # y
                "bo",
                label=min_label,
                linewidth=1.0,
            )
            axs[axs_index].plot(
                plotted_growth_rates[::-1], # x
                max_ratios[::-1], # y
                "ro",
                label=max_label,
                linewidth=1.0,
            )
            import matplotlib
            axs[axs_index].set_title(title, loc="left", fontweight="bold", fontsize=17)
            if figurename_tuple[2] == "A":
                axs[axs_index].set_xlim(0.025, 0.895)
                axs[axs_index].set_ylim(-.000003, 0.00006)
                axs[axs_index].yaxis.set_major_formatter(
                    matplotlib.ticker.FuncFormatter(lambda x, p: '{:.0E}'.format(x))
                )
            elif figurename_tuple[2] == "B":
                axs[axs_index].set_xlim(0.025, 0.895)
                axs[axs_index].set_ylim(-.0000004, 0.000005)
                axs[axs_index].yaxis.set_major_formatter(
                    matplotlib.ticker.FuncFormatter(lambda x, p: '{:.0E}'.format(x))
                )
            elif figurename_tuple[2] == "C":
                axs[axs_index].set_ylim(-.003, 0.05)
                axs[axs_index].yaxis.set_major_formatter(
                    matplotlib.ticker.FuncFormatter(lambda x, p: round(x, 2))
                )
            elif figurename_tuple[2] == "D":
                axs[axs_index].set_ylim(-.04, 1.0)
            axs[axs_index].set_xlabel("Wachstumsrate [1/h]", fontsize=16)
            # axs[axs_index].set_ylabel(r"$\mathrm{\frac{[NADH]/[NAD^{+}]}{[NADPH]/[NADP^{+}]}}$", fontsize=16)
            axs[axs_index].set_ylabel(r"$\mathrm{[NADH]/[NAD^{+}] \ / \ [NADPH]/[NADP^{+}]}$", fontsize=12)
            axs[axs_index].tick_params(labelsize=13)
        fig.legend(loc=(0.16, 0.9525), ncol=2, fontsize=17)
        # fig.subplots_adjust(right=1.25)

        fig.savefig(f"./cosa/full_ratio_ratio_test_figure_{concentration}.png", bbox_inches='tight', pad_inches=0.05)
        plt.close()


# cosa_ratio_ratio_test(anaerobic=False, expanded=False)
# cosa_ratio_ratio_test(anaerobic=True, expanded=False)
# cosa_ratio_ratio_test(anaerobic=False, expanded=False, c_source="acetate")

cosa_create_full_ratio_ratio_test_figure_four_panels()
