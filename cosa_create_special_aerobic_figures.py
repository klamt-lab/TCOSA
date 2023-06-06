"""In this script, the following figures are generated:

1) The figure for the solutions with in vivo concentraitons from Bennett et al., 2009.
2) The figure with results under acetate.

Both figures depict the situation under aerobic conditions as this is the only repective viable environment.
"""

# IMPORTS #
# External
import matplotlib.pyplot as plt
import pandas
# Internal
from helper import json_load


# PUBLIC FUNCTIONS #
def create_in_vivo_concentrations_figure():
    """Creates figure for results with in vivo concentration ranges."""
    ratio_ratio_test_data_aerobic = json_load("cosa/results_aerobic/ratio_ratio_test_data.json")
    concentration = "VIVOCONC"
    output_path = "./cosa/in_vivo_concentrations_figure.png"
    table_path = f"cosa/results_aerobic/optsubmdf_table_{concentration}.csv"
    pad = 4.99

    target = "OPTSUBMDF"

    figurename_tuple = ("aerobic", f"2C_NADH_to_NAD___to___NADPH_to_nadp_{target}_{concentration}.jpg")

    fig, axs = plt.subplots(nrows=1, ncols=2, dpi=500, figsize=(18, 6)) #sharex=True, figsize=(50, 25), dpi=120, facecolor="white")
    fig.tight_layout(pad=pad)

    ########################################################
    ########################################################
    ########################################################
    # 0: Random sampling figure
    str_to_float = lambda x: float(x.replace(",", "."))
    list_to_float = lambda x: [str_to_float(y) for y in x]

    growth_rate_id = "µ [1/h]"
    best_id = "FLEXIBLE"
    in_vivo_id = "WILDTYPE"
    only_one_id = "SINGLE_COFACTOR"

    table = pandas.read_csv(
        table_path,
        sep="\t",
    )

    headers = list(table.keys())
    del(headers[headers.index(best_id)])
    del(headers[headers.index(only_one_id)])
    del(headers[headers.index(in_vivo_id)])
    headers += [only_one_id, in_vivo_id, best_id]

    growth_rates = list_to_float(table[growth_rate_id])
    is_first_random = True
    for header in headers:
        if header == growth_rate_id:
            continue
        elif header == best_id:
            label = "Flexible specificity"
            linestyle = "--"
            color = "yellowgreen"
            linewidth = 2.0
        elif header == in_vivo_id:
            label = "Wild-type specificity"
            linestyle = "-"
            color = "black"
            linewidth = 2.0
        elif header == only_one_id:
            label = "Single cofactor pool"
            linestyle = "-"
            color = "red"
            linewidth = 2.0
        else:
            if is_first_random:
                label = "Random specificities"
                is_first_random = False
            else:
                label = ""
            linestyle = "-"
            color = "paleturquoise"
            linewidth = 1.0

        axs[0].plot(
            growth_rates[:-1], # x
            list_to_float(table[header])[:-1], # y
            label=label,
            linestyle=linestyle,
            color=color,
            linewidth=linewidth,
        )
    axs[0].legend(loc="lower left", fontsize=16)
    axs[0].set_title("A", loc="left", fontweight="bold", fontsize=18)
    axs[0].set_xlabel("Growth rate [1/h]", fontsize=14)
    axs[0].set_ylabel("SubMDF [kJ/mol]", fontsize=14)
    axs[0].set_xlim(min(growth_rates[:-1]), max(growth_rates[:-1]))
    axs[0].tick_params(labelsize=14)


    ########################################################
    ########################################################
    ########################################################
    # 1: Ratio ratio figure
    ratio_ratio_test_data = ratio_ratio_test_data_aerobic

    min_label = "Minimal ratio"
    max_label = "Maximal ratio"

    figurename = figurename_tuple[1]
    plotted_growth_rates = [float(x) for x in ratio_ratio_test_data[figurename]["plotted_growth_rates"][:-1]]
    min_ratios = [float(x) for x in ratio_ratio_test_data[figurename]["min_ratios"][:-1]]
    max_ratios = [float(x) for x in ratio_ratio_test_data[figurename]["max_ratios"][:-1]]
    axs[1].plot(
        plotted_growth_rates[::-1], # x
        min_ratios[::-1], # y
        "bo",
        label=min_label,
        linewidth=1.0,
    )
    axs[1].plot(
        plotted_growth_rates[::-1], # x
        max_ratios[::-1], # y
        "ro",
        label=max_label,
        linewidth=1.0,
    )
    axs[1].legend(loc="upper center", ncol=2, fontsize=16)
    axs[1].set_title("B", loc="left", fontweight="bold", fontsize=18)
    axs[1].set_xlabel("Growth rate [1/h]", fontsize=15)
    # axs[1].set_ylabel(r"$\mathrm{\frac{[NADH]/[NAD^{+}]}{[NADPH]/[NADP^{+}]}}$", fontsize=15)
    axs[1].set_ylabel(r"$\mathrm{([NADH]/[NAD^{+}])/([NADPH]/[NADP^{+}])}$", fontsize=14)
    axs[1].tick_params(labelsize=14)

    axs[1].set_ylim(-.00003, 0.0006)
    axs[1].set_xlim(0.025, 0.89)
    pad_inches = 0.05

    # fig.subplots_adjust(right=1.1)

    fig.savefig(output_path, bbox_inches='tight', pad_inches=pad_inches)
    plt.close()


def create_acetate_figure():
    """Creates figure for results under acetate."""
    ratio_ratio_test_data_aerobic = json_load("cosa/results_aerobic_acetate/ratio_ratio_test_data.json")
    concentration = "STANDARDCONC"
    output_path = "./cosa/acetate_figure.png"
    pad = 5.25

    fig, axs = plt.subplots(nrows=2, ncols=2, dpi=500, figsize=(14, 11)) #sharex=True, figsize=(50, 25), dpi=120, facecolor="white")
    fig.tight_layout(pad=pad)

    for target in ("OPTSUBMDF", "OPTMDF"):
        if target == "OPTSUBMDF":
            xpos = 1
        else:
            xpos = 0
        table_path = f"cosa/results_aerobic_acetate/{target.lower()}_table_{concentration}.csv"

        figurename_tuple = ("aerobic", f"2C_NADH_to_NAD___to___NADPH_to_nadp_{target}_{concentration}.jpg")
        ########################################################
        ########################################################
        ########################################################
        # 0: Random sampling figures
        str_to_float = lambda x: float(x.replace(",", "."))
        list_to_float = lambda x: [str_to_float(y) for y in x]

        growth_rate_id = "µ [1/h]"
        best_id = "FLEXIBLE"
        in_vivo_id = "WILDTYPE"
        only_one_id = "SINGLE_COFACTOR"

        table = pandas.read_csv(
            table_path,
            sep="\t",
        )

        headers = list(table.keys())
        del(headers[headers.index(best_id)])
        del(headers[headers.index(only_one_id)])
        del(headers[headers.index(in_vivo_id)])
        headers += [only_one_id, in_vivo_id, best_id]

        growth_rates = list_to_float(table[growth_rate_id])
        is_first_random = True
        for header in headers:
            if header == growth_rate_id:
                continue
            elif header == best_id:
                label = "Flexible specificity"
                linestyle = "--"
                color = "yellowgreen"
                linewidth = 2.0
            elif header == in_vivo_id:
                label = "Wild-type specificity"
                linestyle = "-"
                color = "black"
                linewidth = 2.0
            elif header == only_one_id:
                label = "Single cofactor pool"
                linestyle = "-"
                color = "red"
                linewidth = 2.0
            else:
                if is_first_random:
                    label = "Random specificities"
                    is_first_random = False
                else:
                    label = ""
                linestyle = "-"
                color = "paleturquoise"
                linewidth = 1.0

            axs[xpos, 0].plot(
                growth_rates[:-1], # x
                list_to_float(table[header])[:-1], # y
                label=label,
                linestyle=linestyle,
                color=color,
                linewidth=linewidth,
            )
        axs[xpos, 0].legend(loc="lower left", fontsize=14)
        axs[xpos, 0].set_title("C Sampling with SubMDF" if "SUB" in target else "A Sampling with MDF", loc="left", fontweight="bold", fontsize=18)
        axs[xpos, 0].set_xlabel("Growth rate [1/h]", fontsize=14)
        axs[xpos, 0].set_ylabel(f"{target.replace('OPT', '').replace('SUB', 'Sub')} [kJ/mol]", fontsize=14)
        axs[xpos, 0].set_xlim(min(growth_rates[:-1]), max(growth_rates[:-1]))
        axs[xpos, 0].tick_params(labelsize=14)


        ########################################################
        ########################################################
        ########################################################
        # 1: Ratio ratio figures
        ratio_ratio_test_data = ratio_ratio_test_data_aerobic

        min_label = "Minimal ratio"
        max_label = "Maximal ratio"

        figurename = figurename_tuple[1]
        plotted_growth_rates = [float(x) for x in ratio_ratio_test_data[figurename]["plotted_growth_rates"][:-1]]
        min_ratios = [float(x) for x in ratio_ratio_test_data[figurename]["min_ratios"][:-1]]
        max_ratios = [float(x) for x in ratio_ratio_test_data[figurename]["max_ratios"][:-1]]
        axs[xpos, 1].plot(
            plotted_growth_rates[::-1], # x
            min_ratios[::-1], # y
            "bo",
            label=min_label,
            linewidth=1.0,
        )
        axs[xpos, 1].plot(
            plotted_growth_rates[::-1], # x
            max_ratios[::-1], # y
            "ro",
            label=max_label,
            linewidth=1.0,
        )
        axs[xpos, 1].legend(loc="upper center", ncol=2, fontsize=15)
        axs[xpos, 1].set_title("D Ratios with SubMDF" if "SUB" in target else "B Ratios with MDF", loc="left", fontweight="bold", fontsize=18)
        axs[xpos, 1].set_xlabel("Growth rate [1/h]", fontsize=15)
        axs[xpos, 1].set_ylabel(r"$\mathrm{[NADH]/[NAD^{+}] \ / \ [NADPH]/[NADP^{+}]}$", fontsize=13)
        axs[xpos, 1].tick_params(labelsize=14)
        axs[xpos, 1].set_ylim(-.0001 if "SUB" in target else -0.005, 0.0015 if "SUB" in target else 0.09)
        axs[xpos, 1].set_xlim(0.045, 0.21)
        pad_inches = 0.05

        # fig.subplots_adjust(right=1.1)

    fig.savefig(output_path, bbox_inches='tight', pad_inches=pad_inches)
    plt.close()


create_in_vivo_concentrations_figure()
create_acetate_figure()
