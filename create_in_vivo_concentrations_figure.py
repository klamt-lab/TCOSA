import matplotlib.pyplot as plt
import pandas
from helper import json_load

def create_in_vivo_concentrations_figure():
    ratio_ratio_test_data_aerobic = json_load("cosa/results_aerobic/ratio_ratio_test_data.json")
    concentration = "VIVOCONC"
    target = "OPTSUBMDF"

    figurename_tuple = ("aerobic", f"2C_NADH_to_NAD___to___NADPH_to_nadp_{target}_{concentration}.jpg")

    fig, axs = plt.subplots(nrows=2, ncols=1, dpi=500, figsize=(6, 8)) #sharex=True, figsize=(50, 25), dpi=120, facecolor="white")
    fig.tight_layout(pad=3.75)

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

    table_path = f"cosa/results_aerobic/optsubmdf_table_VIVOCONC.csv"
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
            growth_rates[:11], # x
            list_to_float(table[header])[:11], # y
            label=label,
            linestyle=linestyle,
            color=color,
            linewidth=linewidth,
        )
    axs[0].legend(loc="lower left")
    axs[0].set_title("A", loc="left", fontweight="bold")
    axs[0].set_xlabel("Growth rate [1/h]")
    axs[0].set_ylabel("OptSubMDF [kJ/mol]")
    axs[0].set_xlim(min(growth_rates[:11]), max(growth_rates[:11]))


    ########################################################
    ########################################################
    ########################################################
    # 1: Ratio ratio figure
    ratio_ratio_test_data = ratio_ratio_test_data_aerobic

    min_label = "Minimal ratio"
    max_label = "Maximal ratio"

    figurename = figurename_tuple[1]
    plotted_growth_rates = ratio_ratio_test_data[figurename]["plotted_growth_rates"][:10]
    min_ratios = ratio_ratio_test_data[figurename]["min_ratios"][:10]
    max_ratios = ratio_ratio_test_data[figurename]["max_ratios"][:10]
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
    axs[1].legend(loc="upper center", ncol=2)
    axs[1].set_title("B", loc="left", fontweight="bold")
    axs[1].set_ylim(-.0005, 0.01)
    axs[1].set_xlabel("Growth rate [1/h]")
    axs[1].set_ylabel(r"$\mathrm{\frac{[NADH]/[NAD^{+}]}{[NADPH]/[NADP^{+}]}}$", fontsize=13)

    fig.savefig(f"./cosa/in_vivo_concentrations_figure.png", bbox_inches='tight', pad_inches=0.05)
    plt.close()

create_in_vivo_concentrations_figure()
