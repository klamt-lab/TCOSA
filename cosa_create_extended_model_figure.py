"""This script creates the 3-cofactor vs. 2-cofactor TCOSA model comparison figure.

This also includes the different set dG0 differences between NAD(P) and NAD(P)H.
"""

# IMPORTS #
# External
import matplotlib.pyplot as plt
# Internal
from helper import ensure_folder_existence


# PRIVATE FUNCTIONS #
def _get_data_list(file_path: str):
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    lines = [x.replace("\n", "").split("\t") for x in lines]
    flexible_location = lines[0].index("FLEXIBLE")
    data_dict = {}
    for line in lines[1:]:
        data_dict[float(line[0].replace(",", "."))] = float(line[flexible_location].replace(",", "."))
    return list(data_dict.values())[:-1], list(data_dict.keys())[:-1]


# LOGIC #
for concentration_ranges in ("STANDARDCONC",): #"VIVOCONC"):
    fig, axs = plt.subplots(nrows=1, ncols=2, dpi=500, figsize=(11, 4)) #sharex=True, figsize=(50, 25), dpi=120, facecolor="white")

    for run in (0, 1, 2):
        for aerobicity in ("aerobic", "anaerobic"):
            for change_state in (0, 30, -30):
                if change_state == 0:
                    extra = ""
                    legend_addition = " (all $\mathrm{ΔE'°}$ at -320 mV)"
                else:
                    extra = f"_nadz_change_{change_state}"
                    # legend_addition = f" (formation energy change: {'' if change_state<0 else '+'}{change_state} kJ/mol)"
                    r"$\mathrm{([NADH]/[NAD^{+}])/([NADPH]/[NADP^{+}])}$"
                    legend_addition = r" ($\mathrm{ΔE'°}$ at "+f"{'-165' if change_state == -30 else '-475'} mV for third cofactor)"
                file_path_optmdf = f"cosa/results_{aerobicity}/optmdf_table_{concentration_ranges}.csv"
                file_path_expanded_optmdf = f"cosa/results_{aerobicity}_expanded{extra}/optmdf_table_{concentration_ranges}.csv"
                file_path_optsubmdf = f"cosa/results_{aerobicity}/optsubmdf_table_{concentration_ranges}.csv"
                file_path_expanded_optsubmdf = f"cosa/results_{aerobicity}_expanded{extra}/optsubmdf_table_{concentration_ranges}.csv"

                data_optmdf, growth_rates = _get_data_list(file_path_optmdf)
                data_optmdf_expanded, _ = _get_data_list(file_path_expanded_optmdf)
                data_optsubmdf, growth_rates = _get_data_list(file_path_optsubmdf)
                data_optsubmdf_expanded, _ = _get_data_list(file_path_expanded_optsubmdf)

                if aerobicity == "aerobic":
                    axs_row = 0
                    title = "a  Aerobic"
                    label_optmdf_two = "MDF with 2 cofactors (all $\mathrm{ΔE'°}$ at -320 mV)"
                    label_optsubmdf_three = "SubMDF with 3 cofactors" + legend_addition
                    label_optmdf_three = "MDF with 3 cofactors" + legend_addition
                    label_optsubmdf_two = "SubMDF with 2 cofactors (all $\mathrm{ΔE'°}$ at -320 mV)"
                else:
                    axs_row = 1
                    title = "b  Anaerobic"
                    label_optmdf_two = None
                    label_optsubmdf_three = None
                    label_optmdf_three = None
                    label_optsubmdf_two = None

                if change_state == 30:
                    linestyle = ":"
                elif change_state == -30:
                    linestyle = "-."
                else:
                    linestyle = "-"

                if (change_state == 0) and (run == 0):
                    axs[axs_row].plot(
                        growth_rates, # x
                        data_optmdf, # y
                        label=label_optmdf_two,
                        linestyle="-",
                        color="salmon",
                        linewidth=2.0,
                    )
                if run == 0:
                    axs[axs_row].plot(
                        growth_rates, # x
                        data_optmdf_expanded, # y
                        label=label_optmdf_three,
                        linestyle=linestyle,
                        color="red",
                        linewidth=2.0,
                    )
                if run == 1:
                    if change_state == 0:
                        axs[axs_row].plot(
                            growth_rates, # x
                            data_optsubmdf, # y
                            label=label_optsubmdf_two,
                            linestyle="-",
                            color="deepskyblue",
                            linewidth=2.0,
                        )
                    axs[axs_row].plot(
                        growth_rates, # x
                        data_optsubmdf_expanded, # y
                        label=label_optsubmdf_three,
                        linestyle=linestyle,
                        color="blue",
                        linewidth=2.0,
                    )

                anaerobic = aerobicity == "anaerobic"
                axs[axs_row].set_title(title, loc="left", fontweight="bold")
                axs[axs_row].set_xlabel("Growth rate [1/h]")
                axs[axs_row].set_ylabel(r"(Sub)MDF [kJ/mol]")
                axs[axs_row].set_xlim(min(growth_rates), max(growth_rates))

                if axs_row == 0:
                    axs[axs_row].set_ylim(0.0, 50)
                else:
                    axs[axs_row].set_ylim(0.0, 40)

    fig.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=2)
    fig.savefig(f"cosa/expanded_comparison_{concentration_ranges}.png", bbox_inches='tight', pad_inches=0.0)
    fig.savefig(f"cosa/Figure5.pdf", bbox_inches='tight', pad_inches=0.0)
    plt.close()
