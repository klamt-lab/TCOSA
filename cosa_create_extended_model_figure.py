import matplotlib.pyplot as plt
from helper import ensure_folder_existence


def _get_data_list(file_path: str):
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    lines = [x.replace("\n", "").split("\t") for x in lines]
    flexible_location = lines[0].index("FLEXIBLE")
    data_dict = {}
    for line in lines[1:]:
        data_dict[float(line[0].replace(",", "."))] = float(line[flexible_location].replace(",", "."))
    return list(data_dict.values()), list(data_dict.keys())


for concentration_ranges in ("STANDARDCONC", "VIVOCONC"):
    fig, axs = plt.subplots(nrows=1, ncols=2, dpi=500, figsize=(11, 4)) #sharex=True, figsize=(50, 25), dpi=120, facecolor="white")
    # fig.tight_layout(pad=7.0)

    for aerobicity in ("aerobic", "anaerobic"):
        file_path_optmdf = f"cosa/results_{aerobicity}/optmdf_table_{concentration_ranges}.csv"
        file_path_expanded_optmdf = f"cosa/results_{aerobicity}_expanded/optmdf_table_{concentration_ranges}.csv"
        file_path_optsubmdf = f"cosa/results_{aerobicity}/optsubmdf_table_{concentration_ranges}.csv"
        file_path_expanded_optsubmdf = f"cosa/results_{aerobicity}_expanded/optsubmdf_table_{concentration_ranges}.csv"

        data_optmdf, growth_rates = _get_data_list(file_path_optmdf)
        data_optmdf_expanded, _ = _get_data_list(file_path_expanded_optmdf)
        data_optsubmdf, growth_rates = _get_data_list(file_path_optsubmdf)
        data_optsubmdf_expanded, _ = _get_data_list(file_path_expanded_optsubmdf)

        if aerobicity == "aerobic":
            axs_row = 0
            title = "A Aerobic"
            label_optmdf_two = "OptMDF of two-cofactor model"
            label_optsubmdf_three = "OptSubMDF of three-cofactor model"
            label_optmdf_three = "OptMDF of three-cofactor model"
            label_optsubmdf_two = "OptSubMDF of two-cofactor model"
        else:
            axs_row = 1
            title = "B Anaerobic"
            label_optmdf_two = None
            label_optsubmdf_three = None
            label_optmdf_three = None
            label_optsubmdf_two = None

        axs[axs_row].plot(
            growth_rates, # x
            data_optmdf, # y
            label=label_optmdf_two,
            linestyle="--",
            color="salmon",
            linewidth=2.0,
        )
        axs[axs_row].plot(
            growth_rates, # x
            data_optmdf_expanded, # y
            label=label_optmdf_three,
            linestyle="-",
            color="red",
            linewidth=2.0,
        )
        axs[axs_row].plot(
            growth_rates, # x
            data_optsubmdf, # y
            label=label_optsubmdf_two,
            linestyle="--",
            color="deepskyblue",
            linewidth=2.0,
        )
        axs[axs_row].plot(
            growth_rates, # x
            data_optsubmdf_expanded, # y
            label=label_optsubmdf_three,
            linestyle="-",
            color="blue",
            linewidth=2.0,
        )

        anaerobic = aerobicity == "anaerobic"
        axs[axs_row].set_title(title, loc="left", fontweight="bold")
        axs[axs_row].set_xlabel("Growth rate [1/h]")
        axs[axs_row].set_ylabel(r"MDF [kJ/mol]", fontsize=13)
        axs[axs_row].set_xlim(min(growth_rates), max(growth_rates))
        axs[axs_row].set_ylim(0.0, 30)

        # plt.legend(loc="center left")
        # plt.title(f"Comparison of two-cofactor and expanded three-cofactor models\nunder {'anaerobic' if anaerobic else 'aerobic'} conditions")
        # plt.xlabel("Growth rate [1/h]")
        # plt.ylabel("OptMDF or OptSubMDF [kJ/mol]")
        # plt.xlim(min(growth_rates), max(growth_rates))
        # plt.savefig(f"cosa/results_{aerobicity}/expanded_comparison_{aerobicity}_{concentration_ranges}.jpg")
        # plt.close()
    fig.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18))
    fig.savefig(f"cosa/expanded_comparison_{concentration_ranges}.png", bbox_inches='tight', pad_inches=0.05)
    plt.close()
