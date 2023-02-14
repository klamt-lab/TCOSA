import copy
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas
from helper import ensure_folder_existence
from statistics import mean, stdev


def create_cosa_figures(data_path: str, figures_path: str, anaerobic: bool) -> None:
    ensure_folder_existence(figures_path)

    str_to_float = lambda x: float(x.replace(",", "."))
    list_to_float = lambda x: [str_to_float(y) for y in x]

    growth_rate_id = "µ [1/h]"
    best_id = "FLEXIBLE"
    in_vivo_id = "WILDTYPE"
    only_one_id = "SINGLE_COFACTOR"

    table_paths = [
        f"{data_path}optmdf_table_STANDARDCONC.csv",
        # f"{data_path}optmdf_table_VIVOCONC.csv",
        f"{data_path}optsubmdf_table_STANDARDCONC.csv",
        # f"{data_path}optsubmdf_table_VIVOCONC.csv",
    ]
    for table_path in table_paths:
        table = pandas.read_csv(
            table_path,
            sep="\t",
        )

        title = ""
        figurename = ""
        ylabel = ""
        if "optsubmdf_" in table_path:
            title += "OptSubMDF"
            figurename += "1A"
            ylabel += "OptSubMDF"
        else:
            title += "OptMDF"
            figurename += "1B"
            ylabel += "OptMDF"
        figurename += "_"
        ylabel += " [kJ/mol]"
        if "_STANDARDCONC" in table_path:
            figurename += "STANDARDCONCS"
            titleaddition = ""
        else:
            figurename += "PAPERCONCS"
            titleaddition = "\nwith adapted measured concentration ranges"
        title += f" under {'anaerobic' if anaerobic else 'aerobic'} conditions"
        figurename += ".jpg"

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
                label = "Flexible (theoretical maximum)"
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

            plt.plot(
                growth_rates, # x
                list_to_float(table[header]), # y
                label=label,
                linestyle=linestyle,
                color=color,
                linewidth=linewidth,
            )
        plt.legend(loc="lower left")
        plt.title(title+titleaddition)
        plt.xlabel("Growth rate [1/h]")
        plt.ylabel(ylabel)
        plt.xlim(min(growth_rates), max(growth_rates))
        plt.savefig(f"{figures_path}{figurename}")
        plt.close()


def create_cosa_dG0_sampling_figures_old(
    data_path: str, figures_path: str,
    anaerobic: bool, num_samplings: int) -> None:
    ensure_folder_existence(figures_path)

    str_to_float = lambda x: float(x.replace(",", "."))
    list_to_float = lambda x: [str_to_float(y) for y in x]

    growth_rate_id = "µ [1/h]"
    best_id = "FLEXIBLE"
    in_vivo_id = "WILDTYPE"
    only_one_id = "SINGLE_COFACTOR"

    table_paths = [
        f"{data_path}optmdf_table_STANDARDCONC.csv",
        f"{data_path}optmdf_table_VIVOCONC.csv",
        f"{data_path}optsubmdf_table_STANDARDCONC.csv",
        f"{data_path}optsubmdf_table_VIVOCONC.csv",
    ]
    for current_random in range(num_samplings):
        for table_path in table_paths:
            table = pandas.read_csv(
                table_path,
                sep="\t",
            )

            title = ""
            figurename = ""
            ylabel = ""
            if "optsubmdf_" in table_path:
                title += "OptSubMDF"
                figurename += "1A"
                ylabel += "OptSubMDF"
            else:
                title += "OptMDF"
                figurename += "1B"
                ylabel += "OptMDF"
            # title += " with "
            figurename += "_"
            ylabel += " [kJ/mol]"
            if "_STANDARDCONC" in table_path:
                figurename += "STANDARDCONCS"
            else:
                figurename += "PAPERCONCS"
            title += f" under {'anaerobic' if anaerobic else 'aerobic'} conditions"
            figurename += f"_sampling_{current_random}.jpg"

            headers = list(table.keys())
            del(headers[headers.index(best_id)])
            del(headers[headers.index(only_one_id)])
            del(headers[headers.index(in_vivo_id)])
            headers += [only_one_id, in_vivo_id, best_id]

            growth_rates = list_to_float(table[growth_rate_id])
            is_first_random_wildtype = True
            is_first_random_flexible = True
            is_first_random_single_cofactor = True
            for header in headers:
                if header == growth_rate_id:
                    continue
                elif header == best_id:
                    label = "Flexible (theoretical maximum)"
                    linestyle = "-"
                    color = "yellowgreen"
                    linewidth = 1.5
                elif header == in_vivo_id:
                    label = "Wild-type specificity"
                    linestyle = "-"
                    color = "black"
                    linewidth = 1.5
                elif header == only_one_id:
                    label = "Single cofactor pool"
                    linestyle = "-"
                    color = "red"
                    linewidth = 1.5
                else:
                    if header.split("_")[-1] != str(current_random):
                        continue

                    if in_vivo_id in header:
                        if is_first_random_wildtype:
                            label = "Random set dG0 wild-type specificity"
                            is_first_random_wildtype = False
                        else:
                            label = ""
                        color = "dimgray"
                        linestyle = "--"
                        linewidth = 4.0
                    elif only_one_id in header:
                        if is_first_random_single_cofactor:
                            label = "Random set dG0 single cofactor pool"
                            is_first_random_single_cofactor = False
                        else:
                            label = ""
                        color = "orangered"
                        linestyle = "--"
                        linewidth = 4.0
                    elif best_id in header:
                        if is_first_random_flexible:
                            label = "Random set dG0 flexible pool"
                            is_first_random_flexible = False
                        else:
                            label = ""
                        color = "lawngreen"
                        linestyle = "--"
                        linewidth = 6.0


                plt.plot(
                    growth_rates, # x
                    list_to_float(table[header]), # y
                    label=label,
                    linestyle=linestyle,
                    color=color,
                    linewidth=linewidth,
                )
            plt.legend(loc="lower left")
            plt.title(title)
            plt.xlabel("Growth rate [1/h]")
            plt.ylabel(ylabel)
            plt.xlim(min(growth_rates), max(growth_rates))
            plt.savefig(f"{figures_path}{figurename}")
            plt.close()



def create_cosa_dG0_sampling_figures(
    data_path: str, figures_path: str,
    anaerobic: bool, num_samplings: int) -> None:
    ensure_folder_existence(figures_path)

    str_to_float = lambda x: float(x.replace(",", "."))
    list_to_float = lambda x: [str_to_float(y) for y in x]

    growth_rate_id = "µ [1/h]"
    best_id = "FLEXIBLE"
    in_vivo_id = "WILDTYPE"
    only_one_id = "SINGLE_COFACTOR"

    table_paths = [
        f"{data_path}optmdf_table_STANDARDCONC.csv",
        f"{data_path}optmdf_table_VIVOCONC.csv",
        f"{data_path}optsubmdf_table_STANDARDCONC.csv",
        f"{data_path}optsubmdf_table_VIVOCONC.csv",
    ]
    for table_path in table_paths:
        table = pandas.read_csv(
            table_path,
            sep="\t",
        )

        title = ""
        figurename = ""
        ylabel = ""
        if "optsubmdf_" in table_path:
            title += "OptSubMDF"
            figurename += "1A"
            ylabel += "OptSubMDF"
        else:
            title += "OptMDF"
            figurename += "1B"
            ylabel += "OptMDF"
        figurename += "_"
        ylabel += " [kJ/mol]"
        if "_STANDARDCONC" in table_path:
            figurename += "STANDARDCONCS"
        else:
            figurename += "PAPERCONCS"
        title += f" under {'anaerobic' if anaerobic else 'aerobic'} conditions"

        headers = list(table.keys())
        del(headers[headers.index(best_id)])
        del(headers[headers.index(only_one_id)])
        del(headers[headers.index(in_vivo_id)])
        headers += [only_one_id, in_vivo_id, best_id]

        growth_rates = list_to_float(table[growth_rate_id])
        get_empty_dict = lambda growth_rates: { x: [] for x in growth_rates }
        random_flexible_dict = get_empty_dict(growth_rates)
        random_single_cofactor_dict = get_empty_dict(growth_rates)
        random_wildtype_dict = get_empty_dict(growth_rates)
        for header in headers:
            if header.startswith("RANDOM_FLEXIBLE_"):
                target_dict = random_flexible_dict
            elif header.startswith("RANDOM_SINGLE_COFACTOR_"):
                target_dict = random_single_cofactor_dict
            elif header.startswith("RANDOM_WILDTYPE_"):
                target_dict = random_wildtype_dict
            else:
                continue
            values = list_to_float(table[header])
            print(values)
            index = 0
            for growth_rate in growth_rates:
                target_dict[growth_rate].append(values[index])
                index += 1
        print(random_flexible_dict)
        print("=")
        print(random_single_cofactor_dict)
        print("=")
        print(random_wildtype_dict)

        means_flexible = [mean(random_flexible_dict[x]) for x in random_flexible_dict.keys()]
        stdevs_flexible = [stdev(random_flexible_dict[x]) for x in random_flexible_dict.keys()]
        means_single_cofactor = [mean(random_single_cofactor_dict[x]) for x in random_single_cofactor_dict.keys()]
        stdevs_single_cofactor = [stdev(random_single_cofactor_dict[x]) for x in random_single_cofactor_dict.keys()]
        means_wildtype = [mean(random_wildtype_dict[x]) for x in random_wildtype_dict.keys()]
        stdevs_wildtype = [stdev(random_wildtype_dict[x]) for x in random_wildtype_dict.keys()]

        is_first_random_wildtype = True
        is_first_random_flexible = True
        is_first_random_single_cofactor = True
        for header in headers:
            print(header)
            y_data = list_to_float(table[header])
            yerr_data = []
            if header == growth_rate_id:
                continue
            elif header == best_id:
                label = "Flexible (theoretical maximum)"
                linestyle = "-"
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
                if in_vivo_id in header:
                    if is_first_random_wildtype:
                        label = "Random set dG0 wild-type specificity"
                        is_first_random_wildtype = False
                        y_data = means_wildtype
                        yerr_data = stdevs_wildtype
                    else:
                        continue
                    color = "dimgray"
                    linestyle = "--"
                    linewidth = 2.0
                elif only_one_id in header:
                    if is_first_random_single_cofactor:
                        label = "Random set dG0 single cofactor pool"
                        is_first_random_single_cofactor = False
                        y_data = means_single_cofactor
                        yerr_data = stdevs_single_cofactor
                    else:
                        continue
                    color = "orangered"
                    linestyle = "--"
                    linewidth = 2.0
                elif best_id in header:
                    if is_first_random_flexible:
                        label = "Random set dG0 flexible pool"
                        is_first_random_flexible = False
                        y_data = means_flexible
                        yerr_data = stdevs_flexible
                    else:
                        continue
                    color = "lawngreen"
                    linestyle = "--"
                    linewidth = 2.0

            plt.plot(
                growth_rates, # x
                y_data, # y
                label=label,
                linestyle=linestyle,
                color=color,
                linewidth=linewidth,
            )
            if yerr_data != []:
                plt.errorbar(
                    growth_rates, # x
                    y_data, # y
                    yerr=yerr_data,
                    linestyle="None",
                    color=color,
                    capsize=3.0,
                )
        plt.legend(loc="lower left")
        plt.title(title)
        plt.xlabel("Growth rate [1/h]")
        plt.ylabel(ylabel)
        plt.xlim(min(growth_rates), max(growth_rates))
        plt.savefig(f"{figures_path}{figurename}")
        plt.close()


def create_total_dG0_sampling_figure() -> None:
    #TTTTTTTTTTTTTTT
    str_to_float = lambda x: float(x.replace(",", "."))
    list_to_float = lambda x: [str_to_float(y) for y in x]

    growth_rate_id = "µ [1/h]"
    best_id = "FLEXIBLE"
    in_vivo_id = "WILDTYPE"
    only_one_id = "SINGLE_COFACTOR"

    for concentrations in ("STANDARDCONC", "VIVOCONC"):
        fig, axs = plt.subplots(nrows=2, ncols=2, dpi=500, figsize=(13, 10)) #sharex=True, figsize=(50, 25), dpi=120, facecolor="white")
        for target in ("OptMDF", "OptSubMDF"):
            table_path_aerobic = f"./cosa/results_aerobic/dG0_sampling_range25/{target.lower()}_table_{concentrations}.csv"
            table_path_anaerobic = f"./cosa/results_anaerobic/dG0_sampling_range25/{target.lower()}_table_{concentrations}.csv"

            # fig.tight_layout(pad=10.0)

            current_table = 0
            for aerobicity in ("aerobic", "anaerobic"):
                current_table += 1
                if aerobicity == "aerobic":
                    table_path = table_path_aerobic
                else:
                    table_path = table_path_anaerobic

                if (aerobicity == "aerobic") and (target == "OptMDF"):
                    ax_x = 0
                    ax_y = 0
                    title = "A Aerobic, OptMDF"
                elif (aerobicity == "anaerobic") and (target == "OptMDF"):
                    ax_x = 1
                    ax_y = 0
                    title = "C Anaerobic, OptMDF"
                elif (aerobicity == "aerobic") and (target == "OptSubMDF"):
                    ax_x = 0
                    ax_y = 1
                    title = "B Aerobic, OptSubMDF"
                elif (aerobicity == "anaerobic") and (target == "OptSubMDF"):
                    ax_x = 1
                    ax_y = 1
                    title = "D Anaerobic, OptSubMDF"

                ylabel = f"{target} [kJ/mol]"

                table = pandas.read_csv(
                    table_path,
                    sep="\t",
                )

                headers = list(table.keys())
                del(headers[headers.index(best_id)])
                del(headers[headers.index(only_one_id)])
                del(headers[headers.index(in_vivo_id)])
                del(headers[headers.index("RANDOM_SINGLE_COFACTOR_1")])
                del(headers[headers.index("RANDOM_WILDTYPE_1")])
                del(headers[headers.index("RANDOM_FLEXIBLE_1")])
                headers += [only_one_id, in_vivo_id, best_id]
                headers = ["RANDOM_SINGLE_COFACTOR_1", "RANDOM_WILDTYPE_1", "RANDOM_FLEXIBLE_1"] + headers

                growth_rates = list_to_float(table[growth_rate_id])
                get_empty_dict = lambda growth_rates: { x: [] for x in growth_rates }
                random_flexible_dict = get_empty_dict(growth_rates)
                random_single_cofactor_dict = get_empty_dict(growth_rates)
                random_wildtype_dict = get_empty_dict(growth_rates)
                for header in headers:
                    if header.startswith("RANDOM_FLEXIBLE_"):
                        target_dict = random_flexible_dict
                    elif header.startswith("RANDOM_SINGLE_COFACTOR_"):
                        target_dict = random_single_cofactor_dict
                    elif header.startswith("RANDOM_WILDTYPE_"):
                        target_dict = random_wildtype_dict
                    else:
                        continue
                    values = list_to_float(table[header])
                    index = 0
                    for growth_rate in growth_rates:
                        target_dict[growth_rate].append(values[index])
                        index += 1

                means_flexible = [mean(random_flexible_dict[x]) for x in random_flexible_dict.keys()]
                stdevs_flexible = [stdev(random_flexible_dict[x]) for x in random_flexible_dict.keys()]
                means_single_cofactor = [mean(random_single_cofactor_dict[x]) for x in random_single_cofactor_dict.keys()]
                stdevs_single_cofactor = [stdev(random_single_cofactor_dict[x]) for x in random_single_cofactor_dict.keys()]
                means_wildtype = [mean(random_wildtype_dict[x]) for x in random_wildtype_dict.keys()]
                stdevs_wildtype = [stdev(random_wildtype_dict[x]) for x in random_wildtype_dict.keys()]

                is_first_random_wildtype = True
                is_first_random_flexible = True
                is_first_random_single_cofactor = True
                for header in headers:
                    y_data = list_to_float(table[header])
                    yerr_data = []
                    if header == growth_rate_id:
                        continue
                    elif header == best_id:
                        label = "Flexible specificity"
                        linestyle = "-"
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
                        if in_vivo_id in header:
                            if is_first_random_wildtype:
                                label = "Random ΔG'° - wild-type specificity"
                                is_first_random_wildtype = False
                                y_data = means_wildtype
                                yerr_data = stdevs_wildtype
                            else:
                                continue
                            color = "dimgray"
                            linestyle = "--"
                            linewidth = 2.0
                        elif only_one_id in header:
                            if is_first_random_single_cofactor:
                                label = "Random ΔG'° - single cofactor pool"
                                is_first_random_single_cofactor = False
                                y_data = means_single_cofactor
                                yerr_data = stdevs_single_cofactor
                            else:
                                continue
                            color = "orangered"
                            linestyle = "--"
                            linewidth = 2.0
                        elif best_id in header:
                            if is_first_random_flexible:
                                label = "Random ΔG'° - flexible specificity"
                                is_first_random_flexible = False
                                y_data = means_flexible
                                yerr_data = stdevs_flexible
                            else:
                                continue
                            color = "lawngreen"
                            linestyle = "--"
                            linewidth = 2.0
                    if current_table > 1:
                        label = None
                    axs[ax_x, ax_y].plot(
                        growth_rates, # x
                        y_data, # y
                        label=label,
                        linestyle=linestyle,
                        color=color,
                        linewidth=linewidth,
                    )
                    if yerr_data != []:
                        axs[ax_x, ax_y].errorbar(
                            growth_rates, # x
                            y_data, # y
                            yerr=yerr_data,
                            linestyle="None",
                            color=color,
                            capsize=3.0,
                        )
                    axs[ax_x, ax_y].set_ylabel(ylabel)
                    axs[ax_x, ax_y].set_xlabel(r"Growth rate [$\mathrm{h^{-1}}$]")
                    axs[ax_x, ax_y].set_xlim(min(growth_rates)-.004, max(growth_rates)+.004) ####
                    axs[ax_x, ax_y].set_title(title, loc="left", fontweight="bold")
            # fig.legend(loc="upper center")
            # plt.legend(loc='upper center')
        legend = [
            Line2D([0], [0], color="red", linestyle="-", linewidth=2.0, label="Original ΔG'°: single cofactor pool"),
            Line2D([0], [0], color="black", linestyle="-", linewidth=2.0, label="Original ΔG'°: wild-type specificity"),
            Line2D([0], [0], color="yellowgreen", linestyle="-", linewidth=2.0, label="Original ΔG'°: flexible specificity"),
            Line2D([0], [0], color="orangered", linestyle="--", linewidth=2.0, label="Random ΔG'°: single cofactor pool"),
            Line2D([0], [0], color="dimgray", linestyle="--", linewidth=2.0, label="Random ΔG'°: wild-type specificity"),
            Line2D([0], [0], color="lawngreen", linestyle="--", linewidth=2.0, label="Random ΔG'°: flexible specificity"),
        ]
        plt.legend(handles=legend, bbox_to_anchor=(0.0, 2.5), ncol=2, loc='upper center')
        plt.savefig(f"./cosa/total_dG0_sampling_figure_{concentrations}.png", bbox_inches='tight', pad_inches=0.05)
        plt.close()


def create_total_dG0_sampling_figures() -> None:
    str_to_float = lambda x: float(x.replace(",", "."))
    list_to_float = lambda x: [str_to_float(y) for y in x]

    growth_rate_id = "µ [1/h]"
    best_id = "FLEXIBLE"
    in_vivo_id = "WILDTYPE"
    only_one_id = "SINGLE_COFACTOR"

    for concentrations in ("STANDARDCONC", "VIVOCONC"):
        for target in ("optmdf", "optsubmdf"):
            table_paths = [
                f"./cosa/results_aerobic/dG0_sampling_range25/{target}_table_{concentrations}.csv",
                f"./cosa/results_anaerobic/dG0_sampling_range25/{target}_table_{concentrations}.csv",
            ]

            fig, axs = plt.subplots(nrows=2, dpi=500, figsize=(8, 10)) #sharex=True, figsize=(50, 25), dpi=120, facecolor="white")
            # fig.tight_layout(pad=10.0)

            current_table = 0
            for table_path in table_paths:
                if current_table == 0:
                    ax_y = 0
                    title = "A Aerobic"
                elif current_table == 1:
                    ax_y = 1
                    title = "B Anaerobic"
                if target == "optsubmdf":
                    ylabel = "OptSubMDF [kJ/mol]"
                else:
                    ylabel = "OptMDF [kJ/mol]"

                current_table += 1
                table = pandas.read_csv(
                    table_path,
                    sep="\t",
                )

                headers = list(table.keys())
                del(headers[headers.index(best_id)])
                del(headers[headers.index(only_one_id)])
                del(headers[headers.index(in_vivo_id)])
                del(headers[headers.index("RANDOM_SINGLE_COFACTOR_1")])
                del(headers[headers.index("RANDOM_WILDTYPE_1")])
                del(headers[headers.index("RANDOM_FLEXIBLE_1")])
                headers += [only_one_id, in_vivo_id, best_id]
                headers = ["RANDOM_SINGLE_COFACTOR_1", "RANDOM_WILDTYPE_1", "RANDOM_FLEXIBLE_1"] + headers

                growth_rates = list_to_float(table[growth_rate_id])
                get_empty_dict = lambda growth_rates: { x: [] for x in growth_rates }
                random_flexible_dict = get_empty_dict(growth_rates)
                random_single_cofactor_dict = get_empty_dict(growth_rates)
                random_wildtype_dict = get_empty_dict(growth_rates)
                for header in headers:
                    if header.startswith("RANDOM_FLEXIBLE_"):
                        target_dict = random_flexible_dict
                    elif header.startswith("RANDOM_SINGLE_COFACTOR_"):
                        target_dict = random_single_cofactor_dict
                    elif header.startswith("RANDOM_WILDTYPE_"):
                        target_dict = random_wildtype_dict
                    else:
                        continue
                    values = list_to_float(table[header])
                    index = 0
                    for growth_rate in growth_rates:
                        target_dict[growth_rate].append(values[index])
                        index += 1

                means_flexible = [mean(random_flexible_dict[x]) for x in random_flexible_dict.keys()]
                stdevs_flexible = [stdev(random_flexible_dict[x]) for x in random_flexible_dict.keys()]
                means_single_cofactor = [mean(random_single_cofactor_dict[x]) for x in random_single_cofactor_dict.keys()]
                stdevs_single_cofactor = [stdev(random_single_cofactor_dict[x]) for x in random_single_cofactor_dict.keys()]
                means_wildtype = [mean(random_wildtype_dict[x]) for x in random_wildtype_dict.keys()]
                stdevs_wildtype = [stdev(random_wildtype_dict[x]) for x in random_wildtype_dict.keys()]

                is_first_random_wildtype = True
                is_first_random_flexible = True
                is_first_random_single_cofactor = True
                for header in headers:
                    y_data = list_to_float(table[header])
                    yerr_data = []
                    if header == growth_rate_id:
                        continue
                    elif header == best_id:
                        label = "Flexible specificity"
                        linestyle = "-"
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
                        if in_vivo_id in header:
                            if is_first_random_wildtype:
                                label = "Random ΔG'° - wild-type specificity"
                                is_first_random_wildtype = False
                                y_data = means_wildtype
                                yerr_data = stdevs_wildtype
                            else:
                                continue
                            color = "dimgray"
                            linestyle = "--"
                            linewidth = 2.0
                        elif only_one_id in header:
                            if is_first_random_single_cofactor:
                                label = "Random ΔG'° - single cofactor pool"
                                is_first_random_single_cofactor = False
                                y_data = means_single_cofactor
                                yerr_data = stdevs_single_cofactor
                            else:
                                continue
                            color = "orangered"
                            linestyle = "--"
                            linewidth = 2.0
                        elif best_id in header:
                            if is_first_random_flexible:
                                label = "Random ΔG'° - flexible specificity"
                                is_first_random_flexible = False
                                y_data = means_flexible
                                yerr_data = stdevs_flexible
                            else:
                                continue
                            color = "lawngreen"
                            linestyle = "--"
                            linewidth = 2.0
                    if current_table > 1:
                        label = None
                    axs[ax_y].plot(
                        growth_rates, # x
                        y_data, # y
                        label=label,
                        linestyle=linestyle,
                        color=color,
                        linewidth=linewidth,
                    )
                    if yerr_data != []:
                        axs[ax_y].errorbar(
                            growth_rates, # x
                            y_data, # y
                            yerr=yerr_data,
                            linestyle="None",
                            color=color,
                            capsize=3.0,
                        )
                    axs[ax_y].set_ylabel(ylabel)
                    axs[ax_y].set_xlabel(r"Growth rate [$\mathrm{h^{-1}}$]")
                    axs[ax_y].set_xlim(min(growth_rates)-.004, max(growth_rates)+.004) ####
                    axs[ax_y].set_title(title, loc="left", fontweight="bold")
            # fig.legend(loc="upper center")
            # plt.legend(loc='upper center')
            legend = [
                Line2D([0], [0], color="red", linestyle="-", linewidth=2.0, label="Original ΔG'°: Single cofactor pool"),
                Line2D([0], [0], color="black", linestyle="-", linewidth=2.0, label="Original ΔG'°: Wild-type specificity"),
                Line2D([0], [0], color="yellowgreen", linestyle="-", linewidth=2.0, label="Original ΔG'°: Flexible specificity"),
                Line2D([0], [0], color="orangered", linestyle="--", linewidth=2.0, label="Random ΔG'°: single cofactor pool"),
                Line2D([0], [0], color="dimgray", linestyle="--", linewidth=2.0, label="Random ΔG'°: wild-type specificity"),
                Line2D([0], [0], color="lawngreen", linestyle="--", linewidth=2.0, label="Random ΔG'°: flexible specificity"),
            ]
            plt.legend(handles=legend, bbox_to_anchor=(0.5, 2.5), ncol=2, loc='upper center')
            plt.savefig(f"./cosa/total_dG0_sampling_figure_{concentrations}_{target}.png", bbox_inches='tight', pad_inches=0.05)
            plt.close()


def create_total_cosa_figure() -> None:
    str_to_float = lambda x: float(x.replace(",", "."))
    list_to_float = lambda x: [str_to_float(y) for y in x]

    growth_rate_id = "µ [1/h]"
    best_id = "FLEXIBLE"
    in_vivo_id = "WILDTYPE"
    only_one_id = "SINGLE_COFACTOR"
    OPTMDF_AEROBIC = 0
    OPTSUBMDF_AEROBIC = 1
    OPTMDF_ANAEROBIC = 2
    OPTSUBMDF_ANAEROBIC = 3

    for concentration in ("STANDARDCONC", "VIVOCONC"):
        table_paths = [
            f"./cosa/results_aerobic/optmdf_table_{concentration}.csv",
            f"./cosa/results_aerobic/optsubmdf_table_{concentration}.csv",
            f"./cosa/results_anaerobic/optmdf_table_{concentration}.csv",
            f"./cosa/results_anaerobic/optsubmdf_table_{concentration}.csv",
        ]

        fig, axs = plt.subplots(nrows=2, ncols=2, dpi=500, figsize=(13, 10)) #sharex=True, figsize=(50, 25), dpi=120, facecolor="white")
        # fig.tight_layout(pad=7.0)

        current_table = 0
        for table_path in table_paths:
            first = False
            if current_table == OPTMDF_AEROBIC:
                first = True
                ax_x = 0
                ax_y = 0
                title = "A Aerobic, OptMDF"
            elif current_table == OPTSUBMDF_AEROBIC:
                ax_x = 0
                ax_y = 1
                title = "B Aerobic, OptSubMDF"
            elif current_table == OPTMDF_ANAEROBIC:
                ax_x = 1
                ax_y = 0
                title = "C Anaerobic, OptMDF"
            elif current_table == OPTSUBMDF_ANAEROBIC:
                ax_x = 1
                ax_y = 1
                title = "D Anaerobic, OptSubMDF"
            current_table += 1

            table = pandas.read_csv(
                table_path,
                sep="\t",
            )

            ylabel = ""
            if "optsubmdf_" in table_path:
                ylabel += "OptSubMDF"
            else:
                ylabel += "OptMDF"
            ylabel += " [kJ/mol]"

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

                if not first:
                    label = None

                values = [x for x in table[header]]
                temp_growth_rates = copy.deepcopy(growth_rates)
                if not "anaerobic" in table_path:
                    current_index = 0
                    for _ in temp_growth_rates:
                        if temp_growth_rates[current_index] < 0.4:
                            temp_growth_rates[current_index] = None
                            values[current_index] = None
                        current_index += 1
                temp_growth_rates = [x for x in temp_growth_rates if x is not None]
                values = [x for x in values if x is not None]

                axs[ax_x, ax_y].plot(
                    temp_growth_rates, # x
                    list_to_float(values), # y
                    label=label,
                    linestyle=linestyle,
                    color=color,
                    linewidth=linewidth,
                )
                axs[ax_x, ax_y].set_title(title, loc="left", fontweight="bold")
                axs[ax_x, ax_y].set_xlabel(r"Growth rate [$\mathrm{h^{-1}}$]")
                axs[ax_x, ax_y].set_ylabel(ylabel)
                axs[ax_x, ax_y].set_xlim(min(temp_growth_rates), max(temp_growth_rates))
                if "OptMDF" in ylabel:
                    axs[ax_x, ax_y].set_ylim(0.0, max(list_to_float(values))+.25)
                else:
                    axs[ax_x, ax_y].set_ylim(0.0, 30.0)
        fig.legend(loc="upper center")
        # plt.legend(loc="lower left")
        # plt.xlabel("Growth rate [1/h]")
        # plt.ylabel(ylabel)
        #
        plt.savefig(f"./cosa/full_random_sampling_figure_{concentration}.png", bbox_inches='tight', pad_inches=0.05)
        plt.close()


def create_total_cosa_figure_optsubmdf_only() -> None:
    str_to_float = lambda x: float(x.replace(",", "."))
    list_to_float = lambda x: [str_to_float(y) for y in x]

    growth_rate_id = "µ [1/h]"
    best_id = "FLEXIBLE"
    in_vivo_id = "WILDTYPE"
    only_one_id = "SINGLE_COFACTOR"
    OPTSUBMDF_AEROBIC = 0
    OPTSUBMDF_ANAEROBIC = 1

    for concentration in ("STANDARDCONC", "VIVOCONC"):
        table_paths = [
            f"./cosa/results_aerobic/optsubmdf_table_{concentration}.csv",
            f"./cosa/results_anaerobic/optsubmdf_table_{concentration}.csv",
        ]

        fig, axs = plt.subplots(nrows=1, dpi=500, figsize=(15, 5)) #sharex=True, figsize=(50, 25), dpi=120, facecolor="white")

        current_table = 0
        for table_path in table_paths:
            first = False
            if current_table == OPTSUBMDF_AEROBIC:
                first = True
                ax_y = 0
                title = "A Aerobic, OptSubMDF"
            elif current_table == OPTSUBMDF_ANAEROBIC:
                ax_y = 1
                title = "B Anaerobic, OptSubMDF"
            current_table += 1

            table = pandas.read_csv(
                table_path,
                sep="\t",
            )

            ylabel = ""
            if "optsubmdf_" in table_path:
                ylabel += "OptSubMDF"
            else:
                ylabel += "OptMDF"
            ylabel += " [kJ/mol]"

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

                if not first:
                    label = None

                values = [x for x in table[header]]
                temp_growth_rates = copy.deepcopy(growth_rates)
                if not "anaerobic" in table_path:
                    current_index = 0
                    for _ in temp_growth_rates:
                        if temp_growth_rates[current_index] < 0.4:
                            temp_growth_rates[current_index] = None
                            values[current_index] = None
                        current_index += 1
                temp_growth_rates = [x for x in temp_growth_rates if x is not None]
                values = [x for x in values if x is not None]

                axs[ax_y].plot(
                    temp_growth_rates, # x
                    list_to_float(values), # y
                    label=label,
                    linestyle=linestyle,
                    color=color,
                    linewidth=linewidth,
                )
                axs[ax_y].set_title(title, loc="left", fontweight="bold")
                axs[ax_y].set_xlabel(r"Growth rate [$\mathrm{h^{-1}}$]")
                axs[ax_y].set_ylabel(ylabel)
                axs[ax_y].set_xlim(min(temp_growth_rates), max(temp_growth_rates))
                axs[ax_y].set_ylim(0, 30.0)
        # fig.legend(loc="upper center", ncol=2)
        legend = [
            Line2D([0], [0], color="paleturquoise", linestyle="-", linewidth=2.0, label="Random specificities"),
            Line2D([0], [0], color="red", linestyle="-", linewidth=2.0, label="Single cofactor pool"),
            Line2D([0], [0], color="black", linestyle="-", linewidth=2.0, label="Wild-type specificity"),
            Line2D([0], [0], color="yellowgreen", linestyle="--", linewidth=2.0, label="Flexible specificity"),
        ]
        plt.legend(handles=legend, bbox_to_anchor=(-0.2, 1.152), ncol=4, loc='upper center')
        # plt.xlabel("Growth rate [1/h]")
        # plt.ylabel(ylabel)
        #
        plt.savefig(f"./cosa/full_random_sampling_figure_optsubmdf_only_{concentration}.png", bbox_inches='tight', pad_inches=0.05)
        plt.close()
