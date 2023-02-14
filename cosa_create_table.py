from helper import get_files, json_zip_load


def create_cosa_tables(data_path: str, output_path: str) -> None:
    filepaths = get_files(data_path)
    filepaths = [data_path+"/"+x for x in filepaths]

    for concentration_scenario in ("STANDARDCONC", "VIVOCONC"):
        optmdf_results = {}
        optsubmdf_results = {}
        for filepath in filepaths:
            concentration_part = f"_{concentration_scenario}_"
            if concentration_part not in filepath:
                continue
            nadx_scenario = filepath.split(concentration_part)[1].split(".")[0]

            jsondata = json_zip_load(filepath.replace(".zip", ""))
            for rounded_used_growth in jsondata.keys():
                if "OPTSUBMDF_" in filepath:
                    if rounded_used_growth not in optsubmdf_results:
                        optsubmdf_results[rounded_used_growth] = {}
                    optsubmdf_results[rounded_used_growth][nadx_scenario] = str(
                        -round(jsondata[rounded_used_growth]["objective_value"], 3)
                    ).replace(".", ",")
                elif "OPTMDF_" in filepath:
                    if rounded_used_growth not in optmdf_results:
                        optmdf_results[rounded_used_growth] = {}
                    optmdf_results[rounded_used_growth][nadx_scenario] = str(
                        -round(jsondata[rounded_used_growth]["objective_value"], 3)
                    ).replace(".", ",")
        print(">Generate OPTMDF table")
        csv_string = "µ [1/h]" + "\t" + "\t".join([str(growth_rate) for growth_rate in optmdf_results[list(optmdf_results.keys())[0]].keys()]) + "\n"
        for growth_rate in optmdf_results.keys():
            print_output = f"{growth_rate}"
            for in_vivo_key in optmdf_results[growth_rate].keys():
                print_output += f"\t{optmdf_results[growth_rate][in_vivo_key]}"
            print_output += "\n"
            csv_string += print_output
        with open(f"{output_path}/optmdf_table_{concentration_scenario}.csv", "w", encoding="utf-8") as f:
            f.write(csv_string)

        print(">Generating SubMDF table")
        print(optsubmdf_results)
        print(concentration_scenario)
        csv_string = "µ [1/h]" + "\t" + "\t".join([str(growth_rate) for growth_rate in optsubmdf_results[list(optsubmdf_results.keys())[0]].keys()]) + "\n"
        for growth_rate in optsubmdf_results.keys():
            print_output = f"{growth_rate}"
            for in_vivo_key in optsubmdf_results[growth_rate].keys():
                print_output += f"\t{optsubmdf_results[growth_rate][in_vivo_key]}"
            print_output += "\n"
            csv_string += print_output
        with open(f"{output_path}/optsubmdf_table_{concentration_scenario}.csv", "w", encoding="utf-8") as f:
            f.write(csv_string)



def create_cosa_dG0_sampling_tables(data_path: str, output_path: str) -> None:
    filepaths = get_files(data_path)
    filepaths = [data_path+"/"+x for x in filepaths]

    for concentration_scenario in ("STANDARDCONC", "VIVOCONC"):
        optmdf_results = {}
        optsubmdf_results = {}
        for filepath in filepaths:
            concentration_part = f"_{concentration_scenario}_"
            if concentration_part not in filepath:
                continue
            nadx_scenario = filepath.split(concentration_part)[1].split(".")[0]

            jsondata = json_zip_load(filepath.replace(".zip", ""))
            for rounded_used_growth in jsondata.keys():
                if "OPTSUBMDF_" in filepath:
                    if rounded_used_growth not in optsubmdf_results:
                        optsubmdf_results[rounded_used_growth] = {}
                    optsubmdf_results[rounded_used_growth][nadx_scenario] = str(
                        -round(jsondata[rounded_used_growth]["objective_value"], 3)
                    ).replace(".", ",")
                elif "OPTMDF_" in filepath:
                    if rounded_used_growth not in optmdf_results:
                        optmdf_results[rounded_used_growth] = {}
                    optmdf_results[rounded_used_growth][nadx_scenario] = str(
                        -round(jsondata[rounded_used_growth]["objective_value"], 3)
                    ).replace(".", ",")
        print(">Generate OPTMDF table")
        csv_string = "µ [1/h]" + "\t" + "\t".join([str(growth_rate) for growth_rate in optmdf_results[list(optmdf_results.keys())[0]].keys()]) + "\n"
        for growth_rate in optmdf_results.keys():
            print_output = f"{growth_rate}"
            for in_vivo_key in optmdf_results[growth_rate].keys():
                print_output += f"\t{optmdf_results[growth_rate][in_vivo_key]}"
            print_output += "\n"
            csv_string += print_output
        with open(f"{output_path}/optmdf_table_{concentration_scenario}.csv", "w", encoding="utf-8") as f:
            f.write(csv_string)

        print(">Generating SubMDF table")
        print(optsubmdf_results)
        print(concentration_scenario)
        csv_string = "µ [1/h]" + "\t" + "\t".join([str(growth_rate) for growth_rate in optsubmdf_results[list(optsubmdf_results.keys())[0]].keys()]) + "\n"
        for growth_rate in optsubmdf_results.keys():
            print_output = f"{growth_rate}"
            for in_vivo_key in optsubmdf_results[growth_rate].keys():
                print_output += f"\t{optsubmdf_results[growth_rate][in_vivo_key]}"
            print_output += "\n"
            csv_string += print_output
        with open(f"{output_path}/optsubmdf_table_{concentration_scenario}.csv", "w", encoding="utf-8") as f:
            f.write(csv_string)
