import pandas
from helper import ensure_folder_existence
from statistics import mean, stdev
from cosa_get_suffix import cosa_get_suffix

def create_wt_vs_random_statistics(anaerobic: bool) -> None:
    report = ""

    str_to_float = lambda x: float(x.replace(",", "."))
    list_to_float = lambda x: [str_to_float(y) for y in x]

    data_path = f"./cosa/results_{'anaerobic' if anaerobic else 'aerobic'}/"

    str_to_float = lambda x: float(x.replace(",", "."))
    list_to_float = lambda x: [str_to_float(y) for y in x]

    in_vivo_id = "WILDTYPE"
    growth_rate_id = "µ [1/h]"
    best_id = "FLEXIBLE"
    only_one_id = "SINGLE_COFACTOR"
    all_excluded_ids = (growth_rate_id, best_id, only_one_id)

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

        headers = list(table.keys())
        del(headers[headers.index(best_id)])
        del(headers[headers.index(only_one_id)])
        del(headers[headers.index(in_vivo_id)])
        headers += [only_one_id, in_vivo_id, best_id]

        growth_rates = list_to_float(table[growth_rate_id])
        random_values_lists = []
        for header in headers:
            if header in all_excluded_ids:
                continue
            elif header == in_vivo_id:
                in_vivo_values = list_to_float(table[header])
            elif "RANDOM" in header:
                random_values_lists.append(list_to_float(table[header]))

        title = "Results with "
        if "optsubmdf_" in table_path:
            title += "OptSubMDF"
        else:
            title += "OptMDF"
        title += " under "
        if "_STANDARDCONC" in table_path:
            title += "standard concentration ranges"
        else:
            title += "paper concentration ranges"

        num_random = len(random_values_lists)
        num_random_text = f"Number random values: {num_random}"

        table = "µ [1/h]\t# worse\t\t# better\t\t# equal\t#both\n"
        index = 0
        for growth_rate in growth_rates:
            in_vivo_value = in_vivo_values[index]
            num_worse = sum([x[index]<in_vivo_value for x in random_values_lists])
            num_better = sum([x[index]>in_vivo_value for x in random_values_lists])
            num_equal = sum([x[index]==in_vivo_value for x in random_values_lists])
            num_equal_and_better = num_equal + num_better

            num_worse_text = f"{num_worse} ({round((num_worse / num_random)*100, 1)} %)"
            num_better_text = f"{num_better} ({round((num_better / num_random)*100, 1)} %)"
            num_equal_text = f"{num_equal} ({round((num_equal / num_random)*100, 1)} %)"
            num_num_equal_and_better_text = f"{num_equal_and_better} ({round((num_equal_and_better / num_random)*100, 1)} %)"

            table += f"{growth_rate}\t{num_worse_text}\t{num_better_text}\t{num_equal_text}\t{num_num_equal_and_better_text}\t\n"
            index += 1
        print(title)
        print(num_random_text)
        print(table)
        print(num_random)
        report += f"{title}\n"
        report += f"{table}\n"

        suffix = cosa_get_suffix(anaerobic, False)
        with open(f"./cosa/results{suffix}/wt_vs_random_statistics.txt", "w", encoding="utf-8") as f:
            f.write(report)

print("AEROBIC")
create_wt_vs_random_statistics(anaerobic=False)
print("================")
print("ANAEROBIC")
print("================")
# create_wt_vs_random_statistics(anaerobic=True)
