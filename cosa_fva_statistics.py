from helper import json_load, get_files

can_be_active = {
    "aerobic": {
        "OptMDF": [],
        "OptSubMDF": [],
    },
    "anaerobic": {
        "OptMDF": [],
        "OptSubMDF": [],
    }
}
for file in get_files("cosa"):
    file = "./cosa/"+file
    if not file.startswith("./cosa/variability__"):
        continue
    all_data = json_load(file)
    if "_aerobic_" in file:
        aerobicity = "aerobic"
    else:
        aerobicity = "anaerobic"
    for target in ("OptMDF", "OptSubMDF"):
        data = all_data[target]
        for element in data:
            name = element[0]
            lb = element[1] if element[1] != float("NaN") else 0
            ub = element[2] if element[2] != float("NaN") else 0

            if lb < 1e-9:
                lb = 0
            if ub < 1e-9:
                ub = 0

            if (lb <= 0.0) and (ub <= 0.0):
                continue

            base_id = name.replace("_ORIGINAL_NADP_TCOSA", "")\
                    .replace("_ORIGINAL_NAD_TCOSA", "")\
                    .replace("_VARIANT_NAD_TCOSA", "")\
                    .replace("_VARIANT_NADP_TCOSA", "")\
                    .replace("_FWD", "")\
                    .replace("_REV", "")
            if base_id not in can_be_active[aerobicity][target]:
                can_be_active[aerobicity][target].append(base_id)

all_active = []
for key1 in can_be_active.keys():
    for key2 in can_be_active[key1].keys():
        print(key1, key2, ":")
        print(len(can_be_active[key1][key2]))
        all_active += can_be_active[key1][key2]
print(all_active)
all_active = list(set(all_active))
print("ALL:", len(all_active))
print("NADK" in all_active)
print("NADTRHD" in all_active)
print("NADPPPS" in all_active)
print("THD2pp" in all_active)
