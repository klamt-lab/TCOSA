import os

path = "./cosa/results_aerobic/runs/"
filenames = os.listdir(path)

filename_end = ".json.zip"
for target in ("OPTMDF", "OPTSUBMDF"):
    for concentration in ("STANDARDCONC", "VIVOCONC"):
        for random in ("RANDOMFIXED", "RANDOMS"):
            filename_start = f"{target}_{concentration}_{random}_"

            counter = 0
            for filename in filenames:
                if filename.startswith(filename_start):
                    counter += 1

            print("~~~")
            print(filename_start)
            print(counter)
            print("~~~")
