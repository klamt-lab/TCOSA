from helper import json_load
from statistics import mean, median, stdev

data = json_load("resources/dG0_iML1515_irreversible_cleaned.json")
dG0s = []
for key in data.keys():
    dG0 = data[key]["dG0"]
    if data[key]["num_compartments"] > 1:
        continue
    if dG0 != 0.0:
        dG0s.append(dG0)
positive_dG0s = [x for x in dG0s if x>0]
negative_dG0s = [x for x in dG0s if x<0]

for values in (("positive dG0s", positive_dG0s), ("negative dG0s", negative_dG0s), ("all dG0s", dG0s)):
    print(f"==Statistics with all {values[0]}==")
    print("MEDIAN:", median(values[1]))
    print("MEAN:", mean(values[1]))
    print("STDEV:", stdev(values[1]))

low_dG0 = -100
negative_l_low = [x for x in negative_dG0s if x<low_dG0]
negative_eg_low = [x for x in negative_dG0s if x>=low_dG0]

print("len(negative_l_low)", len(negative_l_low))
print("len(negative_eg_low)", len(negative_eg_low))
print("Fraction", len(negative_l_low)/(len(negative_l_low)+len(negative_eg_low)))
