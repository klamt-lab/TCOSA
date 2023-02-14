import cobra

cobra_model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")

active_exs = []
closed_exs = []
for reaction in cobra_model.reactions:
    if not reaction.id.startswith("EX_"):
        continue

    lb = reaction.lower_bound
    ub = reaction.upper_bound
    text = f"{reaction.id}, ({lb},{ub})"
    if (lb < -0.0) or (ub > 0.0):
        active_exs.append(text)
    else:
        closed_exs.append(text)

for closed_ex in closed_exs:
    print(closed_ex)
print("====")
for active_ex in active_exs:
    print(active_ex)