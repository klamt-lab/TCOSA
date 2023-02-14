import cobra
from fba import get_fba_base_problem, perform_fba_flux_maximization
from optmdfpathway import get_optmdfpathway_base_problem, get_thermodynamic_bottlenecks, perform_concentration_variability_analysis, perform_concentration_variability_analysis_multi, perform_optmdfpathway_mdf_maximization
from helper import json_load
import time

test_model = cobra.io.read_sbml_model("./resources/astheriscToymodelDouble_irreversible.xml")
test_model.reactions.get_by_id("EX_C_S_exchg_REV").upper_bound = 3.14
R = 8.314e-3  # kJ⋅K⁻1⋅mol⁻1 (standard value is in J⋅K⁻1⋅mol⁻1)
T = 298.15  # K

test_model.objective = "EX_C_P_exchg_FWD"
# test_model.reactions.get_by_id("EXCHG_strain1_B_c_to_B_FWD").upper_bound = float("inf")
# test_model.reactions.get_by_id("EXCHG_strain2_B_c_to_B_REV").upper_bound = float("inf")
test_model.optimize()

print("~~~ORIGINAL FBA~~~")
print(test_model.objective.value)
print("")

print("~~~OWN FBA~~~")
fba_base_problem = get_fba_base_problem(cobra_model=test_model, extra_constraints=[{}])
print(perform_fba_flux_maximization(fba_base_problem, "EX_C_P_exchg_FWD"))
print("")

print("~~~OptMDFpathway~~~")
concentration_values = {
    "DEFAULT": {
        "min": 1.0,
        "max": 10.0,
    },
    "P_exchg": {
        "min": 5.0,
        "max": 10.0,
    }
}
dG0_values = json_load("resources/dG0_astheriscToymodelDouble.json")
for key in [x for x in dG0_values.keys()]:
    dG0_values[key+"_FWD"] = dG0_values[key]
    dG0_values[key+"_REV"] = dG0_values[key]

optmdfpathway_base_problem = get_optmdfpathway_base_problem(
    cobra_model=test_model,
    extra_constraints=[
        {
            "EX_C_P_exchg_FWD": 1,
            "ub": 1,
            "lb": 1,
        }
    ],
    dG0_values=dG0_values,
    metabolite_concentration_values=concentration_values,
    ratio_constraint_data=[],
    R=R,
    T=T,
)

print("A) OptMDFpathway")
optmdfpathway_result = perform_optmdfpathway_mdf_maximization(
    optmdfpathway_base_problem=optmdfpathway_base_problem,
)
print(optmdfpathway_result)
print(optmdfpathway_result["values"]["var_B"])

print("\nB1) Concentration variability analysis - single-threaded")
variability_results = perform_concentration_variability_analysis(
    cobra_model=test_model,
    min_mdf=optmdfpathway_result["values"]["var_B"],
    optmdfpathway_base_problem=optmdfpathway_base_problem,
)
print(variability_results)
print(optmdfpathway_result["values"]["var_B"])

"""
print("\nB2) Concentration variability analysis - multi-processed")
variability_results = perform_concentration_variability_analysis_multi(
    cobra_model=test_model,
    min_mdf=optmdfpathway_result["values"]["var_B"],
    optmdfpathway_base_problem=optmdfpathway_base_problem,
)
print(variability_results)
"""

print("\nC) Check for thermodynamic bottlenecks")
bottlenecks, report = get_thermodynamic_bottlenecks(
    cobra_model=test_model,
    optmdfpathway_result=optmdfpathway_result,
    optmdfpathway_base_problem=optmdfpathway_base_problem
)
print(report)
