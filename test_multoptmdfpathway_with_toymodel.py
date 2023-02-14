import cobra
import copy
from optmdfpathway import get_optmdfpathway_base_problem, get_thermodynamic_bottlenecks, perform_concentration_variability_analysis, perform_concentration_variability_analysis_multi, perform_optmdfpathway_mdf_maximization
from helper import json_load, json_write
from multoptmdfpathway import get_multoptmdfpathway_base_problem
from optimization import perform_variable_maximization, perform_variable_minimization
from fba import add_extra_constraints

test_model = cobra.io.read_sbml_model("./resources/astheriscToymodelDouble_irreversible.xml")
test_model.reactions.get_by_id("EX_C_S_exchg_REV").upper_bound = 3.14
R = 8.314e-3  # kJ⋅K⁻1⋅mol⁻1 (standard value is in J⋅K⁻1⋅mol⁻1)
T = 298.15  # K

test_model_without_community = copy.deepcopy(test_model)

test_model_with_community = copy.deepcopy(test_model)
test_model_with_community.reactions.get_by_id("EXCHG_strain1_B_c_to_B_FWD").upper_bound = float("inf")
test_model_with_community.reactions.get_by_id("EXCHG_strain2_B_c_to_B_REV").upper_bound = float("inf")


print("~~~OptMDFpathway preparation~~~")
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

optmdfpathway_base_problem_with_community = get_optmdfpathway_base_problem(
    cobra_model=test_model_with_community,
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

optmdfpathway_base_problem_without_community = get_optmdfpathway_base_problem(
    cobra_model=test_model_without_community,
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
    optmdfpathway_base_problem=optmdfpathway_base_problem_with_community,
)
print("OptMDF with community:", optmdfpathway_result["values"]["var_B"], "kJ/mol")
optmdfpathway_result = perform_optmdfpathway_mdf_maximization(
    optmdfpathway_base_problem=optmdfpathway_base_problem_without_community,
)
print("OptMDF without community:", optmdfpathway_result["values"]["var_B"], "kJ/mol")

print("B) MultOptMDFpathway")
multoptmdfpathway_base_problem = get_multoptmdfpathway_base_problem(
    cobra_model=test_model_with_community,
    dG0_values=dG0_values,
    metabolite_concentration_values=concentration_values,
    ratio_constraint_data=[],
    R=R,
    T=T,
    sub_network_ids=[],
    conditions=[
        [("EX_C_P_exchg_FWD", 1, 1), ("EXCHG_strain1_B_c_to_B_FWD", 0, float("inf")), ("EXCHG_strain1_B_c_to_B_FWD", 0, float("inf"))],  # With community
        [("EX_C_P_exchg_FWD", 1, 1), ("EXCHG_strain1_B_c_to_B_FWD", 0, 0), ("EXCHG_strain1_B_c_to_B_FWD", 0, 0)],  # Without community
    ],
    condition_names=["FIRST_CONDITION", "SECOND_CONDITION"],
)
mult_variables = multoptmdfpathway_base_problem.variablesDict()
varnames = list(mult_variables.keys())
x = [x for x in varnames if varnames.count(x) > 1]
optmdfpathway_result = perform_variable_maximization(
    base_problem=multoptmdfpathway_base_problem,
    variable_id="var_B_FIRST_CONDITION",
)
json_write("x.json", optmdfpathway_result)
print("OptMDF COND0:", optmdfpathway_result["values"]["var_B_FIRST_CONDITION"], "kJ/mol")

optmdfpathway_result = perform_variable_maximization(
    base_problem=multoptmdfpathway_base_problem,
    variable_id="var_B_SECOND_CONDITION",
)
print("OptMDF COND1:", optmdfpathway_result["values"]["var_B_SECOND_CONDITION"], "kJ/mol")
