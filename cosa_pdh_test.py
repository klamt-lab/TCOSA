import matplotlib.pyplot as plt
import cobra
import copy
import pulp
import math
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from cosa_get_suffix import cosa_get_suffix
from helper import json_write, json_zip_load
from typing import List
from optmdfpathway import (
    STANDARD_R, STANDARD_T, get_optmdfpathway_base_problem,
    add_differential_reactions_constraints, get_z_variable_status,
)
from optimization import perform_variable_minimization, perform_variable_maximization
from cosa_load_model_data import (
    MIN_OPTMDF, load_model_data
)
from typing import Dict
from helper import ensure_folder_existence
from fva import perform_variability_analysis


def cosa_single_swap_test(anaerobic: bool, reac_id: str, mu: float, base_nadx_scenario: str) -> None:
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=False)

    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    suffix = cosa_get_suffix(anaerobic, expanded=False)

    # figures_path = f"./cosa/results{suffix}/figures/"
    # ensure_folder_existence(figures_path)

    report = ""
    original_cobra_model = copy.deepcopy(cobra_model)
    for concentrations in ("STANDARDCONCS",): #, "PAPERCONCS"):
        print(f"=CONCENTRATION RANGES: {concentrations}=")
        report += f"=CONCENTRATION RANGES: {concentrations}=\n"
        if concentrations == "STANDARDCONCS":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
            used_concentration_values = concentration_values_free
            titleaddition = ""
        elif concentrations == "PAPERCONCS":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper
            titleaddition = "\nwith adapted measured concentration ranges"

        if anaerobic:
            targets = ("OPTSUBMDF", "OPTMDF")
        else:
            targets = ("OPTMDF",)
        targets = ("OPTSUBMDF",)
        for target in targets:
            print(f"===OPTIMIZATION TARGET: {target}===")
            report += f"===OPTIMIZATION TARGET: {target}===\n"
            cobra_model = copy.deepcopy(original_cobra_model)
            cobra_model = cosa_get_model_with_nadx_scenario(
                nadx_scenario=base_nadx_scenario,
                cobra_model=cobra_model,
            )

            optmdfpathway_base_problem = get_optmdfpathway_base_problem(
                cobra_model=cobra_model,
                dG0_values=dG0_values,
                metabolite_concentration_values=used_concentration_values,
                ratio_constraint_data=ratio_constraint_data,
                R=STANDARD_R,
                T=STANDARD_T,
                extra_constraints=[],
                sub_network_ids=get_all_tcosa_reaction_ids(cobra_model),
            )
            optmdfpathway_base_variables: Dict[str, pulp.LpVariable] = optmdfpathway_base_problem.variablesDict()

            optmdfpathway_base_variables["var_B"].bounds(
                MIN_OPTMDF,
                1e12,
            )
            optmdfpathway_base_variables[biomass_reaction_id].bounds(
                mu,
                1e12,
            )

            if target == "OPTSUBMDF":
                optmdf_result = perform_variable_maximization(base_problem=optmdfpathway_base_problem, variable_id="var_B2")
                min_optsubmdf = optmdf_result["values"]["var_B2"]
                print("Unswapped OptSubMDF:", min_optsubmdf)
                optmdfpathway_base_variables["var_B"].bounds(
                    MIN_OPTMDF,
                    1e12,
                )
                optmdfpathway_base_variables["var_B2"].bounds(
                    min_optsubmdf,
                    1e12,
                )
            else:
                optmdf_result = perform_variable_maximization(base_problem=optmdfpathway_base_problem, variable_id="var_B")
                min_optmdf = optmdf_result["values"]["var_B"]
                print("Unswapped OptSubMDF:", min_optmdf)
                optmdfpathway_base_variables["var_B"].bounds(
                    min_optmdf,
                    1e12,
                )
                optmdfpathway_base_variables["var_B2"].bounds(
                    -1e12,
                    1e12,
                )

            if reac_id.endswith("_ORIGINAL_NAD_TCOSA"):
                other_id = reac_id.replace("_ORIGINAL_NAD_TCOSA", "_VARIANT_NADP_TCOSA")
            elif reac_id.endswith("_VARIANT_NAD_TCOSA"):
                other_id = reac_id.replace("_VARIANT_NAD_TCOSA", "_ORIGINAL_NADP_TCOSA")
            elif reac_id.endswith("_ORIGINAL_NADP_TCOSA"):
                other_id = reac_id.replace("_ORIGINAL_NADP_TCOSA", "_VARIANT_NAD_TCOSA")
            elif reac_id.endswith("_VARIANT_NADP_TCOSA"):
                other_id = reac_id.replace("_VARIANT_NADP_TCOSA", "_ORIGINAL_NAD_TCOSA")
            else:
                other_id = ""

            tested_vars = []
            tested_vars = [f"f_var_{x}" for x in get_all_tcosa_reaction_ids(cobra_model)] + get_all_tcosa_reaction_ids(cobra_model)
            tested_vars = [x for x in tested_vars if ((("_ORIGINAL_") in x) or (("_VARIANT_") in x)) and (x in optmdf_result["values"].keys())]
            if other_id != "":
                tested_vars += [other_id]
            tested_vars += ["x_nad_tcosa_c", "x_nadp_tcosa_c", "x_nadh_tcosa_c", "x_nadph_tcosa_c"]
            tested_vars += ["var_B2", "var_B", biomass_reaction_id]
            fva_results = {}
            fva_results["unswapped"] = perform_variability_analysis(
                tested_vars=tested_vars,
                base_problem=optmdfpathway_base_problem,
            )

            json_write(f"./cosa/results{suffix}/variability_{reac_id}_swap_{concentrations}_at_max_{target}_{base_nadx_scenario}.json", fva_results)

            if other_id == "":
                continue

            # Do swap
            optmdfpathway_base_variables[reac_id].bounds(
                0.0,
                0.0,
            )
            optmdfpathway_base_variables[other_id].bounds(
                0.0,
                1000.0,
            )
            # end of Do swap

            if target == "OPTSUBMDF":
                optmdf_result = perform_variable_maximization(base_problem=optmdfpathway_base_problem, variable_id="var_B2")
                min_optsubmdf = optmdf_result["values"]["var_B2"]
                print("Swapped OptSubMDF:", min_optsubmdf)
                optmdfpathway_base_variables["var_B"].bounds(
                    MIN_OPTMDF,
                    1e12,
                )
                optmdfpathway_base_variables["var_B2"].bounds(
                    min_optsubmdf,
                    1e12,
                )
            else:
                optmdf_result = perform_variable_maximization(base_problem=optmdfpathway_base_problem, variable_id="var_B")
                min_optmdf = optmdf_result["values"]["var_B"]
                print("Swapped OptSubMDF:", optmdf_result)
                optmdfpathway_base_variables["var_B2"].bounds(
                    -1e12,
                    1e12,
                )
                optmdfpathway_base_variables["var_B"].bounds(
                    min_optmdf,
                    1e12,
                )

            fva_results["swapped"] = perform_variability_analysis(
                tested_vars=tested_vars,
                base_problem=optmdfpathway_base_problem,
            )

            json_write(f"./cosa/results{suffix}/variability_{reac_id}_swap_{concentrations}_at_max_{target}_{base_nadx_scenario}.json", fva_results)

            # Undo swap
            optmdfpathway_base_variables[reac_id].bounds(
                0.0,
                1000.0,
            )
            optmdfpathway_base_variables[other_id].bounds(
                0.0,
                0.0,
            )
            # end of Undo swap


# cosa_single_swap_test(
#     anaerobic=False,
#     reac_id="PDH_ORIGINAL_NAD_TCOSA",
#     mu=0.868,
#     base_nadx_scenario="WILDTYPE"
# )
# cosa_single_swap_test(
#     anaerobic=True,
#     reac_id="PDH_ORIGINAL_NAD_TCOSA",
#     mu=0.371,
#     base_nadx_scenario="WILDTYPE"
# )

cosa_single_swap_test(
    anaerobic=False,
    reac_id="TEST_0_868",
    mu=0.868,
    base_nadx_scenario="WILDTYPE"
)
cosa_single_swap_test(
    anaerobic=False,
    reac_id="TEST_0_818",
    mu=0.818,
    base_nadx_scenario="WILDTPYE"
)
cosa_single_swap_test(
    anaerobic=False,
    reac_id="TEST_0_768",
    mu=0.768,
    base_nadx_scenario="WILDTPYE"
)
cosa_single_swap_test(
    anaerobic=False,
    reac_id="TEST_0_718",
    mu=0.718,
    base_nadx_scenario="WILDTPYE"
)
cosa_single_swap_test(
    anaerobic=False,
    reac_id="TEST_0_668",
    mu=0.668,
    base_nadx_scenario="WILDTPYE"
)
cosa_single_swap_test(
    anaerobic=False,
    reac_id="TEST_0_618",
    mu=0.618,
    base_nadx_scenario="WILDTPYE"
)
cosa_single_swap_test(
    anaerobic=False,
    reac_id="TEST_0_568",
    mu=0.568,
    base_nadx_scenario="WILDTPYE"
)
cosa_single_swap_test(
    anaerobic=False,
    reac_id="TEST_0_518",
    mu=0.518,
    base_nadx_scenario="WILDTPYE"
)
cosa_single_swap_test(
    anaerobic=False,
    reac_id="TEST_0_468",
    mu=0.468,
    base_nadx_scenario="WILDTPYE"
)
cosa_single_swap_test(
    anaerobic=False,
    reac_id="TEST_0_418",
    mu=0.418,
    base_nadx_scenario="WILDTPYE"
)

cosa_single_swap_test(
    anaerobic=True,
    reac_id="TEST_0_418",
    mu=0.418,
    base_nadx_scenario="WILDTPYE"
)
