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
from fva import perform_variability_analysis, perform_fva_multi


core_map_reactions = [
    "EX_glc__D_e",
    "GLCptspp",
    "GLCt2pp",
    "HEX1",
    "G6PDH2r",
    "PGL",
    "F6PA",
    "PGI",
    "GND",
    "PFK",
    "FBP",
    "RPE",
    "RPI",
    "GLYCDx",
    "FBA",
    "TKT2",
    "TKT1",
    "TALA",
    "TPI",
    "G3PD2",
    "G3PD5",
    "GLYK",
    "EX_glyc_e",
    "PGK",
    "GAPD",
    "PGM",
    "ENO",
    "EDA",
    "EDD",
    "EX_h_e",
    "EX_h2o_e",
    "EX_co2_e",
    "PPC",
    "PPCK",
    "PYK",
    "PPS",
    "MGSA",
    "GLYOX3",
    "LDH_D",
    "FHL",
    "PFL",
    "ME2",
    "ME1",
    "PDH",
    "POX",
    "PTAr",
    "ACALD",
    "ACKr",
    "ALCD2x",
    "EX_ac_e",
    "EX_etoh_e",
    "CS",
    "EX_succ_e",
    "SUCCt2_2pp",
    "SUCDi",
    "FRD2",
    "FUM",
    "MDH",
    "MALS",
    "SUCOAS",
    "AKGDH",
    "ICL",
    "ACONTa",
    "ACONTb",
    "ICDHyr",
    "ADK1",
    "ATPM",
    "ATPS4rpp",
    "EX_o2_e",
    "CYTBO3_4pp",
    "NADH16pp",
    "NADH17pp",
    "NADTRHD",
    "THD2pp",
    "EX_lac__D_e",
    "EX_h2_e",
    "EX_for_e",
    "SUCCt2_3pp",
    "SUCCt1pp",
    "BIOMASS_Ec_iML1515_core_75p37M"
]


def cosa_single_swap_test(anaerobic : bool, reac_id: str, mu: float, base_nadx_scenario: str, c_source: str="glucose") -> None:
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=False, c_source=c_source)

    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    suffix = cosa_get_suffix(anaerobic, expanded=False, c_source=c_source)

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
        elif concentrations == "PAPERCONCS":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper

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

        tested_vars = [
            x for x in optmdfpathway_base_variables.keys()
            if x in get_all_tcosa_reaction_ids(cobra_model)
        ]
        for reaction in cobra_model.reactions:
            for core_map_reaction in core_map_reactions:
                if reaction.id.startswith(core_map_reaction):
                    tested_vars.append(reaction.id)
        tested_vars = list(set(tested_vars))

        # tested_vars = [f"f_var_{x}" for x in get_all_tcosa_reaction_ids(cobra_model)] + get_all_tcosa_reaction_ids(cobra_model)
        # tested_vars = [x for x in tested_vars if ((("_ORIGINAL_") in x) or (("_VARIANT_") in x)) and (x in optmdf_result["values"].keys())]
        # tested_vars += ["x_nad_tcosa_c", "x_nadp_tcosa_c", "x_nadh_tcosa_c", "x_nadph_tcosa_c"]
        # tested_vars += ["var_B2", "var_B", biomass_reaction_id]

        optmdfpathway_base_variables[biomass_reaction_id].bounds(
            mu,
            1e12,
        )
        optmdf_result = perform_variable_maximization(base_problem=optmdfpathway_base_problem, variable_id="var_B")

        fva_results = {}
        optmdfpathway_base_variables["var_B"].bounds(
            optmdf_result["values"]["var_B"],
            1e12,
        )
        print(optmdf_result["values"]["var_B"])
        # fva_results["OptMDF"] = perform_variability_analysis(
        #     tested_vars=tested_vars,
        #     base_problem=optmdfpathway_base_problem,
        # )
        fva_results["OptMDF"] = perform_fva_multi(
            var_ids=tested_vars,
            base_problem=optmdfpathway_base_problem,
        )

        optmdfpathway_base_variables["var_B"].bounds(
            MIN_OPTMDF,
            1e12,
        )
        optsubmdf_result = perform_variable_maximization(base_problem=optmdfpathway_base_problem, variable_id="var_B2")
        min_optsubmdf = optsubmdf_result["values"]["var_B2"]
        optmdfpathway_base_variables["var_B2"].bounds(
            min_optsubmdf,
            1e12,
        )
        fva_results["OptSubMDF"] = perform_fva_multi(
            var_ids=tested_vars,
            base_problem=optmdfpathway_base_problem,
        )

        json_write(f"./cosa/variability_{suffix}_{reac_id}_{concentrations}_{base_nadx_scenario}.json", fva_results)
