from tkinter.messagebox import NO
import matplotlib.pyplot as plt
import cobra
import copy
import os
import pulp
from math import exp
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario
from cosa_get_suffix import cosa_get_suffix
from helper import json_load, json_write, json_zip_load
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


def cosa_cva(metabolites: List[str], anaerobic: bool, expanded: bool, growth_epsilon: float = 0.01) -> None:
    suffix = cosa_get_suffix(anaerobic, expanded)
    figures_path = f"./cosa/results{suffix}/figures/"
    ensure_folder_existence(figures_path)
    all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
    standardconc_dG0_values, paperconc_dG0_values,\
    num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
    ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids = load_model_data(anaerobic=anaerobic, expanded=expanded)

    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    report = ""
    original_cobra_model = copy.deepcopy(cobra_model)
    for concentrations in ("STANDARDCONC", "VIVOCONC"):
        print(f"=CONCENTRATION RANGES: {concentrations}=")
        report += f"=CONCENTRATION RANGES: {concentrations}=\n"
        if concentrations == "STANDARDCONC":
            dG0_values = copy.deepcopy(standardconc_dG0_values)
            used_concentration_values = concentration_values_free
        elif concentrations == "VIVOCONC":
            dG0_values = copy.deepcopy(paperconc_dG0_values)
            used_concentration_values = concentration_values_paper

        for target in ("OPTMDF", "OPTSUBMDF"):
            cva_filepath = f"./cosa/results{suffix}/cva_{target}_{concentrations}.json"

            if not os.path.exists(cva_filepath):
                cva_data = {}
            else:
                cva_data = json_load(cva_filepath)

            print(f"===OPTIMIZATION TARGET: {target}===")
            report += f"===OPTIMIZATION TARGET: {target}===\n"
            cobra_model = copy.deepcopy(original_cobra_model)
            cobra_model = cosa_get_model_with_nadx_scenario(
                nadx_scenario="WILDTYPE",
                cobra_model=cobra_model,
            )

            jsondata_invivo = json_zip_load(f"cosa/results{suffix}/runs/{target}_{concentrations}_WILDTYPE.json")

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

            for metabolite in metabolites:
                # if type(metabolite) is tuple:
                #     if metabolite[0] == "SUM":
                #         for
                #     elif metabolit[0] == "RATIO":
                #         pass

                metabolite_var_id = f"x_{metabolite}"

                if metabolite_var_id in cva_data.keys():
                    continue

                growth_rates = jsondata_invivo.keys()

                for growth_rate in growth_rates:
                    if growth_rate == "0,03":
                        continue

                    growth_rate_float = float(growth_rate.replace(",", "."))
                    optmdfpathway_base_variables[biomass_reaction_id].bounds(
                        growth_rate_float-growth_epsilon,
                        1e12,
                    )

                    if target == "OPTMDF":
                        min_target = jsondata_invivo[growth_rate]["values"]["var_B"]
                        optmdfpathway_base_variables["var_B"].bounds(
                            min_target,
                            1e12,
                        )
                    elif target == "OPTSUBMDF":
                        min_target = jsondata_invivo[growth_rate]["values"]["var_B2"]
                        optmdfpathway_base_variables["var_B"].bounds(
                            MIN_OPTMDF,
                            1e12,
                        )
                        optmdfpathway_base_variables["var_B2"].bounds(
                            min_target,
                            1e12,
                        )

                    print(f" @ µ [1/h] of {growth_rate_float} and min {target} of {min_target} kJ/mol")
                    report += f" @ µ [1/h] of {growth_rate_float} and min {target} of {min_target} kJ/mol\n"


                    minimization_result = perform_variable_minimization(
                        optmdfpathway_base_problem,
                        metabolite_var_id,
                    )
                    maximization_result = perform_variable_maximization(
                        optmdfpathway_base_problem,
                        metabolite_var_id,
                    )
                    min_conc = exp(minimization_result["values"][metabolite_var_id])
                    max_conc = exp(maximization_result["values"][metabolite_var_id])

                    if metabolite_var_id not in cva_data.keys():
                        cva_data[metabolite_var_id] = {}

                    cva_data[metabolite_var_id][growth_rate] = {
                        "min": min_conc,
                        "max": max_conc,
                    }
                    json_write(cva_filepath, cva_data)

"""
# "Significant" metabolites
metabolites = {
    "nad_tcosa_c",
    "nadh_tcosa_c",
    "nadp_tcosa_c",
    "nadph_tcosa_c",
    "gthrd_c",
    "utp_c",
    "gtp_c",
    "datp_c",
    "itp_c",
    "dttp_c",
    "ctp_c",
    "dctp_c",
    "glu__L_c",
    "gln__L_c",
}
"""
metabolites = [
    "akg_c",
    "3mob_c",
    # "5oxpro_c", # Does not occur in iML1515
    "prpp_c",
    "6pgc_c",
    "accoa_c",
    "adp_c",
    "r5p_c",
    "amp_c",
    "atp_c",
    "cdp_c",
    "coa_c",
    "ctp_c",
    "datp_c",
    "dcdp_c",
    "dctp_c",
    "fdp_c",
    "f1p_c",
    "f6p_c",
    "dgdp_c",
    "gam6p_c",
    "g6p_c",
    "ru5p__D_c",
    "dtdpglu_c",
    "dttp_c",
    "dump_c",
    "fad_c",
    "gdp_c",
    "gtp_c",
    "imp_c",
    "itp_c",
    "asp__L_c",
    "citr__L_c",
    "glu__L_c",
    "gln__L_c",
    "phe__L_c",
    "ser__L_c",
    "thr__L_c",
    "trp__L_c",
    "nadh_tcosa_c",
    "nad_tcosa_c",
    "nadph_tcosa_c",
    "nadp_tcosa_c",
    "gthox_c",
    "pep_c",
    "gthrd_c",
    "udp_c",
    "ump_c",
    "utp_c",
    # ("RATIO", "nad_tcosa_c", "nadh_tcosa_c"),
    # ("RATIO", "nadp_tcosa_c", "nadph_tcosa_c"),
    # ("SUM", "2pg_c", "3pg_c"),
]

in_vivo_concentrations = json_load("resources/in_vivo_concentration_data/final_concentration_values_paper.json")

metabolites += list(in_vivo_concentrations.keys())
metabolites = list(set(metabolites))
metabolites = [x for x in metabolites if metabolites != "DEFAULT"]

cosa_cva(metabolites=metabolites, anaerobic=False, expanded=False)
cosa_cva(metabolites=metabolites, anaerobic=True, expanded=False)
