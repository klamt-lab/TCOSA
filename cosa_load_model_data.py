import cobra
from cosa_get_all_tcosa_reaction_ids import get_all_tcosa_reaction_ids
from helper import json_load, json_write, pickle_load
from optmdfpathway import STANDARD_R, STANDARD_T, get_optmdfpathway_base_problem, get_thermodynamic_bottlenecks
from optimization import perform_variable_maximization, perform_variable_minimization
from fba import perform_fba_flux_maximization, get_fba_base_problem
import copy
from cosa_get_model_with_nadx_scenario import cosa_get_model_with_nadx_scenario

MIN_OPTMDF = 0.1
LOW_DG0 = -100

def load_model_data(anaerobic: bool, expanded: bool, c_source: str="glucose"):
    biomass_reaction_id = "BIOMASS_Ec_iML1515_core_75p37M"

    #### LOAD AND SET UP MODEL ####
    print("=Loading up iML model=")
    print(">Load SBML file of COSA model")
    if not expanded:
        cobra_model: cobra.Model = cobra.io.read_sbml_model("cosa/iML1515_TCOSA.xml")
    else:
        cobra_model: cobra.Model = cobra.io.read_sbml_model("cosa/iML1515_3TCOSA_expanded.xml")

    cobra_model.solver = "cplex"

    print(">Set aerobicity and C source of model")
    cobra_model.reactions.get_by_id("EX_glc__D_e_REV").upper_bound = 0.0
    if anaerobic:
        cobra_model.reactions.get_by_id("EX_o2_e_REV").lower_bound = 0.0
        cobra_model.reactions.get_by_id("EX_o2_e_REV").upper_bound = 0.0
        if c_source == "glucose":
            cobra_model.reactions.get_by_id("EX_glc__D_e_REV").upper_bound = 20.0
        elif c_source == "acetate":
            cobra_model.reactions.get_by_id("EX_ac_e_REV").upper_bound = 20.0
    else:
        if c_source == "glucose":
            cobra_model.reactions.get_by_id("EX_glc__D_e_REV").upper_bound = 10.0
        elif c_source == "acetate":
            cobra_model.reactions.get_by_id("EX_ac_e_REV").upper_bound = 10.0
    #### END OF LOAD AND SET UP MODEL ####

    ### CALCULATION OF USED GROWTH ###
    # Used growth is: min(max(µ_Single_Cofactor), max(µ_Wildtype))
    #                 without thermodynamic constraints
    print("=CALCULATION OF USED GROWTH=")
    growthtest_cobra_model = copy.deepcopy(cobra_model)
    growthtest_cobra_model = cosa_get_model_with_nadx_scenario(
        nadx_scenario="SINGLE_COFACTOR",
        cobra_model=growthtest_cobra_model,
    )
    fba_base_problem = get_fba_base_problem(
        cobra_model=growthtest_cobra_model,
        extra_constraints=[]
    )
    fba_result = perform_fba_flux_maximization(
        base_problem=fba_base_problem,
        reaction_id=biomass_reaction_id
    )
    precise_max_growth_test_1 = -fba_result["objective_value"]
    print(">Precise max growth with SINGLE_COFACTOR:", precise_max_growth_test_1, "1/h")

    bottlenecktest_cobra_model = copy.deepcopy(cobra_model)
    bottlenecktest_cobra_model = cosa_get_model_with_nadx_scenario(
        nadx_scenario="WILDTYPE",
        cobra_model=bottlenecktest_cobra_model,
    )

    fba_base_problem = get_fba_base_problem(
        cobra_model=bottlenecktest_cobra_model,
        extra_constraints=[]
    )
    fba_result = perform_fba_flux_maximization(
        base_problem=fba_base_problem,
        reaction_id=biomass_reaction_id
    )
    precise_max_growth_test_2 = -fba_result["objective_value"]
    print(">Precise max growth with WILDTYPE:", precise_max_growth_test_2, "1/h")
    print(fba_result["values"]["DHPPDA2"])

    used_growth = min(precise_max_growth_test_1, precise_max_growth_test_2) * 0.99
    print(
        "->Used growth, i.e., 0.99*min(max(µ_Single_Cofactor), max(µ_Wildtype)):",
        used_growth,
        "1/h"
    )
    ### END OF CALCULATION OF USED GROWTH ###

    #### LOAD AND SET CONCENTRATION VALUES ####
    print("=LOAD AND SET CONCENTRATION VALUES=")
    print(">Set concentration ranges from Bennett et al., 2009 (paper concentration ranges)")
    concentration_values_paper = {
        "DEFAULT": {
            "min": 1e-6,
            "max": 0.02,
        },
        "h2o_c": {
            "min": 1.0,
            "max": 1.0,
        },
        "h2o_p": {
            "min": 1.0,
            "max": 1.0,
        },
        "h2o_e": {
            "min": 1.0,
            "max": 1.0,
        },
        "h_c": {
            "min": 1.0,
            "max": 1.0,
        },
        "h_p": {
            "min": 1.0,
            "max": 1.0,
        },
        "h_e": {
            "min": 1.0,
            "max": 1.0,
        },
    }
    paper_concentration_data = json_load(
        "resources/in_vivo_concentration_data/2009_Bennet_full_data.json"
    )
    for key in paper_concentration_data.keys():
        if key in ("nad_c", "nadh_c", "nadp_c", "nadph_c", "h_c", "h2o_c"):
            continue
        concentration_values_paper[key] = {}
        concentration_values_paper[key]["min"] = paper_concentration_data[key]["lb"] / 1
        concentration_values_paper[key]["max"] = paper_concentration_data[key]["ub"] * 1
    json_write("./resources/in_vivo_concentration_data/final_concentration_values_paper.json", concentration_values_paper)

    print(">Set standard concentration ranges")
    concentration_values_free = {
        "DEFAULT": {
            "min": 1e-6,
            "max": 0.02,
        },
        "h2o_c": {
            "min": 1.0,
            "max": 1.0,
        },
        "h2o_p": {
            "min": 1.0,
            "max": 1.0,
        },
        "h2o_e": {
            "min": 1.0,
            "max": 1.0,
        },
        "h_c": {
            "min": 1.0,
            "max": 1.0,
        },
        "h_p": {
            "min": 1.0,
            "max": 1.0,
        },
        "h_e": {
            "min": 1.0,
            "max": 1.0,
        },
    }
    #### END OF LOAD AND SET CONCENTRATION VALUES ####

    #### RATIO CONSTRAINT HANDLING ####
    print("=SET CONCENTRATION RATIOS=")
    """
    Standard ratios in OptMDFpathway paper:
    {
        "c_i": "atp_c",
        "c_j": "adp_c",
        "h_min": 3,
        "h_max": 10,
    },
    {
        "c_i": "adp_c",
        "c_j": "amp_c",
        "h_min": 0.5,
        "h_max": 2,
    },
    {
        "c_i": "nad_c",
        "c_j": "nadh_c",
        "h_min": 3,
        "h_max": 10,
    },
    {
        "c_i": "nadph_c",
        "c_j": "nadp_c",
        "h_min": 3,
        "h_max": 10,
    },
    """
    ratio_constraint_data = [
    ]
    #### END OF RATIO CONSTRAINT HANDLING ####


    #### NEW dG0 DATA HANDLING ####
    #### NEW dG0 DATA HANDLING ####
    print("=SET UP dG0 DATA=")
    print(">Load precomputed eQuilibrator dG0 values")
    dG0_values = json_load("./resources/dG0_iML1515_irreversible_cleaned.json")

    ### print(">Delete precomputed multi-compartmental dG0 values")
    ### dG0_value_ids = list(dG0_values.keys())
    ### for dG0_value_id in dG0_value_ids:
    ###     if dG0_values[dG0_value_id]["num_compartments"] > 1:
    ###         del(dG0_values[dG0_value_id])
    ###

    ###

    print(">Delete all EX_ dG0 values as they have no biological meaning")
    dG0_keys = list(dG0_values.keys())
    for key in dG0_keys:
        if key.startswith("EX_"):
            del(dG0_values[key])

    keys = list(dG0_values.keys())
    reaction_ids = [x.id for x in cobra_model.reactions]
    for key in keys:
        key1 = key + "_ORIGINAL_NAD_TCOSA"
        key2 = key + "_VARIANT_NAD_TCOSA"
        key3 = key + "_ORIGINAL_NADP_TCOSA"
        key4 = key + "_VARIANT_NADP_TCOSA"
        for key_i in [key1, key2, key3, key4]:
            if key_i in reaction_ids:
                dG0_values[key_i] = dG0_values[key]
        if expanded:
            dG0_values[key+"_NADZ_TCOSA"] = copy.deepcopy(dG0_values[key])

    print("Add dG0=0 kJ/mol for all NAD(P) reactions without computed dG0")
    zeroed_reaction_ids = []
    deleted_transporters = ""
    kept_multicompartmentals = ""
    for reaction in cobra_model.reactions:
        if (reaction.id.endswith("_TCOSA")) and (reaction.id not in dG0_values):
            educt_met_ids = [x.id for x in reaction.metabolites.keys() if reaction.metabolites[x] < 0.0]
            product_met_ids = [x.id for x in reaction.metabolites.keys() if reaction.metabolites[x] > 0.0]

            nad_is_educt = "nad_tcosa_c" in educt_met_ids
            nadh_is_educt = "nadh_tcosa_c" in educt_met_ids
            nadp_is_educt = "nadp_tcosa_c" in educt_met_ids
            nadph_is_educt = "nadph_tcosa_c" in educt_met_ids

            nad_is_product = "nad_tcosa_c" in product_met_ids
            nadh_is_product = "nadh_tcosa_c" in product_met_ids
            nadp_is_product = "nadp_tcosa_c" in product_met_ids
            nadph_is_product = "nadph_tcosa_c" in product_met_ids

            if (nad_is_educt) and (nadh_is_product):
                set_dG0 = 15.791
            if (nadh_is_educt) and (nad_is_product):
                set_dG0 = -15.791
            if (nadp_is_educt) and (nadph_is_product):
                set_dG0 = 15.332
            if (nadph_is_educt) and (nadp_is_product):
                set_dG0 = -15.332

            dG0_values[reaction.id] = {}
            dG0_values[reaction.id]["dG0"] = set_dG0
            dG0_values[reaction.id]["uncertainty"] = 0.0
            dG0_values[reaction.id]["num_compartments"] = 1
            # zeroed_reaction_ids.append(reaction.id)

        if ("BIOMASS" not in reaction.id.upper()) and (reaction.id not in dG0_values) and (not reaction.id.startswith("EX_")):
            dG0_values[reaction.id] = {}
            dG0_values[reaction.id]["dG0"] = LOW_DG0
            dG0_values[reaction.id]["uncertainty"] = 0
            dG0_values[reaction.id]["num_compartments"] = 1

        if (("transport" in reaction.name) or ("Transport" in reaction.name) or ("antiport" in reaction.name) or ("symport" in reaction.name) or ("diffusion" in reaction.name) or ("export" in reaction.name) or ("import" in reaction.name) or ("flippase" in reaction.name) or ("permease" in reaction.name) or ("ABC system" in reaction.name) or ("uptake" in reaction.name)) and (reaction.id in dG0_values):
            if dG0_values[reaction.id]["num_compartments"] > 1:
                deleted_transporters += f"{reaction.id} | {dG0_values[reaction.id]['dG0']} kJ/mol | {reaction.name} | {reaction.reaction}\n"
                dG0_values[reaction.id] = {}
                dG0_values[reaction.id]["dG0"] = LOW_DG0
                dG0_values[reaction.id]["uncertainty"] = 0
                dG0_values[reaction.id]["num_compartments"] = 2
        elif (reaction.id in dG0_values):
            if dG0_values[reaction.id]["num_compartments"] > 1:
                kept_multicompartmentals += f"{reaction.id} | {dG0_values[reaction.id]['dG0']} kJ/mol | {reaction.name} | {reaction.reaction}\n"

    with open("deleted_transporters.txt", "w") as f:
        f.write(deleted_transporters)
    with open("kept_multicompartmentals.txt", "w") as f:
        f.write(kept_multicompartmentals)

    ### print(">Delete newly set multi-compartmental dG0 values")
    ### for reaction in cobra_model.reactions:
    ###     reaction: cobra.Reaction = reaction
    ###     if (reaction.id not in dG0_values):
    ###         continue
    ###     compartments = []
    ###     for metabolite in reaction.metabolites.keys():
    ###         compartments.append(metabolite.compartment)
    ###     compartments = set(compartments)
    ###     if len(compartments) > 1:
    ###         if reaction.id.endswith("_TCOSA"):
    ###             print("-", dG0_values[reaction.id], reaction.id, reaction.reaction, compartments)
    ###             continue
    ###         print(reaction.id)
    ###         # if dG0_values[reaction.id] != LOW_DG0:
    ###         #     print("-",reaction.id, reaction.reaction, compartments)
    ###         dG0_values[reaction.id]["dG0"] = LOW_DG0
    del(dG0_values["H2tex_FWD"])
    del(dG0_values["H2tex_REV"])
    del(dG0_values["H2Otex_FWD"])
    del(dG0_values["H2Otex_REV"])
    dG0_values["Single_Cofactor_Pseudoreaction"] = {}
    dG0_values["Single_Cofactor_Pseudoreaction"]["dG0"] = LOW_DG0
    dG0_values["Single_Cofactor_Pseudoreaction"]["uncertainty"] = 0
    ### END OF NEW dG0 HANDLING
    ### END OF NEW dG0 HANDLING

    ### CALCULATION AND SETTING OF BOTTLENECKS ###
    # Calculated at used growth (see above) and
    # for both SINGLE_COFACTOR and WILDTYPE and for each concentration scenatio
    tested_nadx_scenarios = ("SINGLE_COFACTOR", "WILDTYPE")
    print("\n=CALCULATION OF BOTTLENECKS=")
    for concentration_scenario in ("PAPERCONC", "STANDARDCONC"):
        print(f"==USED CONCENTRATION SCENARIO: {concentration_scenario}==")
        if concentration_scenario == "STANDARDCONC":
            used_concentrations = concentration_values_free
        elif concentration_scenario == "PAPERCONC":
            used_concentrations = concentration_values_paper


        """
        ###TTTT
        max_dG0_changes = {}
        for nadx_scenario in tested_nadx_scenarios:
            print(f"===NADX scenario: {nadx_scenario}===")
            bottlenecktest_cobra_model = copy.deepcopy(cobra_model)
            bottlenecktest_cobra_model = cosa_get_model_with_nadx_scenario(
                nadx_scenario=nadx_scenario,
                cobra_model=bottlenecktest_cobra_model,
            )

            print("====OPTMDFPATHWAY ANALYSIS BEFORE DETERMINATION OF BOTTLENECKS====")
            optmdfpathway_base_problem = get_optmdfpathway_base_problem(
                cobra_model=bottlenecktest_cobra_model,
                dG0_values=dG0_values,
                metabolite_concentration_values=used_concentrations,
                ratio_constraint_data=ratio_constraint_data,
                R=STANDARD_R,
                T=STANDARD_T,
                extra_constraints=[],
                sub_network_ids=get_all_tcosa_reaction_ids(cobra_model),
                add_optmdf_bottleneck_analysis=False,
            )
            optmdfpathway_base_variables = optmdfpathway_base_problem.variablesDict()
            optmdfpathway_base_variables[biomass_reaction_id].bounds(
                used_growth,
                1e12
            )

            optmdfpathway_result = perform_variable_maximization(
                optmdfpathway_base_problem,
                "var_B"
            )
            print("Status:", optmdfpathway_result["status"])
            print("OptMDF (var_B):", optmdfpathway_result["values"]["var_B"], "kJ/mol")
            # END OF OPTMDFPATHWAY ANALYSIS *BEFORE* DETERMINATION OF BOTTLENECK

            print("====OPTMDF ANALYSIS FOR THE DETERMINATION OF BOTTLENECKS====")
            optmdfpathway_base_problem = get_optmdfpathway_base_problem(
                cobra_model=bottlenecktest_cobra_model,
                dG0_values=dG0_values,
                metabolite_concentration_values=used_concentrations,
                ratio_constraint_data=ratio_constraint_data,
                R=STANDARD_R,
                T=STANDARD_T,
                extra_constraints=[],
                sub_network_ids=get_all_tcosa_reaction_ids(cobra_model),
                add_optmdf_bottleneck_analysis=True,
            )
            optmdfpathway_base_variables = optmdfpathway_base_problem.variablesDict()
            optmdfpathway_base_variables[biomass_reaction_id].bounds(
                used_growth,
                1e12
            )
            optmdfpathway_base_variables["var_B"].bounds(
                MIN_OPTMDF,
                1e12
            )
            optmdfpathway_result = perform_variable_minimization(
                optmdfpathway_base_problem,
                "zb_sum_var"
            )
            print("Status:", optmdfpathway_result["status"])
            print(
                f"Σ of reaction changes to achieve OptMDF of >= {MIN_OPTMDF} kJ/mol (zb_sum):",
                optmdfpathway_result["values"]["zb_sum_var"],
                "reaction changes"
            )
            print("Reached MDF (lower bound for OptMDF):", optmdfpathway_result["values"]["var_B"], "kJ/mol")

            print(f"->LIST OF FOUND BOTTLENECK CORRECTIONS FOR {nadx_scenario}:")
            for key in optmdfpathway_result["values"].keys():
                dG0_change = optmdfpathway_result["values"][key]
                if key.startswith("zb_var") and (dG0_change > 1e-3):
                    reaction_id = key.replace('zb_var_', '')
                    text = f"{reaction_id}: {dG0_change} kJ/mol"
                    print(text)

                    if reaction_id not in max_dG0_changes.keys():
                        max_dG0_changes[reaction_id] = dG0_change
                    else:
                        max_dG0_changes[reaction_id] = max(max_dG0_changes[reaction_id], dG0_change)
        print(f"-->MAX ΔG'0 OF {tested_nadx_scenarios} WITH {concentration_scenario} AND ANAEROBICITY: {anaerobic}")
        prefix = concentration_scenario.lower()
        print(len(max_dG0_changes.keys()))
        print("Anaerobic:", anaerobic, "Concs:", concentration_scenario)
        print("...as Python code:")
        for key in max_dG0_changes.keys():
            reverted_reaction_ids = []
            added_reaction_ids = []
            newkeys = [key]
            if "_FWD" in key:
                newkeys.append(key.replace("_FWD", "_REV"))
                added_reaction_ids.append(newkeys[-1])
                reverted_reaction_ids.append(newkeys[-1])
            elif "_REV" in key:
                newkeys.append(key.replace("_REV", "_FWD"))
                added_reaction_ids.append(newkeys[-1])
                reverted_reaction_ids.append(newkeys[-1])

            if "_ORIGINAL_NAD_" in key:
                for newkey in copy.deepcopy(newkeys):
                    newkeys.append(newkey.replace("_ORIGINAL_NAD_", "_VARIANT_NADP_"))
                    added_reaction_ids.append(newkeys[-1])
                    if newkey in reverted_reaction_ids:
                        reverted_reaction_ids.append(newkeys[-1])
            elif "_ORIGINAL_NADP_" in key:
                for newkey in copy.deepcopy(newkeys):
                    newkeys.append(newkey.replace("_ORIGINAL_NADP_", "_VARIANT_NAD_"))
                    added_reaction_ids.append(newkeys[-1])
                    if newkey in reverted_reaction_ids:
                        reverted_reaction_ids.append(newkeys[-1])
            elif "_VARIANT_NAD_" in key:
                for newkey in copy.deepcopy(newkeys):
                    newkeys.append(newkey.replace("_VARIANT_NAD_", "_ORIGINAL_NADP_"))
                    added_reaction_ids.append(newkeys[-1])
                    if newkey in reverted_reaction_ids:
                        reverted_reaction_ids.append(newkeys[-1])
            elif "_VARIANT_NADP_" in key:
                for newkey in copy.deepcopy(newkeys):
                    newkeys.append(newkey.replace("_VARIANT_NADP_", "_ORIGINAL_NAD_"))
                    added_reaction_ids.append(newkeys[-1])
                    if newkey in reverted_reaction_ids:
                        reverted_reaction_ids.append(newkeys[-1])

            for newkey in newkeys:
                set_dG0 = LOW_DG0 if newkey not in reverted_reaction_ids else -LOW_DG0
                comment = "" if newkey not in added_reaction_ids else " # added for consistency"
                print(f'{prefix}_dG0_values["{newkey}"]["dG0"] = {set_dG0}{comment}')
        input("X")
        # END OF OPTMDF ANALYSIS *FOR* THE DETERMINATION OF BOTTLENECKS
        #TTTT
        """


        print(f"===TEST OF COMBINED SINGLE_COFACTOR/WILDTYPE BOTTLENECK 'REMOVALS' WITH {concentration_scenario}===")
        print("Set up dG0 fixed bottleneck 'removals' (done for consistency between runs with different solvers)")
        if c_source == "glucose":
            ###GLUCOSE###
            if concentration_scenario == "STANDARDCONC":
                standardconc_dG0_values = copy.deepcopy(dG0_values)

                if not anaerobic: # aerobic
                    standardconc_dG0_values["KDOCT2"]["dG0"] = -100
                    standardconc_dG0_values["MECDPS"]["dG0"] = -100
                    standardconc_dG0_values["DHPPDA2"]["dG0"] = -100
                    standardconc_dG0_values["ATPPRT"]["dG0"] = -100
                    standardconc_dG0_values["IG3PS"]["dG0"] = -100
                    standardconc_dG0_values["MCTP1App"]["dG0"] = -100
                    standardconc_dG0_values["MALCOAMT"]["dG0"] = -100
                    standardconc_dG0_values["AIRC3_REV"]["dG0"] = -100
                    standardconc_dG0_values["AIRC3_FWD"]["dG0"] = 100 # added for consistency
                    standardconc_dG0_values["SHCHD2_ORIGINAL_NAD_TCOSA"]["dG0"] = -100
                    standardconc_dG0_values["SHCHD2_VARIANT_NADP_TCOSA"]["dG0"] = -100 # added for consistency

                    if expanded:
                        standardconc_dG0_values["SHCHD2_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                else: # anaerobic
                    standardconc_dG0_values["KDOCT2"]["dG0"] = -100
                    standardconc_dG0_values["MECDPS"]["dG0"] = -100
                    standardconc_dG0_values["DHPPDA2"]["dG0"] = -100
                    standardconc_dG0_values["ATPPRT"]["dG0"] = -100
                    standardconc_dG0_values["IG3PS"]["dG0"] = -100
                    standardconc_dG0_values["MCTP1App"]["dG0"] = -100
                    standardconc_dG0_values["MALCOAMT"]["dG0"] = -100
                    standardconc_dG0_values["AIRC3_REV"]["dG0"] = -100
                    standardconc_dG0_values["AIRC3_FWD"]["dG0"] = 100 # added for consistency
                    standardconc_dG0_values["SHCHD2_ORIGINAL_NAD_TCOSA"]["dG0"] = -100
                    standardconc_dG0_values["SHCHD2_VARIANT_NADP_TCOSA"]["dG0"] = -100 # added for consistency

                    if expanded:
                        standardconc_dG0_values["SHCHD2_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                test_used_dG0 = copy.deepcopy(standardconc_dG0_values)
            ###GLUCOSE###
            elif concentration_scenario == "PAPERCONC":
                paperconc_dG0_values = copy.deepcopy(dG0_values)

                if not anaerobic: # aerobic
                    paperconc_dG0_values["KDOCT2"]["dG0"] = -100
                    paperconc_dG0_values["MECDPS"]["dG0"] = -100
                    paperconc_dG0_values["DHPPDA2"]["dG0"] = -100
                    paperconc_dG0_values["ATPPRT"]["dG0"] = -100
                    paperconc_dG0_values["IG3PS"]["dG0"] = -100
                    paperconc_dG0_values["MCTP1App"]["dG0"] = -100
                    paperconc_dG0_values["MALCOAMT"]["dG0"] = -100
                    paperconc_dG0_values["AIRC3_REV"]["dG0"] = -100
                    paperconc_dG0_values["AIRC3_FWD"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["SHCHD2_ORIGINAL_NAD_TCOSA"]["dG0"] = -100
                    paperconc_dG0_values["SHCHD2_VARIANT_NADP_TCOSA"]["dG0"] = -100 # added for consistency
                    paperconc_dG0_values["ASAD_REV_VARIANT_NAD_TCOSA"]["dG0"] = -100
                    paperconc_dG0_values["ASAD_FWD_VARIANT_NAD_TCOSA"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["ASAD_REV_ORIGINAL_NADP_TCOSA"]["dG0"] = -100 # added for consistency
                    paperconc_dG0_values["ASAD_FWD_ORIGINAL_NADP_TCOSA"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["GLUDy_REV_VARIANT_NAD_TCOSA"]["dG0"] = -100
                    paperconc_dG0_values["GLUDy_FWD_VARIANT_NAD_TCOSA"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["GLUDy_REV_ORIGINAL_NADP_TCOSA"]["dG0"] = -100 # added for consistency
                    paperconc_dG0_values["GLUDy_FWD_ORIGINAL_NADP_TCOSA"]["dG0"] = 100 # added for consistency

                    if expanded:
                        paperconc_dG0_values["SHCHD2_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                        paperconc_dG0_values["ASAD_FWD_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                        paperconc_dG0_values["ASAD_REV_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                        paperconc_dG0_values["GLUDy_FWD_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                        paperconc_dG0_values["GLUDy_REV_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                else: # anaerobic
                    paperconc_dG0_values["KDOCT2"]["dG0"] = -100
                    paperconc_dG0_values["MECDPS"]["dG0"] = -100
                    paperconc_dG0_values["DHPPDA2"]["dG0"] = -100
                    paperconc_dG0_values["ATPPRT"]["dG0"] = -100
                    paperconc_dG0_values["IG3PS"]["dG0"] = -100
                    paperconc_dG0_values["MCTP1App"]["dG0"] = -100
                    paperconc_dG0_values["MALCOAMT"]["dG0"] = -100
                    paperconc_dG0_values["PTAr_FWD"]["dG0"] = -100
                    paperconc_dG0_values["PTAr_REV"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["AIRC3_REV"]["dG0"] = -100
                    paperconc_dG0_values["AIRC3_FWD"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["SHCHD2_ORIGINAL_NAD_TCOSA"]["dG0"] = -100
                    paperconc_dG0_values["SHCHD2_VARIANT_NADP_TCOSA"]["dG0"] = -100 # added for consistency
                    paperconc_dG0_values["GAPD_FWD_ORIGINAL_NAD_TCOSA"]["dG0"] = -100
                    paperconc_dG0_values["GAPD_REV_ORIGINAL_NAD_TCOSA"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["GAPD_FWD_VARIANT_NADP_TCOSA"]["dG0"] = -100 # added for consistency
                    paperconc_dG0_values["GAPD_REV_VARIANT_NADP_TCOSA"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["ACACT2r_FWD"]["dG0"] = -100
                    paperconc_dG0_values["ACACT2r_REV"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["ACACT4r_FWD"]["dG0"] = -100
                    paperconc_dG0_values["ACACT4r_REV"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["ACACT6r_FWD"]["dG0"] = -100
                    paperconc_dG0_values["ACACT6r_REV"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["ACACT1r_FWD"]["dG0"] = -100
                    paperconc_dG0_values["ACACT1r_REV"]["dG0"] = 100 # added for consistency

                    if expanded:
                        paperconc_dG0_values["SHCHD2_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                        paperconc_dG0_values["ASAD_FWD_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                        paperconc_dG0_values["ASAD_REV_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                        paperconc_dG0_values["GAPD_FWD_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                        paperconc_dG0_values["GAPD_REV_NADZ_TCOSA"]["dG0"] = -100 # added for consistency
                test_used_dG0 = copy.deepcopy(paperconc_dG0_values)
        elif c_source == "acetate":
            standardconc_dG0_values = copy.deepcopy(dG0_values)

            ###ACETATE###
            if concentration_scenario == "STANDARDCONC":
                if not anaerobic: # aerobic
                    standardconc_dG0_values["KDOCT2"]["dG0"] = -100
                    standardconc_dG0_values["MECDPS"]["dG0"] = -100
                    standardconc_dG0_values["DHPPDA2"]["dG0"] = -100
                    standardconc_dG0_values["ATPPRT"]["dG0"] = -100
                    standardconc_dG0_values["IG3PS"]["dG0"] = -100
                    standardconc_dG0_values["MCTP1App"]["dG0"] = -100
                    standardconc_dG0_values["MALCOAMT"]["dG0"] = -100
                    standardconc_dG0_values["AIRC3_REV"]["dG0"] = -100
                    standardconc_dG0_values["AIRC3_FWD"]["dG0"] = 100 # added for consistency
                    standardconc_dG0_values["SHCHD2_ORIGINAL_NAD_TCOSA"]["dG0"] = -100
                    standardconc_dG0_values["SHCHD2_VARIANT_NADP_TCOSA"]["dG0"] = -100 # added for consistency
                    if expanded:
                        pass
                else:
                    pass
                    if expanded:
                        pass
            elif concentration_scenario == "PAPERCONC":
                paperconc_dG0_values = copy.deepcopy(dG0_values)

                if not anaerobic: # aerobic
                    paperconc_dG0_values["KDOCT2"]["dG0"] = -100
                    paperconc_dG0_values["MECDPS"]["dG0"] = -100
                    paperconc_dG0_values["DHPPDA2"]["dG0"] = -100
                    paperconc_dG0_values["ATPPRT"]["dG0"] = -100
                    paperconc_dG0_values["IG3PS"]["dG0"] = -100
                    paperconc_dG0_values["MCTP1App"]["dG0"] = -100
                    paperconc_dG0_values["MALCOAMT"]["dG0"] = -100
                    paperconc_dG0_values["ENO_REV"]["dG0"] = -100
                    paperconc_dG0_values["ENO_FWD"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["AIRC3_REV"]["dG0"] = -100
                    paperconc_dG0_values["AIRC3_FWD"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["SHCHD2_ORIGINAL_NAD_TCOSA"]["dG0"] = -100
                    paperconc_dG0_values["SHCHD2_VARIANT_NADP_TCOSA"]["dG0"] = -100 # added for consistency
                    paperconc_dG0_values["PGCD_ORIGINAL_NAD_TCOSA"]["dG0"] = -100
                    paperconc_dG0_values["PGCD_VARIANT_NADP_TCOSA"]["dG0"] = -100 # added for consistency
                    paperconc_dG0_values["MDH_FWD_ORIGINAL_NAD_TCOSA"]["dG0"] = -100
                    paperconc_dG0_values["MDH_REV_ORIGINAL_NAD_TCOSA"]["dG0"] = 100 # added for consistency
                    paperconc_dG0_values["MDH_FWD_VARIANT_NADP_TCOSA"]["dG0"] = -100 # added for consistency
                    paperconc_dG0_values["MDH_REV_VARIANT_NADP_TCOSA"]["dG0"] = 100 # added for consistency

                    if expanded:
                        pass
                else:
                    pass
                    if expanded:
                        pass

        """
        #FFFF
        ### START OF FINAL TEST ###
        print("=TEST OPTMDFPATHWAYs AFTER BOTTLENECK MITIGATION=")
        for nadx_scenario in tested_nadx_scenarios:
            print(f"===NADX scenario: {nadx_scenario}===")
            test_cobra_model = copy.deepcopy(cobra_model)
            test_cobra_model = cosa_get_model_with_nadx_scenario(
                nadx_scenario=nadx_scenario,
                cobra_model=test_cobra_model,
            )

            print("====OPTMDFPATHWAY ANALYSIS BEFORE DETERMINATION OF BOTTLENECKS====")
            optmdfpathway_base_problem = get_optmdfpathway_base_problem(
                cobra_model=test_cobra_model,
                dG0_values=test_used_dG0,
                metabolite_concentration_values=used_concentrations,
                ratio_constraint_data=ratio_constraint_data,
                R=STANDARD_R,
                T=STANDARD_T,
                extra_constraints=[],
                sub_network_ids=get_all_tcosa_reaction_ids(cobra_model),
                add_optmdf_bottleneck_analysis=False,
            )
            optmdfpathway_base_variables = optmdfpathway_base_problem.variablesDict()
            optmdfpathway_base_variables[biomass_reaction_id].bounds(
                used_growth,
                1e12
            )
            optmdfpathway_base_variables["var_B"].bounds(
                MIN_OPTMDF,
                1e12
            )

            optmdfpathway_result = perform_variable_maximization(
                optmdfpathway_base_problem,
                "var_B"
            )
            print("Status:", optmdfpathway_result["status"])
            print("OptMDF (var_B):", optmdfpathway_result["values"]["var_B"], "kJ/mol")
        input("X")
        ### END OF FINAL TEST ###
        #FFFF
        """


    print("Used growth [1/h]:", used_growth)

    json_write("./cosa/dG0_values.json", dG0_values)
    json_write("./cosa/standardconc_dG0_values.json", standardconc_dG0_values)
    json_write("./cosa/paperconc_dG0_values.json", paperconc_dG0_values)

    ### LEGACY VARIABLES (NOT USED ANYMORE) ###
    def get_sorted_base_ids(reaction_ids): return sorted(list(set(
        [x.replace("_NADX", "").replace("_NADY", "").replace("_NADZ", "") for x in reaction_ids]
    )))

    nad_base_ids = get_sorted_base_ids([x.id.replace("_ORIGINAL_NAD_TCOSA", "") for x in cobra_model.reactions if x.id.endswith("_ORIGINAL_NAD_TCOSA")])
    nadp_base_ids = get_sorted_base_ids([x.id.replace("_ORIGINAL_NADP_TCOSA", "") for x in cobra_model.reactions if x.id.endswith("_ORIGINAL_NADP_TCOSA")])
    all_base_ids = sorted(nad_base_ids + nadp_base_ids)
    num_nad_base_ids = len(nad_base_ids)
    num_nadp_base_ids = len(nadp_base_ids)
    num_nad_and_nadp_reactions = num_nad_base_ids + num_nadp_base_ids
    ### END OF UNUSED LEGACY VARIABLES ###

    return all_base_ids, cobra_model, concentration_values_free, concentration_values_paper,\
        standardconc_dG0_values, paperconc_dG0_values,\
        num_nad_and_nadp_reactions, num_nad_base_ids, num_nadp_base_ids,\
        ratio_constraint_data, nad_base_ids, nadp_base_ids, used_growth, zeroed_reaction_ids


if __name__ == "__main__":
    print("AEROBIC - NOT EXPANDED")
    # load_model_data(anaerobic=False, expanded=True)
    load_model_data(anaerobic=False, expanded=False, c_source="acetate")
    print("=============")
    # load_model_data(anaerobic=True, expanded=False)
    # print("~~~~~~~")
    # print("~~~~~~~")
    # print("~~~~~~~")
    # print("~~~~~~~")
    # print("~~~~~~~")
    # print("ANAEROBIC - NOT EXPANDED")
    # load_model_data(anaerobic=True, expanded=False)
    # print("======")
    # print("======")
    # print("======")
    # print("======")
    # print("AEROBIC - EXPANDED")
    # load_model_data(anaerobic=False, expanded=True)
    # print("=====")
    # print("=====")
    # print("ANAEROBIC - EXPANDED")
    # load_model_data(anaerobic=True, expanded=True)
