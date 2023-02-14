import cobra
from typing import Dict, List


def cosa_get_model_with_nadx_scenario(
        nadx_scenario: str, cobra_model: cobra.Model,
        randoms_random_base_lists: Dict[str, Dict[str, float]]={},
        randomfixed_random_base_lists: Dict[str, Dict[str, float]]={},
        nad_base_ids: List[str]=[],
        nadp_base_ids: List[str]=[],
    ) -> cobra.Model:
    if nadx_scenario == "WILDTYPE":
        for reaction in cobra_model.reactions:
            reaction: cobra.Reaction = reaction
            key = reaction.id
            if (("_VARIANT_" in key) or ("_NADZ_" in key)) and (key.endswith("_TCOSA")):
                reaction.lower_bound = 0.0
                reaction.upper_bound = 0.0
    elif nadx_scenario == "SINGLE_COFACTOR":
        for reaction in cobra_model.reactions:
            # if reaction.id == "NADK":
            #     continue
            reaction: cobra.Reaction = reaction
            met_ids = [x.id for x in reaction.metabolites.keys()]
            if (("nadp_tcosa_c" in met_ids) or ("nadph_tcosa_c" in met_ids) or ("nadz_tcosa_c" in met_ids) or ("nadzh_tcosa_c" in met_ids)) and (not reaction.id.startswith("BIOMASS_")):
                reaction.lower_bound = 0.0
                reaction.upper_bound = 0.0
    elif not (nadx_scenario == "FLEXIBLE"):
        base_list_i = nadx_scenario.split("_")[1]
        if nadx_scenario.startswith("RANDOMS_"):
            base_list = randoms_random_base_lists[base_list_i]
        else:
            base_list = randomfixed_random_base_lists[base_list_i]
        reac_i = 0
        for base_id_key in base_list.keys():
            if base_list[base_id_key] <= 0.0:
                if base_id_key in nad_base_ids:
                    id_addition = "_ORIGINAL_"
                elif base_id_key in nadp_base_ids:
                    id_addition = "_VARIANT_"
                else:
                    print(base_id_key)
                    raise(AttributeError)
                cobra_model.reactions.get_by_id(
                    f"{base_id_key}{id_addition}NAD_TCOSA").lower_bound = 0.0
                cobra_model.reactions.get_by_id(
                    f"{base_id_key}{id_addition}NAD_TCOSA").upper_bound = 0.0
            else:
                if base_id_key in nadp_base_ids:
                    id_addition = "_ORIGINAL_"
                elif base_id_key in nad_base_ids:
                    id_addition = "_VARIANT_"
                else:
                    print(base_id_key)
                    raise(AttributeError)
                cobra_model.reactions.get_by_id(
                    f"{base_id_key}{id_addition}NADP_TCOSA").lower_bound = 0.0
                cobra_model.reactions.get_by_id(
                    f"{base_id_key}{id_addition}NADP_TCOSA").upper_bound = 0.0
            reac_i += 1
    if nadx_scenario == "SINGLE_COFACTOR":
        cobra_model.reactions.get_by_id("Single_Cofactor_Pseudoreaction").lower_bound = 0.0
        cobra_model.reactions.get_by_id("Single_Cofactor_Pseudoreaction").upper_bound = 1000.0
    else:
        cobra_model.reactions.get_by_id("Single_Cofactor_Pseudoreaction").lower_bound = 0.0
        cobra_model.reactions.get_by_id("Single_Cofactor_Pseudoreaction").upper_bound = 0.0
    return cobra_model
