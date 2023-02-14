#!/usr/bin/env python3
"""[summary]"""

## IMPORTS ##
# External
import cobra
from equilibrator_api import ComponentContribution , Q_, Reaction
from typing import Any, Dict, List, Tuple


## CONSTANTS ##
STANDARD_INNER_TO_OUTER_COMPARTMENTS: List[str] = [
    "c",
    "p",
    "e",
]

STANDARD_IONIC_STRENGTHS: Dict[str, float] = {  # In mM
    "c": 250,  # Source: eQuilibrator standard
    "p": 250,
    "e": 250,
}

STANDARD_PHS: Dict[str, float] = {  # Unitless
    "c": 7.5,  # Source: Bionumbers ID 105980
    "p": 7.5,  # Source: Bionumbers ID 105980
    "e": 7.5,
}

STANDARD_PMGS: Dict[str, float] = {  # Unitless
    "c": 2.5,  # Source: eQuilibrator standard
    "p": 2.5,
    "e": 2.5,
}

STANDARD_POTENTIAL_DIFFERENCES: Dict[Tuple[str, str], float] = {  # In V
    ("c", "p"): 0.15,  # Source: eQuilibrator standard
    ("p", "e"): 0.15,
}

USED_IDENTIFIERS = [
    "inchi",
    "inchi_key",
    "metanetx.chemical",
    "bigg.metabolite",
    "kegg.compound",
    "chebi",
    "sabiork.compound",
    "metacyc.compound",
    "hmdb",
    "swisslipid",
    "reactome",
    "lipidmaps",
    "seed",
]


## PUBLIC FUNCTIONS ##
def get_model_dG0_values(cobra_model: cobra.Model,
        inner_to_outer_compartments: List[str],
        phs: Dict[str, float],
        pmgs: Dict[str, float],
        ionic_strengths: Dict[str, float],
        potential_differences: Dict[Tuple[str, str], float]
    ) -> Dict[str, Dict[str, float]]:
    reaction_dG0s: Dict[str, Dict[str, float]] = {}
    single_compartment_cc_reactions: List[Reaction] = []
    single_compartment_reaction_ids: List[str] = []
    cc = ComponentContribution()
    for reaction_x in cobra_model.reactions:
        reaction: cobra.Reaction = reaction_x

        stoichiometries: List[float] = []
        compartments: List[str] = []
        identifiers: List[str] = []
        identifier_keys: List[str] = []
        for metabolite_x in reaction.metabolites.keys():
            metabolite: cobra.Metabolite = metabolite_x
            stoichiometries.append(reaction.metabolites[metabolite])
            compartments.append(metabolite.compartment)
            identifier = ""
            for used_identifier in USED_IDENTIFIERS:
                if used_identifier not in metabolite.annotation.keys():
                    continue
                metabolite_identifiers = metabolite.annotation[used_identifier]
                if type(metabolite_identifiers) is list:
                    identifier_temp = metabolite_identifiers[0]
                elif type(metabolite_identifiers) is str:
                    identifier_temp = metabolite_identifiers
                if used_identifier == "inchi":
                    compound = cc.get_compound_by_inchi(identifier_temp)
                elif used_identifier == "inchi_key":
                    compound_list = cc.search_compound_by_inchi_key(identifier_temp)
                    if len(compound_list) > 0:
                        compound = compound_list[0]
                    else:
                        compound = None
                else:
                    identifier_temp = used_identifier + ":" + identifier_temp
                    compound = cc.get_compound(identifier_temp)
                if compound is not None:
                    identifier_key = used_identifier
                    identifier = identifier_temp
                    break
            if identifier == "":
                break
            identifier_keys.append(identifier_key)
            identifiers.append(identifier)

        if identifier == "":
            print("ERROR: Metabolite has no identifier of the given types!")
            continue

        # Check for three cases:
        # 1: Single-compartment reaction
        # 2: Double-compartment reaction
        # 3: Multi-compartment reaction (not possible)
        unique_reaction_compartments = list(set(compartments))
        num_compartments = len(unique_reaction_compartments)
        if num_compartments  == 1:
            # Set compartment conditions
            compartment = unique_reaction_compartments[0]
            cc.p_h = Q_(phs[compartment])
            cc.p_mg = Q_(pmgs[compartment])
            cc.ionic_strength = Q_(str(ionic_strengths[compartment])+"mM")

            # Build together reaction
            reaction_dict: Dict[Any, float] = {}
            for i in range(len(stoichiometries)):
                identifier_string = identifiers[i]
                identifier_key = identifier_keys[i]
                stoichiometry = stoichiometries[i]
                if identifier_key == "inchi":
                    compound = cc.get_compound_by_inchi(identifier_string)
                elif identifier_key == "inchi_key":
                    compound = cc.search_compound_by_inchi_key(identifier_string)[0]
                else:
                    compound = cc.get_compound(identifier_string)
                reaction_dict[compound] = stoichiometry
            cc_reaction = Reaction(reaction_dict)

            # Check whether or not the reaction is balanced and...
            if not cc_reaction.is_balanced():
                print(f"ERROR: Reaction {reaction.id} is not balanced")
                continue
            print(f"No error with reaction {reaction.id}, ΔG'0 succesfully calculated!")
            # ...if it is balanced, calculate dG0 later on
            single_compartment_cc_reactions.append(cc_reaction)
            single_compartment_reaction_ids.append(reaction.id)
        elif num_compartments  == 2:
            index_zero = inner_to_outer_compartments.index(unique_reaction_compartments[0])
            index_one = inner_to_outer_compartments.index(unique_reaction_compartments[1])

            if index_one > index_zero:
                outer_compartment = unique_reaction_compartments[1]
                inner_compartment = unique_reaction_compartments[0]
            else:
                outer_compartment = unique_reaction_compartments[0]
                inner_compartment = unique_reaction_compartments[1]

            ph_inner = Q_(phs[inner_compartment])
            ph_outer = Q_(phs[outer_compartment])
            ionic_strength_inner = Q_(str(ionic_strengths[inner_compartment])+" mM")
            ionic_strength_outer = Q_(str(ionic_strengths[outer_compartment])+" mM")
            pmg_inner = Q_(pmgs[inner_compartment])
            pmg_outer = Q_(pmgs[outer_compartment])

            if (inner_compartment, outer_compartment) in potential_differences.keys():
                potential_difference = Q_(str(potential_differences[(inner_compartment, outer_compartment)])+" V")
            elif (outer_compartment, inner_compartment) in potential_differences.keys():
                potential_difference = Q_(str(potential_differences[(outer_compartment, inner_compartment)])+" V")
            else:
                print("ERROR")
                continue

            inner_reaction_dict: Dict[str, float] = {}
            outer_reaction_dict: Dict[str, float] = {}
            for i in range(len(stoichiometries)):
                key = identifiers[i]
                stoichiometry = stoichiometries[i]
                try:
                    compound_key = cc.get_compound(key)
                except Exception:  # sqlalchemy.orm.exc.MultipleResultsFound
                    print("ERROR")
                    continue

                if compound_key is None:
                    print("NONE in compound")
                    continue

                if compartments[i] == inner_compartment:
                    inner_reaction_dict[compound_key] = stoichiometry
                else:
                    outer_reaction_dict[compound_key] = stoichiometry

            cc_inner_reaction = Reaction(inner_reaction_dict)
            cc_outer_reaction = Reaction(outer_reaction_dict)

            cc.p_h = ph_inner
            cc.ionic_strength = ionic_strength_inner
            cc.p_mg = pmg_inner
            try:
                standard_dg_prime = cc.multicompartmental_standard_dg_prime(
                    cc_inner_reaction,
                    cc_outer_reaction,
                    e_potential_difference=potential_difference,
                    p_h_outer=ph_outer,
                    p_mg_outer=pmg_outer,
                    ionic_strength_outer=ionic_strength_outer,
                )
                uncertainty = standard_dg_prime.error.m_as("kJ/mol")
                if uncertainty < 1000:
                    dG0 = standard_dg_prime.value.m_as("kJ/mol")
                    reaction_dG0s[reaction.id] = {}
                    reaction_dG0s[reaction.id]["dG0"] = dG0
                    reaction_dG0s[reaction.id]["uncertainty"] = abs(uncertainty)
                    reaction_dG0s[reaction.id]["num_compartments"] = 2
            except ValueError:
                print("ERROR: Multi-compartmental reaction is not balanced")
                continue
        else:
            # TODO: More extended errors
            print("ERROR: More than two compartments are not possible")
            continue


    # Calculate - as batch - the rest of the single-compartment reactions
    print("Start batch ΔG'0 calculation...")
    (
        standard_dg_primes, dg_uncertainties
    ) = cc.standard_dg_prime_multi(single_compartment_cc_reactions, uncertainty_representation="sqrt")
    print("Done!")

    for single_compartment_index in range(len(standard_dg_primes)):
        uncertainty = dg_uncertainties[single_compartment_index][0].m_as("kJ/mol")
        if True:
            dG0 = standard_dg_primes[single_compartment_index].m_as("kJ/mol")
            reaction_id = single_compartment_reaction_ids[single_compartment_index]
            reaction_dG0s[reaction_id] = {}
            reaction_dG0s[reaction_id]["dG0"] = dG0
            reaction_dG0s[reaction_id]["uncertainty"] = uncertainty
            reaction_dG0s[reaction_id]["num_compartments"] = 1
    return reaction_dG0s
