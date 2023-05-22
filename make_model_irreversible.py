#!/usr/bin/env python3
"""Includes a function with which reversible reactions of a cobrapy model can be splitted."""

## IMPORTS ##
# External
import cobra
import copy


## PUBLIC FUNCTIONS ##
def make_model_irreversible(cobra_model: cobra.Model,
                            forward_id:str = "_FWD", reverse_id: str="_REV") -> cobra.Model:
    """Returns a cobrapy model where all reversible reactions are splitted into irreversible ones.

    Args:
        cobra_model (cobra.Model): The cobrapy model for which all reversible reactions shall be splitted.
        forward_id (str, optional): The forward reaction ID suffix used in new splitted reactions. Defaults to "_FWD".
        reverse_id (str, optional): The reverse reaction ID suffix used in new splitted reactions. Defaults to "_REV".

    Returns:
        cobra.Model: The cobrapy model with splitted reversible reactions.
    """
    cobra_reaction_ids = [x.id for x in cobra_model.reactions]
    for cobra_reaction_id in cobra_reaction_ids:
        cobra_reaction: cobra.Reaction = cobra_model.reactions.get_by_id(cobra_reaction_id)

        create_forward = False
        create_reverse = False
        if cobra_reaction.lower_bound < 0:
            create_reverse = True
        elif (cobra_reaction.lower_bound == 0) and (cobra_reaction.upper_bound == 0):
            create_reverse = True
            create_forward = True
        elif (cobra_reaction.id.startswith("EX_")):
            create_reverse = True
            create_forward = True
        else:
            continue

        if cobra_reaction.upper_bound > 0:
            create_forward = True

        if create_forward:
            forward_reaction_id = cobra_reaction.id + forward_id
            forward_reaction = copy.deepcopy(cobra_reaction)
            forward_reaction.id = forward_reaction_id
            if cobra_reaction.lower_bound >= 0:
                forward_reaction.lower_bound = cobra_reaction.lower_bound
            else:
                forward_reaction.lower_bound = 0
            cobra_model.add_reactions([forward_reaction])

        if create_reverse:
            reverse_reaction_id = cobra_reaction.id + reverse_id
            reverse_reaction = copy.deepcopy(cobra_reaction)
            reverse_reaction.id = reverse_reaction_id

            metabolites_to_add = {}
            for metabolite in reverse_reaction.metabolites:
                metabolites_to_add[metabolite] = reverse_reaction.metabolites[metabolite] * -2
            reverse_reaction.add_metabolites(metabolites_to_add)

            if cobra_reaction.upper_bound < 0:
                reverse_reaction.lower_bound = -cobra_reaction.upper_bound
            else:
                reverse_reaction.lower_bound = 0
            reverse_reaction.upper_bound = -cobra_reaction.lower_bound

            cobra_model.add_reactions([reverse_reaction])

        if create_forward or create_reverse:
            cobra_model.remove_reactions([cobra_reaction])

    return cobra_model
