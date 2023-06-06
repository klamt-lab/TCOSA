"""Contains function for adding promiscuity constraints."""

# IMPORTS #
# External
import cobra
import copy
import pulp
from typing import Any, Dict

# PUBLIC FUNCTIONS #
def cosa_add_promiscuity_constraints(
    optmdfpathway_base_problem: pulp.LpProblem, optmdfpathway_base_variables: Dict[Any, Any],
    cobra_model: cobra.Model, dG0_values: Dict[Any, Any]) -> pulp.LpProblem:
    """This function adds the TCOSA promiscuity constraints to the pulp problem.

    The promuscuity constraints ensure that out of all generated forward and reverse
    and TCOSA reactions, only one of the descendants of an original iML1515 reaction
    can be active at the same time.

    Args:
        optmdfpathway_base_problem (pulp.LpProblem): The pulp problem.
        optmdfpathway_base_variables (Dict[Any, Any]): The pulp problem's base variables dict.
        cobra_model (cobra.Model): The cobra model on which the pulp problem is based upon.
        dG0_values (Dict[Any, Any]): The used dG0 values.

    Returns:
        pulp.LpProblem: The pulp problem enhanced by the promiscuity constraints.
    """

    for reaction in cobra_model.reactions:
        key = reaction.id

        if key not in dG0_values.keys():
            continue

        newkeys = [key]
        if "_FWD" in key:
            newkeys.append(key.replace("_FWD", "_REV"))
        elif "_REV" in key:
            newkeys.append(key.replace("_REV", "_FWD"))

        if "_ORIGINAL_NAD_" in key:
            for newkey in copy.deepcopy(newkeys):
                newkeys.append(newkey.replace("_ORIGINAL_NAD_", "_VARIANT_NADP_"))
        elif "_ORIGINAL_NADP_" in key:
            for newkey in copy.deepcopy(newkeys):
                newkeys.append(newkey.replace("_ORIGINAL_NADP_", "_VARIANT_NAD_"))
        elif "_VARIANT_NAD_" in key:
            for newkey in copy.deepcopy(newkeys):
                newkeys.append(newkey.replace("_VARIANT_NAD_", "_ORIGINAL_NADP_"))
        elif "_VARIANT_NADP_" in key:
            for newkey in copy.deepcopy(newkeys):
                newkeys.append(newkey.replace("_VARIANT_NADP_", "_ORIGINAL_NAD_"))

        #####
        if len(newkeys) <= 1:
            continue

        current_z_cluster_sum: pulp.LpAffineExpression = 0.0
        for newkey in newkeys:
            current_z_variable = optmdfpathway_base_variables[f"z_var_{newkey}"]
            current_z_cluster_sum += current_z_variable

        var_current_z_cluster_sum = pulp.LpVariable(
            name=f"var_z_cluster_sum_{newkeys[0]}",
            lowBound=0.0,
            upBound=1.0,
            cat=pulp.LpContinuous,
        )
        optmdfpathway_base_problem += current_z_cluster_sum - var_current_z_cluster_sum == 0
        optmdfpathway_base_problem += var_current_z_cluster_sum <= 1

    return optmdfpathway_base_problem


def cosa_add_stoichiometric_promiscuity_constraints(
    optmdfpathway_base_problem: pulp.LpProblem, optmdfpathway_base_variables: Dict[Any, Any],
    cobra_model: cobra.Model, dG0_values: Dict[Any, Any]) -> pulp.LpProblem:
    """This function adds the TCOSA promiscuity constraints to the pulp problem in the case that no thermodynamic constraints are used.

    The promuscuity constraints ensure that out of all generated forward and reverse
    and TCOSA reactions, only one of the descendants of an original iML1515 reaction
    can be active at the same time.

    Args:
        optmdfpathway_base_problem (pulp.LpProblem): The pulp problem.
        optmdfpathway_base_variables (Dict[Any, Any]): The pulp problem's base variables dict.
        cobra_model (cobra.Model): The cobra model on which the pulp problem is based upon.
        dG0_values (Dict[Any, Any]): The used dG0 values.

    Returns:
        pulp.LpProblem: The pulp problem enhanced by the promiscuity constraints.
    """

    current_cluster = 0
    for reaction in cobra_model.reactions:
        key = reaction.id

        newkeys = [key]
        if "_FWD" in key:
            newkeys.append(key.replace("_FWD", "_REV"))
        elif "_REV" in key:
            newkeys.append(key.replace("_REV", "_FWD"))

        if "_ORIGINAL_NAD_" in key:
            for newkey in copy.deepcopy(newkeys):
                newkeys.append(newkey.replace("_ORIGINAL_NAD_", "_VARIANT_NADP_"))
        elif "_ORIGINAL_NADP_" in key:
            for newkey in copy.deepcopy(newkeys):
                newkeys.append(newkey.replace("_ORIGINAL_NADP_", "_VARIANT_NAD_"))
        elif "_VARIANT_NAD_" in key:
            for newkey in copy.deepcopy(newkeys):
                newkeys.append(newkey.replace("_VARIANT_NAD_", "_ORIGINAL_NADP_"))
        elif "_VARIANT_NADP_" in key:
            for newkey in copy.deepcopy(newkeys):
                newkeys.append(newkey.replace("_VARIANT_NADP_", "_ORIGINAL_NAD_"))

        #####
        if len(newkeys) <= 1:
            continue

        current_cluster += 1
        current_z_cluster_sum: pulp.LpAffineExpression = 0.0
        for newkey in newkeys:
            new_z_cluster_variable = pulp.LpVariable(
                name=f"var_z_cluster_var_{newkey}_{current_cluster}",
                lowBound=0.0,
                upBound=1.0,
                cat=pulp.LpBinary,
            )
            optmdfpathway_base_problem += optmdfpathway_base_variables[newkey] <= 1_500 * new_z_cluster_variable
            current_z_cluster_sum += new_z_cluster_variable

        var_current_z_cluster_sum = pulp.LpVariable(
            name=f"var_z_cluster_sum_{newkeys[0]}_{current_cluster}",
            lowBound=0.0,
            upBound=1.0,
            cat=pulp.LpBinary,
        )
        optmdfpathway_base_problem += current_z_cluster_sum - var_current_z_cluster_sum == 0
        optmdfpathway_base_problem += var_current_z_cluster_sum <= 1

    return optmdfpathway_base_problem
