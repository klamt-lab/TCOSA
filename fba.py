#!/usr/bin/env python3
"""This module contains all base functions for creating an FBA pulp problem from a cobrapy model.

In addition, it contains functions for maximization and minimization FBA tasks.
"""

## IMPORTS ##
# External
import cobra
import pulp
from typing import Any, Dict, List
# Internal
from optimization import perform_variable_maximization, perform_variable_minimization


## CONSTANTS ##
MAX_FLUX = 1_000  # The standard maximal flux
"""The standard maximal flux in mmol/(gDW*h)."""


## PUBLIC FUNCTIONS ##
def add_extra_constraints(problem: pulp.LpProblem, extra_constraints: List[Dict[str, float]], variables_dict=None) -> pulp.LpProblem:
    """Adds additional linear constraints to the pulp LpProblem.

    Args:
        problem (pulp.LpProblem): The pulp LpProblem to which the additional linear constraints are added.
        extra_constraints (List[Dict[str, float]]): A dictionary describing the additional linear constraints. As keys,
        it contains the names of the affected variables and 'lb' (the expression's lower bound) and/or 'ub' (upper bound).
        The dictionary's values are the linear coefficients of the variable and for 'lb' or 'ub', its respective value.
        variables_dict (_type_, optional): In order to save time to call problem.variablesDict(), its result
        can be send beforehand. Defaults to None.

    Returns:
        pulp.LpProblem: The LpProblem with the added additional linear constraints.
    """
    if variables_dict is None:
        variables_dict = problem.variablesDict()

    # Set extra constraint expressions
    current_extra_constraint = 0
    for constraint in extra_constraints:
        expression = pulp.LpAffineExpression(name=str(current_extra_constraint))
        lb = None
        ub = None
        for key in constraint.keys():
            # Lower bound
            if key == "lb":
                lb = constraint[key]
            # Upper bound
            elif key == "ub":
                ub = constraint[key]
            # If the key is (hopefully) a valid variable ID
            else:
                expression += constraint[key] * variables_dict[key]
        # Add lower bound if given
        if lb is not None:
            lb_constraint: pulp.LpConstraint = lb <= expression
            lb_constraint.name = f"lb_{current_extra_constraint}"
            problem += lb_constraint
        # Add upper bound if given
        if ub is not None:
            ub_constraint: pulp.LpConstraint = expression <= ub
            ub_constraint.name = f"ob_{current_extra_constraint}"
            problem += ub_constraint
        current_extra_constraint += 1

    return problem


def get_fba_base_problem(cobra_model: cobra.Model, extra_constraints: List[Dict[str, float]], name: str="Flux_Balance_Analysis") -> pulp.LpProblem:
    """Returns an FBA 'base problem', i.e., N*v=0 with N as stoichiometric matrix and v as the reaction flux vector with upper and lower flux bounds.

    Note that the typical FBA objective is missing and has to be added afterwards using pulp or this module's optimization functions.

    Args:
        cobra_model (cobra.Model): The cobrapy model from which the stoichiometric matrix and the flux bounds are build up.
        extra_constraints (List[Dict[str, float]]): [Optional] Additional linear constraints. See the comment to this module's add_extra_constraints(...) for more.
        name (str, optional): [Optional] An internal name for the returned pulp LpProblem. Defaults to "Flux_Balance_Analysis".

    Returns:
        pulp.LpProblem: A pulp LpProblem describing the FBA 'base problem'.
    """
    # Set problem instance
    problem = pulp.LpProblem(name=name, sense=pulp.LpMinimize)

    # Add reaction flux variables and set lower and upper fluxes for them
    flux_variables: Dict[str, pulp.LpVariable] = {}
    metabolite_expressions: Dict[str, List[pulp.LpAffineExpression]] = {}
    for reaction in cobra_model.reactions:
        lowBound = max(-MAX_FLUX, reaction.lower_bound)
        lowBound = min(lowBound, MAX_FLUX)
        upBound = min(MAX_FLUX, reaction.upper_bound)
        upBound = max(-MAX_FLUX, upBound)

        flux_variables[reaction.id] =  pulp.LpVariable(
            name=reaction.id,
            lowBound=lowBound,
            upBound=upBound,
            cat=pulp.LpContinuous,
        )

        for metabolite in reaction.metabolites.keys():
            stoichiometry = reaction.metabolites[metabolite]
            if metabolite.id not in metabolite_expressions.keys():
                metabolite_expressions[metabolite.id] = []
            metabolite_expressions[metabolite.id] += stoichiometry * flux_variables[reaction.id]

    # Set steady-state expressions
    for metabolite_id in metabolite_expressions.keys():
        constraint: pulp.LpConstraint = metabolite_expressions[metabolite_id] == 0
        constraint.name = metabolite_id
        problem += constraint

    # Add additional linear constraints
    problem = add_extra_constraints(problem, extra_constraints, flux_variables)

    return problem


def perform_fba_flux_maximization(base_problem: pulp.LpProblem, reaction_id: str, **kwargs) -> Dict[str, Any]:
    """Maximizes the flux of the given reaction in the given FBA base problem.

    The base problem is expected to have been created with this module's get_fba_base_problem().

    Args:
        base_problem (pulp.LpProblem): The FBA base problem.
        reaction_id (str): The ID of the reaction whose flux shall be maximized.

    Returns:
        Dict[str, Any]: The optimization results as given in perform_variable_maximization(...).
    """
    return perform_variable_maximization(
        base_problem=base_problem,
        variable_id=reaction_id,
        **kwargs,
    )


def perform_fba_flux_minimization(base_problem: pulp.LpProblem, reaction_id: str, **kwargs) -> Dict[str, Any]:
    """Minimizes the flux of the given reaction in the given FBA base problem.

    The base problem is expected to have been created with this module's get_fba_base_problem().

    Args:
        base_problem (pulp.LpProblem): The FBA base problem.
        reaction_id (str): The ID of the reaction whose flux shall be minimized.

    Returns:
        Dict[str, Any]: The optimization results as given in perform_variable_maximization(...).
    """
    return perform_variable_minimization(
        base_problem=base_problem,
        variable_id=reaction_id,
        **kwargs,
    )
