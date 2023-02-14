#!/usr/bin/env python3
"""[summary]"""

## IMPORTS ##
# External
import cobra
import pulp
from optimization import perform_variable_maximization, perform_variable_minimization
from typing import Any, Dict, List


## CONSTANTS ##
MAX_FLUX = 1_000


## PUBLIC FUNCTIONS ##
def add_extra_constraints(problem: pulp.LpProblem, extra_constraints: List[Dict[str, float]], variables_dict=None) -> pulp.LpProblem:
    if variables_dict is None:
        variables_dict = problem.variablesDict()

    # Set extra constraint expressions
    current_extra_constraint = 0
    for constraint in extra_constraints:
        expression = pulp.LpAffineExpression(name=str(current_extra_constraint))
        lb = None
        ub = None
        for key in constraint.keys():
            if key == "lb":
                lb = constraint[key]
            elif key == "ub":
                ub = constraint[key]
            else:  # If the key is (hoperandoms) a valid variable ID
                expression += constraint[key] * variables_dict[key]
        if lb is not None:
            lb_constraint: pulp.LpConstraint = lb <= expression
            lb_constraint.name = f"lb_{current_extra_constraint}"
            problem += lb_constraint
        if ub is not None:
            ub_constraint: pulp.LpConstraint = expression <= ub
            ub_constraint.name = f"ob_{current_extra_constraint}"
            problem += ub_constraint
        current_extra_constraint += 1

    return problem


def get_fba_base_problem(cobra_model: cobra.Model, extra_constraints: List[Dict[str, float]], name: str="Flux_Balance_Analysis") -> pulp.LpProblem:
    """[summary]

    Args:
        cobra_model (cobra.Model): [description]
        extra_constraints (List[Dict[str, float]]): [description]

    Returns:
        pulp.LpProblem: [description]
    """
    # Set problem instance
    problem = pulp.LpProblem(name=name, sense=pulp.LpMinimize)

    # Get flux variables
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

    problem = add_extra_constraints(problem, extra_constraints, flux_variables)

    return problem


def perform_fba_flux_maximization(base_problem: pulp.LpProblem, reaction_id: str, **kwargs) -> Dict[str, Any]:
    """[summary]

    Args:
        base_problem (cvxpy.Problem): [description]
        reaction_id (str): [description]

    Returns:
        Dict[str, Any]: [description]
    """
    return perform_variable_maximization(
        base_problem=base_problem,
        variable_id=reaction_id,
        **kwargs,
    )


def perform_fba_flux_minimization(base_problem: pulp.LpProblem, reaction_id: str, **kwargs) -> Dict[str, Any]:
    """[summary]

    Args:
        base_problem (cvxpy.Problem): [description]
        reaction_id (str): [description]

    Returns:
        Dict[str, Any]: [description]
    """
    return perform_variable_minimization(
        base_problem=base_problem,
        variable_id=reaction_id,
        **kwargs,
    )
