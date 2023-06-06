#!/usr/bin/env python3
"""Generic functions for performing optimizations, i.e., maximizations and/or optimizations.

Basically, these are pulp wrappers.
"""

## IMPORTS ##
import ray
import pulp
from typing import Any, Dict, List, Tuple


## CONSTANTS ##
SOLVER = "CPLEX_PY"  # Can be whatever pulp supports.
TIMELIMIT = 180  # In seconds
IS_VERBOSE = True
WARMSTART = False


## PUBLIC FUNCTIONS ##
def set_timelimit(timelimit: int):
    """Sets the solver timelimit to the given value."""
    global TIMELIMIT
    TIMELIMIT = timelimit


def add_objective(base_problem: pulp.LpProblem, objective: List[Tuple[float, str]], direction: str, variables_dict=None) -> pulp.LpProblem:
    """Adds the given objective to the given pulp problem

    Args:
        base_problem (pulp.LpProblem): The pulp problem.
        objective (List[Tuple[float, str]]): The objective in the form of tuples where element 0 is the coefficient of the
        objective functions and element 1 the variable's ID.
        direction (str): The objective's direction, i.e., 'min' or 'max'.
        variables_dict (_type_, optional): Optional already calculated variablesDict from pulp. Defaults to None.

    Returns:
        pulp.LpProblem: _description_
    """
    if variables_dict is None:
        base_problem_variables: Dict[str,
                                     pulp.LpVariable] = base_problem.variablesDict()
    else:
        base_problem_variables = variables_dict

    objective_expression: pulp.LpAffineExpression = 0.0
    for objective_part in objective:
        multiplier = objective_part[0]
        variable_id = objective_part[1]

        if variable_id == "NO_OPTIMIZATION":
            break

        objective_expression += multiplier * \
            base_problem_variables[variable_id]

    if variable_id != "NO_OPTIMIZATION":
        final_objective: pulp.LpAffineExpression = None
        if direction == "min":
            final_objective = objective_expression
        elif direction == "max":
            final_objective = -objective_expression
        final_objective.name = "objective"
        base_problem.objective = final_objective

    return base_problem


def init_ray_multiprocessing(num_cpus=16, ignore_reinit_error=True, **kwargs):
    """Initialize ray with the given settings."""
    ray.init(num_cpus=num_cpus, ignore_reinit_error=ignore_reinit_error, **kwargs)


def perform_optimization(base_problem: pulp.LpProblem, objective: List[Tuple[float, str]],
                         direction: str, variables_dict=None, **kwargs) -> Dict[str, Any]:
    """Performs the given optimization.

    Args:
        base_problem (pulp.LpProblem): The pulp problem.
        objective (List[Tuple[float, str]]): The objective in the form of tuples where element 0 is the coefficient of the
        objective functions and element 1 the variable's ID.
        direction (str): The objective's direction, i.e., 'min' or 'max'.
        variables_dict (_type_, optional): Optional already calculated variablesDict from pulp. Defaults to None.

    Returns:
        Dict[str, Any]: Returns a dict with 'status' as key (there, the solver status is the value) as well as
        'objective_value' (its namesake) and 'values' (another dictionary with all variable IDs as keys and their
        solution values as values).
    """
    if variables_dict is None:
        base_problem_variables: Dict[str,
                                     pulp.LpVariable] = base_problem.variablesDict()
    else:
        base_problem_variables = variables_dict

    add_objective(
        base_problem=base_problem,
        objective=objective,
        direction=direction,
        variables_dict=variables_dict,
    )

    base_problem = solve_current_problem(base_problem=base_problem)

    results: Dict[str, Any] = {}
    results["status"] = pulp.LpStatus[base_problem.status]
    results["objective_value"] = base_problem.objective.value()
    results["values"] = {}
    for variable_id in base_problem_variables.keys():
        variable_value = base_problem_variables[variable_id].value()
        results["values"][variable_id] = variable_value

    return results


def perform_optimization_with_given_objective(base_problem: pulp.LpProblem, variables_dict=None, **kwargs) -> Dict[str, Any]:
    """Performs optimization with pre-added objective.

    Args:
        base_problem (pulp.LpProblem): The pulp problem.
        variables_dict (_type_, optional): Optional already calculated variablesDict from pulp. Defaults to None.

    Returns:
        Dict[str, Any]: Returns a dict with 'status' as key (there, the solver status is the value) as well as
        'objective_value' (its namesake) and 'values' (another dictionary with all variable IDs as keys and their
        solution values as values).
    """
    return perform_optimization(
        base_problem=base_problem,
        objective=[(0.0, "NO_OPTIMIZATION")],
        direction="min",
        variables_dict=variables_dict,
        **kwargs
    )


def perform_stepwise_variable_optimization(
    base_problem: pulp.LpProblem,
    optimized_variable_id: str,
    start_value: float = 0.0,
    step_sizes: List[float] = [0.5, 0.2, 0.1, 0.05, 0.01],
    variables_dict=None,
    **kwargs
) -> Dict[str, Any]:
    """The same as perform_optimization(), but only performs simplex phase I by trying to go closer to the objective with the given step sizes.

    Args:
        base_problem (pulp.LpProblem): The pulp problem.
        optimized_variable_id (str): The optimized variable ID.
        start_value (float, optional): The initial value for the optimized variable. Defaults to 0.0.
        step_sizes (List[float], optional): The phase I approximation step sizes. Defaults to [0.5, 0.2, 0.1, 0.05, 0.01].

    Returns:
        Dict[str, Any]: Returns a dict with 'status' as key (there, the solver status is the value) as well as
        'objective_value' (its namesake) and 'values' (another dictionary with all variable IDs as keys and their
        solution values as values).
    """
    if variables_dict is None:
        base_problem_variables = base_problem.variablesDict()
    else:
        base_problem_variables = variables_dict
    variable: pulp.LpVariable = base_problem_variables[optimized_variable_id]

    last_working_value = start_value
    for step_size in step_sizes:
        current_status = "Optimal"
        current_min_value = last_working_value + step_size
        while current_status == "Optimal":
            variable.lowBound = current_min_value
            result = perform_optimization_with_given_objective(
                base_problem=base_problem, variables_dict=base_problem_variables, **kwargs)
            current_status = result["status"]

            if current_status == "Optimal":
                last_working_value = current_min_value

            current_min_value += step_size

    variable.lowBound = last_working_value
    output_result = perform_optimization_with_given_objective(
        base_problem=base_problem, variables_dict=base_problem_variables, **kwargs)

    return output_result


def perform_variable_maximization(base_problem: pulp.LpProblem, variable_id: str, variables_dict=None, **kwargs) -> Dict[str, Any]:
    """Maximizes the given variable with perform_optimization(...)."""
    return perform_optimization(
        base_problem=base_problem,
        objective=[(1.0, variable_id)],
        direction="max",
        variables_dict=variables_dict,
        **kwargs
    )


def perform_variable_minimization(base_problem: pulp.LpProblem, variable_id: str, variables_dict=None, **kwargs) -> Dict[str, Any]:
    """Minimizes the given variable with perform_optimization(...)."""
    return perform_optimization(
        base_problem=base_problem,
        objective=[(1.0, variable_id)],
        direction="min",
        variables_dict=variables_dict,
        **kwargs
    )


def solve_current_problem(base_problem: pulp.LpProblem) -> pulp.LpProblem:
    """Solves the given problem using pulp.

    Args:
        base_problem (pulp.LpProblem): The pulp problem.

    Returns:
        pulp.LpProblem: The solved pulp problem.
    """
    # solver = pulp.getSolver(SOLVER, msg=IS_VERBOSE,
    #                        warmStart=WARMSTART, timeLimit=TIMELIMIT)
    # base_problem.solve(solver)
    solver = pulp.CPLEX_PY()
    solver.buildSolverModel(base_problem)
    solver.solverModel.parameters.mip.display.set(0)
    solver.solverModel.parameters.paramdisplay.set(0)
    solver.solverModel.parameters.simplex.tolerances.feasibility.set(1e-9)
    # solver.solverModel.parameters.simplex.tolerances.optimality.set(1e-6)
    solver.solverModel.parameters.emphasis.numerical.set(1)
    solver.solverModel.parameters.mip.tolerances.absmipgap.set(1e-9)
    solver.solverModel.parameters.mip.tolerances.integrality.set(0.0)
    solver.solverModel.parameters.timelimit.set(TIMELIMIT)
    # solver.solverModel.parameters.read.scale.set(-1)
    solver.callSolver(base_problem)
    status = solver.findSolutionValues(base_problem)

    return base_problem
