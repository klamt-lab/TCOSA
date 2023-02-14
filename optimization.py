#!/usr/bin/env python3
"""[summary]"""

## IMPORTS ##
import ray
import pulp
from typing import Any, Dict, List, Tuple


## CONSTANTS ##
SOLVER = "CPLEX_PY"
TIMELIMIT = 180
IS_VERBOSE = True
WARMSTART = False


## PUBLIC FUNCTIONS ##
def set_timelimit(timelimit: int):
    global TIMELIMIT
    TIMELIMIT = timelimit


def add_objective(base_problem: pulp.LpProblem, objective: List[Tuple[float, str]], direction: str, variables_dict=None) -> pulp.LpProblem:
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
    ray.init(num_cpus=num_cpus, ignore_reinit_error=ignore_reinit_error, **kwargs)


def perform_optimization(base_problem: pulp.LpProblem, objective: List[Tuple[float, str]],
                         direction: str, variables_dict=None, **kwargs) -> Dict[str, Any]:
    """[summary]

    Args:
        base_problem (pulp.LpProblem): [description]
        objective (List[Tuple[float, str]]): [description]
        direction (str): [description]

    Returns:
        Dict[str, Any]: [description]
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
    """[summary]

    Args:
        base_problem (pulp.LpProblem): [description]

    Returns:
        Dict[str, Any]: [description]
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
    """[summary]

    Args:
        base_problem (pulp.LpProblem): [description]
        optimized_variable_id (str): [description]
        start_value (float, optional): [description]. Defaults to 0.0.
        step_sizes (List[float], optional): [description]. Defaults to [0.5, 0.2, 0.1, 0.05, 0.01].

    Returns:
        Dict[str, Any]: [description]
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
    """[summary]

    Args:
        base_problem (pulp.LpProblem): [description]
        variable_id (str): [description]

    Returns:
        Dict[str, Any]: [description]
    """
    return perform_optimization(
        base_problem=base_problem,
        objective=[(1.0, variable_id)],
        direction="max",
        variables_dict=variables_dict,
        **kwargs
    )


def perform_variable_minimization(base_problem: pulp.LpProblem, variable_id: str, variables_dict=None, **kwargs) -> Dict[str, Any]:
    """[summary]

    Args:
        base_problem (pulp.LpProblem): [description]
        variable_id (str): [description]

    Returns:
        Dict[str, Any]: [description]
    """
    return perform_optimization(
        base_problem=base_problem,
        objective=[(1.0, variable_id)],
        direction="min",
        variables_dict=variables_dict,
        **kwargs
    )


def solve_current_problem(base_problem: pulp.LpProblem) -> pulp.LpProblem:
    # solver = pulp.getSolver(SOLVER, msg=IS_VERBOSE,
    #                        warmStart=WARMSTART, timeLimit=TIMELIMIT)
    # base_problem.solve(solver)
    solver = pulp.CPLEX_PY()
    solver.buildSolverModel(base_problem)
    solver.solverModel.parameters.mip.display.set(1)
    solver.solverModel.parameters.paramdisplay.set(1)
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
