#!/usr/bin/env python3
"""[summary]"""

## IMPORTS ##
# External
from typing import List
import cobra
import concurrent.futures
import pulp
import ray
from typing import Any, Dict
# Internal imports
from optimization import perform_variable_maximization, perform_variable_minimization, solve_current_problem


## PUBLIC FUNCTIONS ##
def perform_variability_analysis(tested_vars: List[str],
                                 base_problem: pulp.LpProblem,
                                 **kwargs) -> Dict[str, Dict[str, float]]:
    """[summary]

    Args:
        cobra_model (cobra.Model): [description]
        base_problem (pulp.LpProblem): [description]

    Returns:
        Dict[str, Dict[str, float]]: [description]
    """
    solve_current_problem(base_problem=base_problem)
    variables_dict = base_problem.variablesDict()

    flux_bounds: Dict[str, Dict[str, float]] = {}
    for tested_var in tested_vars:
        # print(reaction.id, "...")
        flux_bounds[tested_var] = {}

        # Minimization
        min_result = perform_variable_minimization(
            base_problem=base_problem,
            variable_id=tested_var,
            variables_dict=variables_dict,
            **kwargs,
        )
        if min_result["status"] == "Optimal":
            min_value = min_result["values"][tested_var]
        else:
            min_value = float("NaN")
        flux_bounds[tested_var]["min"] = min_value

        # Maximization
        max_result = perform_variable_maximization(
            base_problem=base_problem,
            variable_id=tested_var,
            variables_dict=variables_dict,
            **kwargs,
        )
        if max_result["status"] == "Optimal":
            max_value = max_result["values"][tested_var]
        else:
            max_value = float("NaN")
        flux_bounds[tested_var]["max"] = max_value

    return flux_bounds


def perform_fva(cobra_model: cobra.Model,
                base_problem: pulp.LpProblem,
                test,
                **kwargs) -> Dict[str, Dict[str, float]]:
    """[summary]

    Args:
        cobra_model (cobra.Model): [description]
        base_problem (pulp.LpProblem): [description]

    Returns:
        Dict[str, Dict[str, float]]: [description]
    """
    solve_current_problem(base_problem=base_problem)
    variables_dict = base_problem.variablesDict()

    flux_bounds: Dict[str, Dict[str, float]] = {}
    for reaction in cobra_model.reactions:
        if not (("_COSA_" in reaction.id) or ("Biomass" in reaction.id)):
            continue
        if reaction.id not in test:
            continue
        # print(reaction.id, "...")
        flux_bounds[reaction.id] = {}

        # Minimization
        min_result = perform_variable_minimization(
            base_problem=base_problem,
            variable_id=reaction.id,
            variables_dict=variables_dict,
            **kwargs,
        )
        if min_result["status"] == "Optimal":
            min_value = min_result["values"][reaction.id]
        else:
            min_value = float("NaN")
        flux_bounds[reaction.id]["min"] = min_value

        # Maximization
        max_result = perform_variable_maximization(
            base_problem=base_problem,
            variable_id=reaction.id,
            variables_dict=variables_dict,
            **kwargs,
        )
        if max_result["status"] == "Optimal":
            max_value = max_result["values"][reaction.id]
        else:
            max_value = float("NaN")
        flux_bounds[reaction.id]["max"] = max_value

    return flux_bounds


@ray.remote
def _perform_fva_multi(var_id: str, base_problem_dict: Dict[Any, Any]):
    variables_dict, base_problem = pulp.LpProblem.from_dict(base_problem_dict)

    # Minimization
    min_result = perform_variable_minimization(
        base_problem=base_problem,
        variable_id=var_id,
        variables_dict=variables_dict,
    )
    if min_result["status"] == "Optimal":
        min_value = min_result["values"][var_id]
    else:
        min_value = float("NaN")

    # Maximization
    max_result = perform_variable_maximization(
        base_problem=base_problem,
        variable_id=var_id,
        variables_dict=variables_dict,
    )
    if max_result["status"] == "Optimal":
        max_value = max_result["values"][var_id]
    else:
        max_value = float("NaN")

    return (var_id, min_value, max_value)


def perform_fva_multi(var_ids: List[str],
                      base_problem: pulp.LpProblem) -> Dict[str, Dict[str, float]]:
    solve_current_problem(base_problem=base_problem)
    base_problem_dict = base_problem.to_dict()
    futures = [
        _perform_fva_multi.remote(
            var_id, base_problem_dict
        ) for var_id in var_ids
    ]
    results = ray.get(futures)

    return results

###########################################
###########################################


"""
def _mt_perform_fva_multi(reaction: cobra.Reaction, base_problem_dict: Dict[Any, Any]):
    variables_dict, base_problem = pulp.LpProblem.from_dict(base_problem_dict)

    # Minimization
    min_result = perform_variable_minimization(
        base_problem=base_problem,
        variable_id=reaction.id,
        variables_dict=variables_dict,
    )
    if min_result["status"] == "Optimal":
        min_value = min_result["values"][reaction.id]
    else:
        min_value = float("NaN")

    # Maximization
    max_result = perform_variable_maximization(
        base_problem=base_problem,
        variable_id=reaction.id,
        variables_dict=variables_dict,
    )
    if max_result["status"] == "Optimal":
        max_value = max_result["values"][reaction.id]
    else:
        max_value = float("NaN")

    return (reaction.id, min_value, max_value)


@ray.remote
def _ray_perform_fva_multi(reactions: List[cobra.Reaction], base_problem_dict: Dict[Any, Any]):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(
            _mt_perform_fva_multi, reaction, base_problem_dict) for reaction in reactions]
        results = [future.result() for future in futures]

    return results


def perform_fva_multi(cobra_model: cobra.Model,
                      base_problem: pulp.LpProblem) -> Dict[str, Dict[str, float]]:
    base_problem_dict = base_problem.to_dict()

    batch_size = 10
    loop_values = []
    start_value = 0
    while start_value < len(cobra_model.reactions):
        loop_values.append(start_value)
        start_value += batch_size
    futures = [_ray_perform_fva_multi.remote(
        cobra_model.reactions[loop_value:loop_value+batch_size], base_problem_dict) for loop_value in loop_values]
    results = ray.get(futures)

    full_results = []
    for result in results:
        full_results += result

    return full_results
"""
