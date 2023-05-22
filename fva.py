#!/usr/bin/env python3
"""Functions for performing parallelized and serial variability analyses.

Note that these functions are heavily unoptimized. If you want to perform
much faster variability analyses, use other packages.
"""

## IMPORTS ##
# External
from typing import List
import cobra
import concurrent.futures
import pulp
import ray
from typing import Any, Dict, Tuple
# Internal imports
from optimization import perform_variable_maximization, perform_variable_minimization, solve_current_problem


## PUBLIC FUNCTIONS ##
def perform_variability_analysis(tested_vars: List[str],
                                 base_problem: pulp.LpProblem,
                                 **kwargs) -> Dict[str, Dict[str, float]]:
    """Performs a serial (non-parallelized) variability analysis for the given variables in the given problem.

    Each given variable is maximized and minimized.

    Args:
        tested_vars (List[str]): The list of variable IDs for which a variability analysis shall be performed.
        base_problem (pulp.LpProblem): The pulp problem in which the analysis shall happen.

    Returns:
        Dict[str, Dict[str, float]]: Contains as keys the tested variable IDs, and as value another dictionary
        containing the 'min' and the 'max' values with these keys.
    """
    # Test run of problem
    solve_current_problem(base_problem=base_problem)
    variables_dict = base_problem.variablesDict()

    # Go through each variable
    flux_bounds: Dict[str, Dict[str, float]] = {}
    for tested_var in tested_vars:
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


@ray.remote
def _perform_fva_multi(var_id: str, base_problem_dict: Dict[Any, Any]):
    """Ray subfunction of perform_fva_multi(...) (see there). """
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
                      base_problem: pulp.LpProblem) -> List[Tuple[str, float, float]]:
    """Performs a (hopefully) parallelized Flux Variability Analysis using ray for the given variables.

    Args:
        var_ids (List[str]): The list of (reaction) variable IDs for which the FVA shall be performed.
        base_problem (pulp.LpProblem): The pulp problem in which the FVA shall be performed.

    Returns:
        List[Tuple[str, float, float]]: Each tuplke contains at index 0 the variable ID, at 1 the respective minimal
        and at 2 the maximal value given by the FVA.
    """
    solve_current_problem(base_problem=base_problem)
    base_problem_dict = base_problem.to_dict()
    futures = [
        _perform_fva_multi.remote(
            var_id, base_problem_dict
        ) for var_id in var_ids
    ]
    results = ray.get(futures)

    return results
