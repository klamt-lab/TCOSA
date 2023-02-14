#!/usr/bin/env python3
"""[summary]"""

## IMPORTS ##
# External
import cobra
import copy
import pulp
from math import exp, log
from typing import Any, List, Dict, Tuple

## CONSTANTS ##
M = 10_000
STANDARD_R = 8.314e-3  # kJ⋅K⁻1⋅mol⁻1 (standard value is in J⋅K⁻1⋅mol⁻1)
STANDARD_T = 298.15  # K


## PUBLIC FUNCTIONS ##
# [
#   [("var_B", 0.01, 0.11), ...],
#   ...
# ]
def get_multoptmdfpathway_base_problem(cobra_model: cobra.Model, dG0_values: Dict[str, Dict[str, float]],
                                       metabolite_concentration_values: Dict[str, Dict[str, float]],
                                       ratio_constraint_data: List[Dict[str, Any]],
                                       R: float, T: float,
                                       sub_network_ids: List[str],
                                       conditions: List[List[Tuple[str, float, float]]],
                                       condition_names: List[str]=[],
                                       add_promiscuity_constraints: bool=False) -> pulp.LpProblem:
    # Concentrations in M
    # T in K
    # R in kJ⋅K⁻1⋅mol⁻1 (standard value is in J⋅K⁻1⋅mol⁻1)
    # Get FBA base
    # Set problem instance
    base_problem = pulp.LpProblem(name="MultOptMDF", sense=pulp.LpMinimize)

    ### SET BASE FBA ###
    condition_counter = 0
    flux_variables: Dict[str, pulp.LpVariable] = {}
    concentration_vars: Dict[str, pulp.LpVariable] = {}
    global_z_vars: Dict[str, pulp.LpVariable] = {}
    local_z_vars: Dict[str, pulp.LpVariable] = {}
    dG0_vars: Dict[str, pulp.LpVariable] = {}

    global_tcosa_z_var_sum_expression = 0.0
    for condition in conditions:
        metabolite_expressions: Dict[str, List[pulp.LpAffineExpression]] = {}

        if condition_names == []:
            suffix = f"_COND{condition_counter}"
        else:
            suffix = "_"+condition_names[condition_counter]
        condition_counter += 1

        # Get flux variables
        for reaction in cobra_model.reactions:
            lowBound = max(-1_000, reaction.lower_bound)
            lowBound = min(lowBound, 1_000)
            upBound = min(1_000, reaction.upper_bound)
            upBound = max(-1_000, upBound)

            flux_variables[reaction.id+suffix] =  pulp.LpVariable(
                name=reaction.id+suffix,
                lowBound=lowBound,
                upBound=upBound,
                cat=pulp.LpContinuous,
            )

            for metabolite in reaction.metabolites.keys():
                stoichiometry = reaction.metabolites[metabolite]
                if (metabolite.id+suffix) not in metabolite_expressions.keys():
                    metabolite_expressions[metabolite.id+suffix] = 0.0
                metabolite_term = stoichiometry * flux_variables[reaction.id+suffix]
                metabolite_expressions[metabolite.id+suffix] += metabolite_term

        # Set steady-state expressions
        for metabolite_id in metabolite_expressions.keys():
            constraint: pulp.LpConstraint = metabolite_expressions[metabolite_id] == 0
            constraint.name = metabolite_id
            base_problem += constraint

        # Min/max concentrations
        for metabolite in cobra_model.metabolites:
            if metabolite.id in metabolite_concentration_values.keys():
                id_key = metabolite.id
            else:
                id_key = "DEFAULT"

            min_concentration = log(metabolite_concentration_values[id_key]["min"])
            max_concentration = log(metabolite_concentration_values[id_key]["max"])

            concentration_vars[metabolite.id+suffix] = pulp.LpVariable(
                name="x_"+metabolite.id+suffix,
                lowBound=min_concentration,
                upBound=max_concentration,
                cat=pulp.LpContinuous,
            )

        # Concentration ratios
        ratio_counter = 0
        for ratio_constraint in ratio_constraint_data:
            # c_i / c_j <= h_max AND c_i / c_j >= h_min
            # <=> x_i - x_j <= ln(h_max) AND x_i - x_j >= ln(h_min)
            # <=> (A) x_i - x_j - ln(h_max) <= 0 AND (B) -x_i + x_j - ln(h_min) <= 0
            c_i_id = ratio_constraint["c_i"]
            c_j_id = ratio_constraint["c_j"]

            ln_h_min = log(ratio_constraint["h_min"])
            ln_h_max = log(ratio_constraint["h_max"])

            # (A) x_i - x_j - ln(h_max) <= 0
            max_ratio_constraint: pulp.LpConstraint = concentration_vars[c_i_id+suffix] - concentration_vars[c_j_id+suffix] - ln_h_max <= 0
            # (B) -x_i + x_j - ln(h_min) <= 0
            min_ratio_constraint: pulp.LpConstraint = -concentration_vars[c_i_id+suffix] + concentration_vars[c_j_id+suffix] - ln_h_min <= 0

            base_problem += max_ratio_constraint
            base_problem += min_ratio_constraint

            ratio_counter += 1

        # df variables (f_r)
        var_B = pulp.LpVariable(
            name="var_B"+suffix,
            lowBound=-float("inf"),
            upBound=float("inf"),
            cat=pulp.LpContinuous,
        )

        has_subnetwork = len(sub_network_ids) > 0
        if has_subnetwork:
            var_B2 = pulp.LpVariable(
                name="var_B2"+suffix,
                lowBound=-float("inf"),
                upBound=float("inf"),
                cat=pulp.LpContinuous,
            )

        for reaction in cobra_model.reactions:
            if reaction.id not in dG0_values.keys():
                continue

            current_f_variable = pulp.LpVariable(
                name=f"f_var_{reaction.id}{suffix}",
                lowBound=-float("inf"),
                upBound=float("inf"),
                cat=pulp.LpContinuous,
            )

            # f_r = -(ΔG'0_r + RT * S^(t)_.,r * x)
            # <=> -ΔG'0_r - - f_r RT * S^(t)_.,r * x = 0
            dG0_value = dG0_values[reaction.id]["dG0"]

            if reaction.id not in dG0_vars.keys():
                dG0_vars[reaction.id] = pulp.LpVariable(
                    name=f"dG0_{reaction.id}",
                    cat=pulp.LpContinuous,
                )
                dG0_vars[reaction.id].setInitialValue(dG0_value)
                dG0_vars[reaction.id].fixValue()
            dG0_randomfixed_variable = dG0_vars[reaction.id]

            # e.g.,
            #     RT * ln(([S]*[T])/([A]²*[B]))
            # <=> RT*ln([S]) + RT*ln([T]) - 2*RT*ln([A]) - RT*ln([B])
            reaction_f_expression = -current_f_variable
            for metabolite in reaction.metabolites:
                stoichiometry = reaction.metabolites[metabolite]
                reaction_f_expression -= stoichiometry * R * T * concentration_vars[metabolite.id+suffix]
            reaction_f_expression -= dG0_randomfixed_variable
            reaction_f_constraint: pulp.LpConstraint = reaction_f_expression == 0

            base_problem += reaction_f_constraint

            # OptMDFpathway binary variable introduction
            global_z_varname = f"z_var_{reaction.id}"
            local_z_varname = global_z_varname+suffix
            if global_z_varname not in global_z_vars.keys():
                global_z_vars[global_z_varname] = pulp.LpVariable(
                    name=global_z_varname,
                    cat=pulp.LpBinary,
                )
            if local_z_varname not in local_z_vars.keys():
                local_z_vars[local_z_varname] = pulp.LpVariable(
                    name=local_z_varname,
                    cat=pulp.LpBinary,
                )
            current_global_z_variable = global_z_vars[global_z_varname]
            current_local_z_variable = local_z_vars[local_z_varname]

            # Find flux variable
            current_flux_variable = flux_variables[reaction.id+suffix]

            # z_r = 0 -> v_r = 0
            # In Big M form: v_r <= M*z_r <=> v_r - M*z_r <= 0
            z_zero_constraint: pulp.LpConstraint = current_flux_variable - M*current_local_z_variable <= 0
            # z_r = 1 -> f_r >= B == f_r - B >= 0
            # In Big M form: f_r + (1-z_i)*M >= B
            z_one_constraint = var_B <= current_f_variable + (1-current_local_z_variable)*M

            base_problem += z_zero_constraint
            base_problem += z_one_constraint

            if has_subnetwork:
                if reaction.id in sub_network_ids:
                    z_one_constraint_2: pulp.LpConstraint = var_B2 <= current_f_variable + (1-current_local_z_variable)*M
                    base_problem += z_one_constraint_2

            # Make local z var dependent on global z var
            base_problem += current_local_z_variable <= current_global_z_variable

        counter = 0
        all_variables = base_problem.variablesDict()
        for var_setting in condition:
            counter += 1
            var_id = var_setting[0] + suffix
            low_bound = max(-1000, var_setting[1])
            up_bound = min(1000, var_setting[2])
            all_variables[var_id].lowBound = low_bound
            all_variables[var_id].upBound = up_bound

        if add_promiscuity_constraints:
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

                current_local_z_cluster_sum: pulp.LpAffineExpression = 0.0
                for newkey in newkeys:
                    current_local_z_variable = local_z_vars[f"z_var_{newkey}{suffix}"]
                    current_local_z_cluster_sum += current_local_z_variable

                var_current_local_z_cluster_sum = pulp.LpVariable(
                    name=f"var_local_z_cluster_sum_{newkeys[0]}{suffix}",
                    lowBound=0.0,
                    upBound=1.0,
                    cat=pulp.LpContinuous,
                )
                base_problem += current_local_z_cluster_sum - var_current_local_z_cluster_sum == 0
                base_problem += var_current_local_z_cluster_sum <= 1

    global_z_var_sum = pulp.LpVariable(
        name="global_z_var_sum",
        lowBound=-float("inf"),
        upBound=float("inf"),
    )
    global_z_var_sum_expression = 0.0
    for global_z_var in global_z_vars.values():
        global_z_var_sum_expression += global_z_var
        if has_subnetwork:
            reac_id = global_z_var.name.replace("z_var_", "")
            if reac_id in sub_network_ids:
                global_tcosa_z_var_sum_expression += global_z_var
    base_problem += global_z_var_sum == global_z_var_sum_expression

    global_tcosa_z_var_sum = pulp.LpVariable(
        name="global_tcosa_z_var_sum",
        lowBound=-float("inf"),
        upBound=float("inf"),
    )
    base_problem += global_tcosa_z_var_sum == global_tcosa_z_var_sum_expression

    return base_problem
