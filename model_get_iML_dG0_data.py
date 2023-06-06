"""This script finds the dG0 values using the equilibrator-API wrapper from the equilibrator submodule.

Here, all relevant thermodynamic settings are also given.
"""

# IMPORT SECTION #
# External imports
import cobra
from typing import Dict, List, Tuple
# Internal imports
from equilibrator import get_model_dG0_values
from helper import json_write

## User variables
inner_to_outer_compartments: List[str] = [
    "c",
    "p",
    "e",
]
phs: Dict[str, float] = {  # Unitless
    "c": 7.5,  # Source: Bionumbers ID 105980
    "p": 7.0,
    "e": 7.0,
}
pmgs: Dict[str, float] = {  # Unitless
    "c": 2.5,  # Source: eQuilibrator standard
    "p": 2.5,
    "e": 2.5,
}
ionic_strengths: Dict[str, float] = {  # In mM
    "c": 250,  # Source: eQuilibrator standard
    "p": 250,
    "e": 250,
}
potential_differences: Dict[Tuple[str, str], float] = {  # In V
    ("c", "p"): -0.15,  # Source: eQuilibrator standard
    ("p", "e"): -0.15,
}

# Load iML1515 and perform the dG0 determination
cobra_model = cobra.io.read_sbml_model("resources/iML1515_irreversible_cleaned.xml")
dG0_values = get_model_dG0_values(
    cobra_model,
    inner_to_outer_compartments,
    phs,
    pmgs,
    ionic_strengths,
    potential_differences
)
json_write("./resources/dG0_iML1515_irreversible_cleaned.json", dG0_values)
