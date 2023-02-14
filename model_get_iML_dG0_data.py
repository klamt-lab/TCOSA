import cobra
from equilibrator import get_model_dG0_values
from helper import json_write
from typing import Dict, List, Tuple

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
