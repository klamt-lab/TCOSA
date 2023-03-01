import cobra
from typing import List

def get_all_tcosa_reaction_ids(cobra_model: cobra.Model) -> List[str]:
    return [x.id for x in cobra_model.reactions if x.id.endswith("_TCOSA")] + ["NADK", "NADTRHD", "NADPPPS", "THD2pp"]
