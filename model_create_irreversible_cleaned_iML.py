#!/usr/bin/env python3
#
# Copyright 2020-2021 PSB
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""This scripts loads and cleans the iML1515 SBML (Monk et al., 2017) as given in the BiGG database (Schellenberger et al., 2010) and saves it again with cobrapy.

This is done in order to gain a version which is not altered by cobrapy while loading it.
The "cleaning" contains steps such as setting 1000 flux bounds to inf.

References:
* Monk, J. M., Lloyd, C. J., Brunk, E., Mih, N., Sastry, A., King, Z., ... & Feist, A. M. (2017).
  iML1515, a knowledgebase that computes Escherichia coli traits.
  Nature biotechnology, 35(10), 904-908.
* Schellenberger, J., Park, J. O., Conrad, T. M., & Palsson, B. Ø. (2010).
  BiGG: a Biochemical Genetic and Genomic knowledgebase of large scale metabolic
  reconstructions. BMC bioinformatics, 11(1), 213.
"""
# IMPORT SECTION
import cobra
from make_model_irreversible import make_model_irreversible


# ACTUAL ROUTINE SECTION
print("=>Loading and saving of original iML1515 model using cobrapy and cleaning of wrong metabolite and reaction ID parts")

print("Load original iML1515 SBML...")
model = cobra.io.read_sbml_model("./resources/iML1515.xml")

print("Set -inf/+inf bounds to -1000/1000...")
for reaction in model.reactions:
    if reaction.upper_bound >= 1000:
        reaction.upper_bound = 1000
    if reaction.lower_bound <= -1000:
        reaction.lower_bound = -1000

print("Cleaning reaction ID parts...")
for reaction in model.reactions:
    if "_DASH_" in reaction.id:
        reaction.id = reaction.id.replace("_DASH_", "__")
    if "_LPAREN_" in reaction.id:
        reaction.id = reaction.id.replace("_LPAREN_", "_")
    if "_RPAREN_" in reaction.id:
        reaction.id = reaction.id.replace("_RPAREN_", "")

print("Delete boundary-condition-free exchange metabolites...")
metabolite_ids = [x.id for x in model.metabolites]
for metabolite_id in metabolite_ids:
    if metabolite_id.endswith("_ex"):
        metabolite_instance = model.metabolites.get_by_id(metabolite_id)
        model.remove_metabolites([metabolite_instance])

print("Remove exchange metabolite reactions...")
reaction_ids = [x.id for x in model.reactions]
for reaction_id in reaction_ids:
    reaction = model.reactions.get_by_id(reaction_id)
    if reaction.metabolites == {}:
        model.remove_reactions([reaction])
        continue

print("Rename wrong metabolite name parts...")
for metabolite in model.metabolites:
    if "_DASH_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_DASH_", "__")
    if "_LPAREN_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_LPAREN_", "_")
    if "_RPAREN_" in metabolite.id:
        metabolite.id = metabolite.id.replace("_RPAREN_", "")

print("Deactivate all C sources except of D-glucose...")
model.reactions.EX_succ_e.lower_bound = 0
model.reactions.EX_glyc_e.lower_bound = 0
model.reactions.EX_ac_e.lower_bound = 0

print("Make iML1515 irreversible...")
model = make_model_irreversible(model)

print("Save iML1515 irreversible model version...")
cobra.io.write_sbml_model(model, "./resources/iML1515_irreversible_cleaned.xml")
