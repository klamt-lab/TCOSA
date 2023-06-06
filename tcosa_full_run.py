"""
This script re-runs all calculations performed for the TCOSA publication.
In order to really run all calculations again, delete the "cosa" subfolder
beforehand as all OptMDFpathway-based calculations for the random sampling
are cached and would not be re-run.
All scripts run for aerobic and anaerobic conditions and can, where useful,
use the expanded model.
These "imports" will directly run the top-level code of the modules
thereby starting the calculations.
"""

# Create an irreversible and cleaned-up version of iML which
# will be later converted to a TCOSA model.
import model_create_irreversible_cleaned_iML

# Calculate the ΔG'° values using the eQuilibrator
import model_get_iML_dG0_data

# Get ΔG'° statistics (these are implicitly used later already) for ...
# ...all reactions
import model_dG0_statistics
# ...NAD(P)(H) reactions only
import cosa_dG0_TCOSA_reaction_statistics

# Transform the concentration data form Bennett et al., 2009
# into a machine-readable format.
import model_in_vivo_concentration_data_setup

# Create the actual TCOSA models, including iML1515_TCOSA.
import cosa_create_model_and_data

# Perform the random specificity sampling where the in vivo distribution
# is compared with random distributions and other special distributions. Herein, the CSV tables
# and the figures for this analysis are also generated.
# >Aerobic with glucose, including 3-cofactor analyses
import cosa_random_sampling_aerobic
# >Anaerobic with glucose, including 3-cofactor analyses
import cosa_random_sampling_anaerobic
# >Aerobic with acetate
import cosa_random_sampling_acetate
# >Create figures
import cosa_random_sampling_create_figures

# Create the extended vs. two-cofactor model comparison figure.
import cosa_create_extended_model_figure

# Run the ratio ratio range variability analysis and generate its figures.
import cosa_ratio_ratio_test

# Create in vivo concentrations and acetate figures
import cosa_create_special_aerobic_figures

# Run the number of minimal changes to reach theoretical optimality analysis
# and generate a text report about it.
import cosa_minimal_changes_test

# Perform dGf or redox potential sampling (called dG0 sampling for historical reasons).
import cosa_dG0_sampling_aerobic
import cosa_dG0_sampling_anaerobic

# Single swap analyses
import cosa_single_swap_effect_analysis_aerobic
import cosa_single_swap_effect_analysis_anaerobic
import cosa_single_swap_effect_analysis_statistics

# Thermodynamic Flux Variability Analyses to find out minimal and maximal active reactions
import cosa_fva_aerobic
import cosa_fva_anaerobic
import cosa_fva_statistics

# Perform stoichiometric FVAs
import cosa_stoichiometric_fva_aerobic
import cosa_stoichiometric_fva_anaerobic

# Test full switch of NAD and NADP in a source-code selected solution
import cosa_switch_test

# Perform Concentration Variability Analyses
import cosa_cva

# Create Supplementary Table 3
import cosa_fva_excel_comparison_files
