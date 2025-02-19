"""Performs the dGf sampling under anaerobic conditions."""
from cosa_dG0_sampling import cosa_dG0_sampling
from cosa_random_sampling_figures import create_cosa_dG0_sampling_figures, create_total_dG0_sampling_figure

cosa_dG0_sampling(anaerobic=True, expanded=False, num_samplings=100, change_range=5)
create_total_dG0_sampling_figure(change_range=5)
