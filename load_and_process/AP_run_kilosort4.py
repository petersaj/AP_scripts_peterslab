from kilosort.utils import download_probes
from kilosort import run_kilosort
import numpy as np

# Download probe channel map library
# (unused now that probe dictionary is always created below)
download_probes()

# Print pykilosort starting notice
print('PYTHON: Running kilosort 4 ', data_filename)

# Set settings
# NOTE: 'n_chan_bin' is a required setting, and should reflect the total number
#       of channels in the binary file. For information on other available
#       settings, see `kilosort.run_kilosort.default_settings`.
settings = {'filename': data_filename, 'n_chan_bin': n_chan}

# Create probe dictionary for recorded probe geometry
probe_info = {
    'chanMap': np.array(chanMap,dtype=int),
    'xc': np.array(xc),
    'yc': np.array(yc),
    'kcoords': np.array(kcoords),
    'n_chan': n_chan
}

# Run Kilsort
ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = \
    run_kilosort(settings=settings, results_dir=kilosort_output_path ,probe=probe_info)

# Print pykilosort ending notice
print('PYTHON: Finished kilosort 4: ', kilosort_output_path)
