from kilosort.utils import download_probes
from kilosort import run_kilosort

# Download probe channel map library
download_probes()

# Print pykilosort starting notice
print('PYTHON: Running kilosort 4 ', data_filename)

# Set settings
# NOTE: 'n_chan_bin' is a required setting, and should reflect the total number
#       of channels in the binary file. For information on other available
#       settings, see `kilosort.run_kilosort.default_settings`.
settings = {'filename': data_filename, 'n_chan_bin': 384}

ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate = \
    run_kilosort(settings=settings, results_dir=kilosort_output_path ,probe_name='neuropixPhase3B2_kilosortChanMap.mat')

# Print pykilosort ending notice
print('PYTHON: Finished kilosort 4: ', kilosort_output_path)