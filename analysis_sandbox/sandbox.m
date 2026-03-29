%% SANDBOX
%
% dev/test code


%% MCMS API test

mcms_token = plab.mcms.login;

weights = plab.mcms.query(mcms_token,'weight','AP032');



%% Find all animals that have done specific task
% (this takes a long time because it looks through all server folders)

% workflow = '*opacity*';
% workflow = '*no_change*';
% workflow = '*opacity*';
% workflow = '*angle*';
% workflow = 'stim_wheel_right_stage1';
workflow = 'ImageDisplay';

task_dir = dir(fullfile(plab.locations.server_data_path, ...
    '**','bonsai','**',strcat(workflow,'.bonsai')));

animals = unique(extractBetween({task_dir.folder}, ...
    [plab.locations.server_data_path,filesep], ...
    filesep));




%% 4-shank processing test

oe_metadata = jsondecode(fileread("D:\4shank_test\experiment1\recording1\structure.oebin"));

n_channels = oe_metadata.continuous(1).num_channels;

channel_info = horzcat(oe_metadata.continuous(1).channels.channel_metadata);

electrode_idx = vertcat(channel_info(2,:).value);

% need to make probe dictionary? 
% https://kilosort.readthedocs.io/en/latest/tutorials/make_probe.html

% X/Y pos are in settings.xml


curr_ephys_settings_file = "D:\4shank_test\settings.xml";
curr_ephys_settings = readstruct(curr_ephys_settings_file,'filetype','xml');

% (wow, this is really buried...)
channels_str = curr_ephys_settings.SIGNALCHAIN.PROCESSOR(1).EDITOR.CUSTOM_PARAMETERS.NP_PROBE.CHANNELS;
x_str = curr_ephys_settings.SIGNALCHAIN.PROCESSOR(1).EDITOR.CUSTOM_PARAMETERS.NP_PROBE.ELECTRODE_XPOS;
y_str = curr_ephys_settings.SIGNALCHAIN.PROCESSOR(1).EDITOR.CUSTOM_PARAMETERS.NP_PROBE.ELECTRODE_YPOS;

channel_idx = cellfun(@(x) sscanf(x,'CH%dAttribute'),fieldnames(x_str));

x_pos = struct2array(x_str);
y_pos = struct2array(y_str);

shank = cellfun(@(x) sscanf(x,'%*d:%d'),struct2array(channels_str));

% KS probe dictionary
% these need to be numpy ndarrays


chanMap = (1:length(channel_idx))-1;
xc = x_pos;
yc = y_pos;
k_coords = shank;
n_chan = length(channel_idx);


% probe = {
%     'chanMap': chanMap,
%     'xc': xc,
%     'yc': yc,
%     'kcoords': kcoords,
%     'n_chan': n_chan
% }

% In kilosort, will need "probe" argument
% is this all available for 1-shank? if so, just do this every time?

% p = load_probe('.../test_prb.prb')
% results = run_kilosort(..., probe=p)

curr_ephys_settings_file = "\\qnap-ap001.dpag.ox.ac.uk\APlab\Data\PG003\2026-03-26\ephys\settings.xml";
curr_ephys_settings = readstruct(curr_ephys_settings_file,'filetype','xml');



oe_settings_fn = "D:\4shank_test\settings.xml";
oe_settings_fn = "\\qnap-ap001.dpag.ox.ac.uk\APlab\Data\PG003\2026-03-26\ephys\settings.xml";

probe = plab.ephys.oe_probe_geometry(oe_settings_fn);






