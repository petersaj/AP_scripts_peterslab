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




%% Move pupil SLEAP files

sleap_dir = dir(fullfile(plab.locations.server_path,'Users','Peter_Gorman','Pupils_temporary','Outputs_Sorted','lcr_passive','**','*.h5'));
sleap_fileparts = regexp({sleap_dir.name},'(?<animal>.*?)_(?<rec_day>.*?)_(?<rec_time>Recording_.*?)_mousecam','names');

for curr_file_idx = 1:length(sleap_dir)  
    curr_fn = fullfile(sleap_dir(curr_file_idx).folder,sleap_dir(curr_file_idx).name);
    target_fn = fullfile(plab.locations.server_data_path, ...
        sleap_fileparts{curr_file_idx}.animal, ...
        sleap_fileparts{curr_file_idx}.rec_day, ...
        sleap_fileparts{curr_file_idx}.rec_time, ...
        'mousecam','sleap','pupil_v1',sleap_dir(curr_file_idx).name);

    mkdir(fileparts(target_fn));
    copyfile(curr_fn,target_fn);
    fprintf('Copied: %s\n',sleap_dir(curr_file_idx).name)
end

%%

% timestamps: convert to timestamps + time jump


plab.ephys.ks2oe_timestamps(ks_spike_times_fn,oe_samples_fn, ...
    oe_metadata(1).continuous(1).sample_rate);

spike_times_openephys = readNPY(spike_times_openephys_filename);


open_ephys_ttl_timestamps = ...
    double(vertcat(open_ephys_ttl_sample_numbers{:}))/oe_ap_samplerate;


% Create array of timestamps

% maybe do this in plab.ephys.ks2oe_timestamps?
% the problem is that the ttl timestamps will be wrong though

oe_bad_samples = plab.ephys.find_dropped_ephys(open_ephys_path{1});





