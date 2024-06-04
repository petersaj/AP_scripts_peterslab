%% Batch widefield (master U)

% Loop through animals
animals = {'AP021','AP022'};

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.03;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

day_V_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
   
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
%     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    n_components = 2000;
    day_V = nan(n_components,length(t_centers),length(recordings));

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % If widefield isn't aligned, skip
        if ~load_parts.widefield_align
            continue
        end

        % Get quiescent trials
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        if isfield(trial_events.values,'TrialStimX')
            stim_x = vertcat(trial_events.values.TrialStimX);
            % temp - what happened here, not all trials shown?
            stim_x = stim_x(1:length(stimOn_times));
            align_times = stimOn_times(stim_x == 90 & quiescent_trials);
        else
            % Stim times (task)
            align_times = stimOn_times;
        end

        % Align ROI trace to align times    
        peri_event_t = align_times + t_centers;
        
        aligned_V = interp1(wf_t,wf_V',peri_event_t);

        aligned_V_baselinesub = aligned_V - ...
            mean(aligned_V(:,t_centers < 0,:),2);

        % Store mean aligned ROI trace
        day_V(:,:,curr_recording) = permute(nanmean(aligned_V_baselinesub,1),[3,2,1]);

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_V_all{curr_animal} = day_V;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Load master U
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');

load(master_U_fn);









