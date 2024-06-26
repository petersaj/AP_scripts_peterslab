%% Behavior (quick)

trial_modality = [trial_events.values(1:n_trials).TaskType]';

% Get association significance
trial_events_vis = trial_events;
trial_events_vis.timestamps = trial_events_vis.timestamps(trial_modality == 0);
trial_events_vis.values = trial_events_vis.values(trial_modality == 0);

trial_events_aud = trial_events;
trial_events_aud.timestamps = trial_events_aud.timestamps(trial_modality == 1);
trial_events_aud.values = trial_events_aud.values(trial_modality == 1);

vis_p = AP_stimwheel_association_pvalue(stimOn_times(trial_modality == 0),trial_events_vis,stim_to_move(trial_modality == 0));
aud_p = AP_stimwheel_association_pvalue(stimOn_times(trial_modality == 1),trial_events_aud,stim_to_move(trial_modality == 1));

% Align wheel movement to stim onset
surround_time = [-5,5];
surround_sample_rate = 100;
surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

align_times = stimOn_times(1:n_trials);
pull_times = align_times + surround_time_points;

event_aligned_wheel_vel = interp1(timelite.timestamps, ...
    wheel_velocity,pull_times);
event_aligned_wheel_move = interp1(timelite.timestamps, ...
    +wheel_move,pull_times,'previous');

% Plot wheel move average by modality
stim_wheel_vel_med = ap.groupfun(@median,event_aligned_wheel_vel,trial_modality,[]);
stim_wheel_move_avg = ap.groupfun(@mean,event_aligned_wheel_move,trial_modality,[]);

figure; tiledlayout('flow'); 
nexttile; hold on
plot(surround_time_points,stim_wheel_vel_med','linewidth',2);
ylabel('Median velocity');

nexttile; hold on
plot(surround_time_points,stim_wheel_move_avg','linewidth',2);
ylabel('Frac moving')
legend({'Vis','Aud'},'location','nw');
title(sprintf('Vis p = %.2f\nAud p = %.2f',vis_p,aud_p));

nexttile; hold on;
histogram(stim_to_move(trial_modality == 0),[0:0.01:0.5]);
histogram(stim_to_move(trial_modality == 1),[0:0.01:0.5]);
xlabel('Reaction time');
ylabel('Count');



%% PSTH - units (split by visual/auditory)

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% PSTH for all conditions
if contains(bonsai_workflow,'lcr')
    % (visual passive)
    stim_type = vertcat(trial_events.values.TrialStimX);
    align_times = cellfun(@(x) stimOn_times(stim_type == x),num2cell(unique(stim_type)),'uni',false);
elseif contains(bonsai_workflow,'hml')
    % (auditory passive)
    stim_type = vertcat(trial_events.values.StimFrequence);
    align_times = cellfun(@(x) stimOn_times(stim_type == x),num2cell(unique(stim_type)),'uni',false);
elseif contains(bonsai_workflow,'stim_wheel')
    % (task)
    trial_modality = [trial_events.values(1:n_trials).TaskType];
    align_times = { ...
        stimOn_times(trial_modality==0),stimOn_times(trial_modality==1), ...
        stim_move_time,reward_times};
end

n_units = size(templates,1);
unit_psth = nan(n_units,length(t_bins)-1,2);
for curr_align = 1:length(align_times)
    use_align = align_times{curr_align};
    t_peri_event = bsxfun(@plus,use_align,t_bins);
    for curr_unit = 1:n_units
        curr_spikes = spike_times_timelite(spike_templates == curr_unit);

        curr_spikes_binned = cell2mat(arrayfun(@(x) ...
            histcounts(curr_spikes,t_peri_event(x,:)), ...
            [1:size(t_peri_event,1)]','uni',false))./psth_bin_size;
        curr_mean_psth = mean(curr_spikes_binned,1);

        unit_psth(curr_unit,:,curr_align) = curr_mean_psth;
    end
end

smooth_size = 100;
unit_psth_smooth = smoothdata(unit_psth,2,'gaussian',smooth_size);

% Normalize to baseline
unit_baseline = nanmean(nanmean(unit_psth(:,t_bins(2:end) < 0,:),2),3);
unit_psth_smooth_norm = (unit_psth_smooth-unit_baseline)./(unit_baseline+1);

% Plot depth-sorted
[~,sort_idx] = sort(template_depths);
AP_imscroll(unit_psth_smooth_norm(sort_idx,:,:));
clim([-2,2]);
colormap(AP_colormap('BWR'));



%% ~~~~~~~~~~~ BATCH 

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
            stim_type = vertcat(trial_events.values.TrialStimX);
            % temp - what happened here, not all trials shown?
            stim_type = stim_type(1:length(stimOn_times));
            align_times = stimOn_times(stim_type == 90 & quiescent_trials);
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









