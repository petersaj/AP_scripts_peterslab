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

figure; tiledlayout(1,3); 
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



%% ~~~~~~~~~~~ BATCH WIDEFIELD

%% Batch widefield (master U) - passive V/A 

% Loop through animals
animals = {'AP021','AP022','DS001', ... & V>A
    'DS000','DS003','DS004'};           % A>V

% Set times for PSTH        
raster_window = [-0.5,1];
psth_bin_size = 0.03;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

rec_days_all = cell(length(animals),1);
task_workflow_all = cell(length(animals),1);
day_V_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
   
    animal = animals{curr_animal};

    vis_workflow = 'lcr_passive';
    aud_workflow = 'hml_passive_audio';

    vis_recordings = plab.find_recordings(animal,[],vis_workflow);
    aud_recordings = plab.find_recordings(animal,[],aud_workflow);

    rec_days = unique(horzcat(...
        {vis_recordings(cellfun(@any,{vis_recordings.widefield})).day},...
        {aud_recordings(cellfun(@any,{aud_recordings.widefield})).day}));

    task_workflow = cell(length(rec_days),1);
    day_V = cell(length(rec_days),2);

    for curr_day = 1:length(rec_days)

        % Get the task workflow
        task_workflow_name = 'stim_wheel*';
        day_recordings = plab.find_recordings(animal,rec_days{curr_day},task_workflow_name);
        if ~isempty(day_recordings)
            task_workflow{curr_day} = day_recordings.workflow{end};
        end

        for curr_modality = 1:2
            switch curr_modality
                case 1
                    curr_workflow = vis_workflow;
                case 2
                    curr_workflow = aud_workflow;
            end

            % Grab pre-load vars
            preload_vars = who;

            % Load data
            try
                curr_recording = plab.find_recordings(animal,rec_days{curr_day},curr_workflow);
            catch me
                continue
            end

            rec_day = curr_recording.day;
            rec_time = curr_recording.recording{end};

            load_parts = struct;
            load_parts.widefield = true;
            load_parts.widefield_master = true;
            ap.load_recording;

            % If widefield isn't aligned, skip
            if ~load_parts.widefield_align
                continue
            end

            % Stim times (quiescent trials only)
            switch curr_modality
                case 1
                    stim_type = vertcat(trial_events.values.TrialStimX);             
                case 2
                    stim_type = vertcat(trial_events.values.StimFrequence);
            end
            
            stim_window = [0,0.5];
            quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
                timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
                timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
                1:length(stimOn_times))';
      
            % sometimes not all trials shown?
            stim_type = stim_type(1:length(stimOn_times));
            align_times = cellfun(@(x) stimOn_times(stim_type == x & quiescent_trials), ...
                num2cell(unique(stim_type)),'uni',false);

            % Align ROI trace to align times
            aligned_V_mean = nan(size(wf_V,1),length(t_centers),length(align_times));
            for curr_align = 1:length(align_times)

                % (skip if <5 usable trials)
                if length(align_times{curr_align}) < 5
                    continue
                end

                peri_event_t = align_times{curr_align} + t_centers;

                aligned_V = interp1(wf_t,wf_V',peri_event_t);

                aligned_V_baselinesub = aligned_V - ...
                    mean(aligned_V(:,t_centers < 0,:),2);

                aligned_V_mean(:,:,curr_align) = ...
                    permute(nanmean(aligned_V_baselinesub,1),[3,2,1]);
            end

            % Store mean aligned ROI trace
            day_V{curr_day,curr_modality} = aligned_V_mean;

            % Prep for next loop
            clearvars('-except',preload_vars{:});

        end

        ap.print_progress_fraction(curr_day,length(rec_days));
    end

    rec_days_all{curr_animal} = rec_days;
    task_workflow_all{curr_animal} = task_workflow;
    day_V_all{curr_animal} = day_V;
    disp(['Done ' animal]);

end


% Load master U, convert V to px
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);

% Average and plot
use_modality = 1;
use_stim = 3;
% use_workflow = 'stim_wheel_right_stage2_audio_volume';
use_workflow = 'stim_wheel_right_stage2';
% use_workflow = 'stim_wheel_right_stage2_mixed_VA';

x = cellfun(@(x,task) cat(4,x{strcmp(task,use_workflow),use_modality}),day_V_all,task_workflow_all,'uni',false);
x2 = cellfun(@(x) plab.wf.svd2px(U_master,squeeze(x(:,:,use_stim,:))),x,'uni',false);
x3 = cellfun(@(x) nanmean(x,4),x2,'uni',false);
x4 = nanmean(cat(4,x3{:}),4);

ap.imscroll(x4);
axis image;
colormap(AP_colormap('PWG',[],1));
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');

%% Batch fraction of responsive units 

% Loop through animals
animals = {'DS004','DS007'};

frac_responsive_units = cell(length(animals),1);
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};

    recordings = plab.find_recordings(animal);
    ephys_recordings = recordings([recordings.ephys]);

    frac_responsive_units{curr_animal} = nan(length(ephys_recordings),2);

    for curr_day = 1:length(ephys_recordings)
        for curr_modality = 1:2
            switch curr_modality
                case 1
                    curr_workflow = 'lcr_passive';
                case 2
                    curr_workflow = 'hml_passive_audio';
            end

            % Grab pre-load vars
            preload_vars = who;

            % Load data
            curr_recording = plab.find_recordings(animal,ephys_recordings(curr_day).day,curr_workflow);
    
            rec_day = curr_recording.day;
            rec_time = curr_recording.recording{end};
            ap.load_recording;

            % Set event to get response
            % (get quiescent trials)
            stim_window = [0,0.5];
            quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
                timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
                timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
                (1:length(stimOn_times))');

            if contains(bonsai_workflow,'lcr')
                % (vis passive)
                stim_type = vertcat(trial_events.values.TrialStimX);
                use_align = stimOn_times(stim_type(1:length(stimOn_times)) == 90 & quiescent_trials);
            elseif contains(bonsai_workflow,'hml')
                % (aud passive)
                stim_type = vertcat(trial_events.values.StimFrequence);
                use_align = stimOn_times(stim_type == 8000 & quiescent_trials);
            elseif contains(bonsai_workflow,'stim_wheel')
                % (task)
                use_align = stimOn_times(stim_to_move > 0.15);
            end

            baseline_t = [-0.2,0];
            response_t = [0,0.2];

            baseline_bins = use_align + baseline_t;
            response_bins = use_align + response_t;

            event_bins = [baseline_bins,response_bins];
            spikes_binned_continuous = histcounts2(spike_times_timelite,spike_templates, ...
                reshape([baseline_bins,response_bins]',[],1),1:size(templates,1)+1);

            event_spikes = permute(reshape(spikes_binned_continuous(1:2:end,:),2, ...
                size(event_bins,1),[]),[2,1,3]);

            event_response = squeeze(mean(diff(event_spikes,[],2),1));

            n_shuff = 1000;
            event_response_shuff = cell2mat(arrayfun(@(shuff) ...
                squeeze(mean(diff(ap.shake(event_spikes,2),[],2),1)), ...
                1:n_shuff,'uni',false));

            event_response_rank = tiedrank(horzcat(event_response,event_response_shuff)')';
            event_response_p = event_response_rank(:,1)./(n_shuff+1);
            responsive_units = event_response_p < 0.05 | event_response_p > 0.95;

            frac_responsive_units{curr_animal}(curr_day,curr_modality) = ...
                mean(responsive_units);

        end
    end

    ap.print_progress_fraction(curr_animal,length(animals));

end



%% ~~~~~~~~~~~ BATCH EPHYS

%% Plot all recording locations
animals = {'AP022','DS000','DS007','DS010', ...
    'AP019','AP021','DS001','DS003','DS004'};

animal_col = [brewermap(8,'Dark2');brewermap(8,'Set2')];
ccf_draw = ap.ccf_draw;

for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    probe_color = animal_col(curr_animal,:);
    ccf_draw.draw_probes_nte;
    

    % Trajectory explorer probe position files
    nte_filepattern = plab.locations.filename('server',animal,'*',[],'ephys','*probe_positions.mat');
    nte_fns = dir(nte_filepattern);

    for curr_recording = 1:length(nte_fns)

        load(fullfile(nte_fns(curr_recording).folder, nte_fns(curr_recording).name));

        % Loop through probes and draw
        for curr_probe = 1:length(probe_positions_ccf)

            % Probe line is directly stored
            probe_line = probe_positions_ccf{curr_probe}';

            % Draw probe in 3D view
            line(ccf_3d_axes,probe_line(:,1),probe_line(:,3),probe_line(:,2), ...
                'linewidth',2,'color',probe_color)

            % Draw probes on coronal + saggital
            line(ccf_axes(1),probe_line(:,3),probe_line(:,2),'linewidth',2,'color',probe_color);
            line(ccf_axes(2),probe_line(:,3),probe_line(:,1),'linewidth',2,'color',probe_color);
            line(ccf_axes(3),probe_line(:,2),probe_line(:,1),'linewidth',2,'color',probe_color);

            % Draw probe start/end on horizontal
            plot(ccf_axes(2), probe_line(1,3),probe_line(1,1), ...
                'o','MarkerSize',5,'color',probe_color);
            plot(ccf_axes(2), probe_line(end,3),probe_line(end,1), ...
                '.','MarkerSize',20,'color',probe_color);

        end
    end
end

%% Do something

animals = {'AP022','DS000','DS007','DS010', ...
    'AP019','AP021','DS001','DS003','DS004'};

for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    workflow = 'stim_wheel_right_stage2_mixed_VA';
    recordings = plab.find_recordings(animal,[],workflow);

    recordings_ephys = recordings([recordings.ephys]);


end


















