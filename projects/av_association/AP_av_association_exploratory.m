%% Mixed-task behavior

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

%% Grab average passive widefield for one animal

animal = 'AP018';

% task_workflow = 'stim_wheel_right_stage\d_size_up';
% task_workflow = 'stim_wheel_right_stage\d_angle_size60';
task_workflow = 'stim_wheel_right_stage\d_audio_volume';
% task_workflow = 'stim_wheel_right_stage\d';

% passive_workflow = 'lcr_passive';
passive_workflow = 'hml_passive_audio';

% Set times for PSTH        
raster_window = [-0.5,1];
psth_bin_size = 0.03;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

task_recordings = plab.find_recordings(animal,[],task_workflow);
passive_recordings = plab.find_recordings(animal,[],passive_workflow);



vis_workflow = 'lcr_passive';
aud_workflow = 'hml_passive_audio';
task_workflow = 'stim_wheel*';

vis_recordings = plab.find_recordings(animal,[],vis_workflow);
aud_recordings = plab.find_recordings(animal,[],aud_workflow);
task_recordings = plab.find_recordings(animal,[],task_workflow);

rec_days = {task_recordings(cellfun(@any,{task_recordings.widefield})).day};

task_workflow = cell(length(rec_days),1);
day_V = cell(length(rec_days),2);

for curr_day = 1:length(rec_days)

    % Get the task workflow
    task_workflow_name = 'stim_wheel*';
    day_recordings = plab.find_recordings(animal,rec_days{curr_day},task_workflow_name);
    task_workflow{curr_day} = day_recordings.workflow{end};

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
        curr_recording = plab.find_recordings(animal,rec_days{curr_day},curr_workflow);
        if isempty(curr_recording)
            continue
        end

        rec_day = curr_recording.day;
        rec_time = curr_recording.recording{end};

        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;

        % (Try/catch: rare catastrophic dropped frames)
        try
            ap.load_recording;
        catch me
            continue
        end

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





%% ~~~~~~~~~~~ BATCH WIDEFIELD

%% Batch widefield (master U) - passive V/A 

% Loop through animals
animals = {'AP021','AP022','DS001','DS007','DS010','DS011', ... & V>A
    'DS000','DS003','DS004','DS013','DS014','DS015','DS016'};           % A>V

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
    task_workflow = 'stim_wheel*';

    vis_recordings = plab.find_recordings(animal,[],vis_workflow);
    aud_recordings = plab.find_recordings(animal,[],aud_workflow);
    task_recordings = plab.find_recordings(animal,[],task_workflow);

    rec_days = {task_recordings(cellfun(@any,{task_recordings.widefield})).day};

    task_workflow = cell(length(rec_days),1);
    day_V = cell(length(rec_days),2);

    for curr_day = 1:length(rec_days)

        % Get the task workflow
        task_workflow_name = 'stim_wheel*';
        day_recordings = plab.find_recordings(animal,rec_days{curr_day},task_workflow_name);
        task_workflow{curr_day} = day_recordings.workflow{end};

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
            curr_recording = plab.find_recordings(animal,rec_days{curr_day},curr_workflow);
            if isempty(curr_recording)
                continue
            end

            rec_day = curr_recording.day;
            rec_time = curr_recording.recording{end};

            load_parts = struct;
            load_parts.widefield = true;
            load_parts.widefield_master = true;
            
            % (Try/catch: rare catastrophic dropped frames)
            try
                ap.load_recording;
            catch me
                continue
            end

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

            % Align V to align times
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

% Save
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data';
save_fn = fullfile(data_path,'widefield_passive');
save(save_fn,'rec_days_all','task_workflow_all','day_V_all')
fprintf('Saved %s\n',save_fn);


% Load master U, convert V to px
U_master = plab.wf.load_master_U;

% Average and plot
use_modality = 2;
use_stim = 2;
% use_workflow = 'stim_wheel_right_stage1';
% use_workflow = 'stim_wheel_right_stage2';
use_workflow = 'stim_wheel_right_stage1_audio_volume';
% use_workflow = 'stim_wheel_right_stage2_audio_volume';
% use_workflow = 'stim_wheel_right_stage2_mixed_VA';

% use_animals = 1:6; % VA
use_animals = 7:13; % AV
x = cellfun(@(x,task) cat(4,x{strcmp(task,use_workflow),use_modality}), ...
    day_V_all(use_animals),task_workflow_all(use_animals),'uni',false);
x2 = cellfun(@(x) plab.wf.svd2px(U_master,squeeze(x(:,:,use_stim,:))),x,'uni',false);
x3 = cellfun(@(x) nanmean(x,4),x2,'uni',false);
x4 = nanmean(cat(4,x3{:}),4);

ap.imscroll(x4,t_centers);
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

animals = {'AP022','DS007','DS010','DS011', ... % V-A, performed both
    'DS000','DS013','DS014','DS015','DS016', ... % A-V, performed both
    'DS005', ... % A-V, frequency task
    'AP021','DS001', ... % V-A, didn't perform A mixed
    'DS003','DS004', ... % A-v, didn't perform A mixed
    };

animal_col = [brewermap(8,'Dark2');brewermap(8,'Set2')];
ccf_draw = ap.ccf_draw;
ccf_draw.draw_name('Caudoputamen');

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    probe_color = animal_col(curr_animal,:);
    ccf_draw.draw_probes_nte(animal,probe_color);
end

%% Get recording locations/days

animals = {'AP022','DS007','DS010','DS011', ... % V-A, performed both
    'DS000','DS013','DS014','DS015','DS016', ... % A-V, performed both
    'DS005', ... % A-V, frequency task
    'AP021','DS001', ... % V-A, didn't perform A mixed
    'DS003','DS004', ... % A-v, didn't perform A mixed
    };

% Find anterior/posterior striatum recordings (from NTE)
str_rec_days = cell(length(animals),2); % [ant,pos]
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};

    nte_filepattern = plab.locations.filename('server',animal,'*',[],'ephys','*probe_positions.mat');
    nte_fns = dir(nte_filepattern);

    % Get days for all NTE files
    day_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
    rec_days = extract({nte_fns.folder},day_pattern);

    % Load NTE files
    nte_all = cell(size(nte_fns));
    for curr_recording = 1:length(nte_fns)
        nte_all{curr_recording} = load(fullfile(nte_fns(curr_recording).folder, ...
            nte_fns(curr_recording).name));
    end

    % Get NTE files with A/P striatum
    ap_striatum_threshold = 550;
    anterior_recordings = cellfun(@(x) x.probe_positions_ccf{1}(1,2),nte_all) < ap_striatum_threshold;
    striatum_recordings = cellfun(@(x) any(strcmp(x.probe_areas{1}.name,'Caudoputamen')),nte_all);

    str_ant_rec = rec_days(striatum_recordings & anterior_recordings);
    str_pos_rec = rec_days(striatum_recordings & ~anterior_recordings);

    str_rec_days{curr_animal,1} = str_ant_rec;
    str_rec_days{curr_animal,2} = str_pos_rec;

end

data_dir = 'C:\Users\petersa\Documents\PetersLab\analysis\ant_pos_striatum\data';
data_fn = fullfile(data_dir,'str_rec_days');
save(data_fn,'str_rec_days');

%% Get association performance

animals = {'AP022','DS007','DS010','DS011', ... % V-A, performed both
    'DS000','DS013','DS014','DS015','DS016', ... % A-V, performed both
    'DS005', ... % A-V, frequency task
    'AP021','DS001', ... % V-A, didn't perform A mixed
    'DS003','DS004', ... % A-v, didn't perform A mixed
    };

% Load recording days
data_dir = 'C:\Users\petersa\Documents\PetersLab\analysis\ant_pos_striatum\data';
load(fullfile(data_dir,'str_rec_days'));

% Find days with sensorimotor associations
bhv_p = cell(size(str_rec_days));
for curr_animal = 1:length(animals)
    for curr_ap = 1:2
        for curr_day = 1:length(str_rec_days{curr_animal,curr_ap})

            % Load data
            animal = animals{curr_animal};
            rec_day = str_rec_days{curr_animal,curr_ap}{curr_day};

            workflow = 'stim_wheel_right*mixed_VA';
            recordings = plab.find_recordings(animal,rec_day,workflow);
            rec_time = recordings.recording{end};

            load_parts.ephys = false;
            ap.load_recording;

            % Get association significance by modality
            trial_modality = [trial_events.values(1:n_trials).TaskType]';

            trial_events_vis = trial_events;
            trial_events_vis.timestamps = trial_events_vis.timestamps(trial_modality == 0);
            trial_events_vis.values = trial_events_vis.values(trial_modality == 0);

            trial_events_aud = trial_events;
            trial_events_aud.timestamps = trial_events_aud.timestamps(trial_modality == 1);
            trial_events_aud.values = trial_events_aud.values(trial_modality == 1);

            vis_p = AP_stimwheel_association_pvalue(stimOn_times(trial_modality == 0), ...
            trial_events_vis,stim_to_move(trial_modality == 0));
            aud_p = AP_stimwheel_association_pvalue(stimOn_times(trial_modality == 1), ...
                trial_events_aud,stim_to_move(trial_modality == 1));

            bhv_p{curr_animal,curr_ap}{curr_day} = [vis_p,aud_p];
        end
    end
    ap.print_progress_fraction(curr_animal,length(animals));
end

data_dir = 'C:\Users\petersa\Documents\PetersLab\analysis\ant_pos_striatum\data';
data_fn = fullfile(data_dir,'bhv_p');
save(data_fn,'bhv_p');


%% Striatum MUA (whole striatum together)

animals = {'AP022','DS007','DS010','DS011', ... % V-A, performed both
    'DS000','DS013','DS014','DS015','DS016', ... % A-V, performed both
    'DS005', ... % A-V, frequency task
    'AP021','DS001', ... % V-A, didn't perform A mixed
    'DS003','DS004', ... % A-v, didn't perform A mixed
    };

% Load recording days
data_dir = 'C:\Users\petersa\Documents\PetersLab\analysis\ant_pos_striatum\data';
load(fullfile(data_dir,'str_rec_days'));

% Get multiunit activity in striatum
striatum_mua_all = cell(length(animals),2);
for curr_animal = 1:length(animals)
    for curr_ap = 1:2
        for curr_day = 1:length(str_rec_days{curr_animal,curr_ap})

            % Load data
            animal = animals{curr_animal};
            rec_day = str_rec_days{curr_animal,curr_ap}{curr_day};

%             workflow = 'stim_wheel_right_stage2_mixed_VA';
            workflow = 'lcr_passive';
%             workflow = 'hml_passive_audio';
            recordings = plab.find_recordings(animal,rec_day,workflow);
            rec_time = recordings.recording{end};

            ap.load_recording;

            % Identify striatum spikes
            striatum_idx = strcmp(probe_areas{1}.name,'Caudoputamen');
            striatum_depth = prctile(probe_areas{1}.probe_depth(striatum_idx,:),[0,100],'all');
            striatum_spikes = spike_depths >= striatum_depth(1) & ...
                spike_depths <= striatum_depth(2);

            % Set times for PSTH
            raster_window = [-0.5,1];
            psth_bin_size = 0.001;
            t_bins = raster_window(1):psth_bin_size:raster_window(2);
            t_centers = conv2(t_bins,[1,1]/2,'valid');

            % Set align times:
            if contains(bonsai_workflow,'passive')
                % Passive
                % Get quiescent trials
                stim_window = [0,0.5];
                quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
                    timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
                    timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
                    (1:length(stimOn_times))');

                if isfield(trial_events.values,'TrialStimX')
                    align_category = vertcat(trial_events.values.TrialStimX);
                elseif isfield(trial_events.values,'StimFrequence')
                    align_category = vertcat(trial_events.values.StimFrequence);
                end
                align_times = cellfun(@(x) ...
                    stimOn_times(align_category(1:length(stimOn_times)) == x & ...
                    quiescent_trials),num2cell(unique(align_category)),'uni',false);

            elseif contains(bonsai_workflow,'stim_wheel')
                % Task
                trial_type = [trial_events.values.TaskType];
                align_times = cell(4,1);
                % Stim on (by type)
                align_times{1} = stimOn_times(trial_type(1:n_trials) == 0);
                align_times{2} = stimOn_times(trial_type(1:n_trials) == 1);
                % Stim move
                align_times{3} = stim_move_time;
                % Reward
                align_times{4} = reward_times_task;
%                 % ITI move
%                 wheel_starts = timelite.timestamps(diff([0;wheel_move]) == 1);
%                 wheel_stops = timelite.timestamps(diff([0;wheel_move]) == -1);
%                 iti_move_idx = interp1(photodiode_times, ...
%                     photodiode_values,wheel_starts,'previous') == 0;
%                 align_times{5} =  wheel_starts(iti_move_idx);
            end
           
            % Get multiunit PSTH
            striatum_psth = nan(length(align_times),length(t_bins)-1);
            for curr_align = 1:length(align_times)

                t_peri_event = align_times{curr_align} + t_bins;

                use_spikes_time = spike_times_timelite >= min(t_peri_event,[],'all') & ...
                    spike_times_timelite <= max(t_peri_event,[],'all');

                curr_use_spikes = striatum_spikes & use_spikes_time;

                spikes_binned_continuous = histcounts(spike_times_timelite(curr_use_spikes), ...
                    reshape(t_peri_event',[],1))./psth_bin_size;

                use_continuous_bins = reshape(padarray(true(size(t_peri_event(:,1:end-1)')),[1,0],false,'post'),[],1);
                spikes_binned = reshape(spikes_binned_continuous(use_continuous_bins),size(t_peri_event(:,1:end-1)'))';
                striatum_psth(curr_align,:) = nansum(spikes_binned,1);
            end

            striatum_mua_all{curr_animal,curr_ap}{curr_day} = striatum_psth;
        end
    end
    ap.print_progress_fraction(curr_animal,length(animals));
end

%%% ---> TODO: Save data



figure;
for str_ap = 1:2
    mua_avg = smoothdata( ...
        cell2mat(permute(cat(2,striatum_mua_all{:,str_ap}),[1,3,2])), ...
        2,'gaussian',100);
    mua_avg_norm = (mua_avg-nanmean(mua_avg(:,1:500,:),2))./nanmean(mua_avg(:,1:500,:),2);
    nexttile;plot(nanmean(mua_avg_norm,3)');
    switch str_ap
        case 1
            title('Anterior');
        case 2
            title('Posterior');
    end
end
linkaxes(gcf().Children.Children);



%% Striatum unit PSTHs: passive

animals = {'AP022','DS007','DS010','DS011', ... % V-A, performed both
    'DS000','DS013','DS014','DS015','DS016', ... % A-V, performed both
    'DS005', ... % A-V, frequency task
    'AP021','DS001', ... % V-A, didn't perform A mixed
    'DS003','DS004', ... % A-V, didn't perform A mixed
    };

% Load recording days
data_dir = 'C:\Users\petersa\Documents\PetersLab\analysis\ant_pos_striatum\data';
load(fullfile(data_dir,'str_rec_days'));

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Get PSTHs and responsive units
unit_psth_all = cell(length(animals),2);
responsive_units_all = cell(length(animals),2);
for curr_animal = 1:length(animals)
    for curr_ap = 1:2
        for curr_day = 1:length(str_rec_days{curr_animal,curr_ap})
            for curr_modality = 1:2

            % Load data
            animal = animals{curr_animal};
            rec_day = str_rec_days{curr_animal,curr_ap}{curr_day};

            switch curr_modality
                case 1
                    workflow = 'lcr_passive';
                case 2
                    workflow = 'hml_passive_audio';
            end
            recordings = plab.find_recordings(animal,rec_day,workflow);
            rec_time = recordings.recording{end};

            ap.load_recording;

            % Get quiescent trials
            stim_window = [0,0.5];
            quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
                timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
                timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
                (1:length(stimOn_times))');

            % Get stim times to use
            if isfield(trial_events.values,'TrialStimX')
                align_category = vertcat(trial_events.values.TrialStimX);
            elseif isfield(trial_events.values,'StimFrequence')
                align_category = vertcat(trial_events.values.StimFrequence);
            end
            align_times = cellfun(@(x) ...
                stimOn_times(align_category(1:length(stimOn_times)) == x & ...
                quiescent_trials),num2cell(unique(align_category)),'uni',false);

            % Get PSTH by 2D histogram
            n_units = size(templates,1);
            unit_psth = nan(n_units,length(t_bins)-1,length(align_times));
            for curr_align = 1:length(align_times)
                t_peri_event = align_times{curr_align} + t_bins;

                use_spikes = spike_times_timelite >= min(t_peri_event,[],'all') & ...
                    spike_times_timelite <= max(t_peri_event,[],'all');

                spikes_binned_continuous = histcounts2(spike_times_timelite(use_spikes),spike_templates(use_spikes), ...
                    reshape(t_peri_event',[],1),1:size(templates,1)+1)./psth_bin_size;

                use_continuous_bins = reshape(padarray(true(size(t_peri_event(:,1:end-1)')),[1,0],false,'post'),[],1);
                spikes_binned = permute(reshape(spikes_binned_continuous(use_continuous_bins,:), ...
                    size(t_peri_event,2)-1,size(t_peri_event,1),size(templates,1)),[3,1,2]);

                unit_psth(:,:,curr_align) = nanmean(spikes_binned,3);
            end

            % Get responsive units
            switch curr_modality
                case 1
                    response_align = align_times{3}; % response = R vis
                case 2
                    response_align = align_times{2}; % response = M aud
            end

            baseline_t = [-0.2,0];
            response_t = [0,0.2];

            baseline_bins = response_align + baseline_t;
            response_bins = response_align + response_t;

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

            responsive_units = event_response_p > 0.95;

            % Get striatum templates
            striatum_idx = strcmp(probe_areas{1}.name,'Caudoputamen');
            striatum_depth = prctile(probe_areas{1}.probe_depth(striatum_idx,:),[0,100],'all');
            striatum_templates = template_depths >= striatum_depth(1) & ...
                template_depths <= striatum_depth(2);

            % Store unit PSTHs and responsive units (striatum units only)
            unit_psth_all{curr_animal,curr_ap}{curr_day,curr_modality} = unit_psth(striatum_templates,:,:);
            responsive_units_all{curr_animal,curr_ap}{curr_day,curr_modality} = responsive_units(striatum_templates);

            end
        end
    end
    ap.print_progress_fraction(curr_animal,length(animals));
end

data_dir = 'C:\Users\petersa\Documents\PetersLab\analysis\ant_pos_striatum\data';
data_fn = fullfile(data_dir,'unit_data');
save(data_fn,'unit_psth_all','responsive_units_all');


%% ^^ Unit data analysis
% (some dirty plots)

animals = {'AP022','DS007','DS010','DS011', ... % V-A, performed both
    'DS000','DS013','DS014','DS015','DS016', ... % A-V, performed both
    'DS005', ... % A-V, frequency task
    'AP021','DS001', ... % V-A, didn't perform A mixed
    'DS003','DS004', ... % A-V, didn't perform A mixed
    };

% Load unit and performance data
data_dir = 'C:\Users\petersa\Documents\PetersLab\analysis\ant_pos_striatum\data';
load(fullfile(data_dir,'unit_data'));
load(fullfile(data_dir,'bhv_p'));

% (get days with auditory performance?)
cell2mat(horzcat(bhv_p{:,1})') < 0.05;
cell2mat(horzcat(bhv_p{:,2})') < 0.05;

% (manually define V>A)
training_order = ~ismember(animals', ...
    {'AP022','DS007','DS010','DS011','AP021','DS001'});

% Get total fraction of responsive units by modality
responsive_units_frac = cellfun(@(x) arrayfun(@(curr_day) ...
    [mean(x{curr_day,1}), ... % V
    mean(x{curr_day,2})], ... % A
    (1:size(x,1))','uni',false), ...
    responsive_units_all,'uni',false);

r = cellfun(@(x) nanmean(cell2mat(x),1),responsive_units_frac,'uni',false);
figure;
for str_ap = 1:2
    r2 = cell2mat(r(:,str_ap));
    nexttile;
    bar(categorical({'V>A','A>V'}),ap.groupfun(@nanmean,r2,training_order,[]))
    legend({'V','A'});
    ylabel('Fraction responsive units');
end
linkaxes(gcf().Children.Children);


% Get overlap in responsive units by modality
responsive_units_n = cellfun(@(x) arrayfun(@(curr_day) ...
    [sum(x{curr_day,1} & x{curr_day,2}), ... % V+A
    sum(x{curr_day,1} & ~x{curr_day,2}), ... % V
    sum(~x{curr_day,1} & x{curr_day,2})], ...% A
    (1:size(x,1))','uni',false), ...
    responsive_units_all,'uni',false);

responsive_units_overlap_frac = ...
    cellfun(@(x) cellfun(@(x) x./sum(x),x,'uni',false), ...
    responsive_units_n,'uni',false);

a = cell2mat(vertcat(responsive_units_overlap_frac{:,1}));

r = cellfun(@(x) nanmean(cell2mat(x),1),responsive_units_overlap_frac,'uni',false);
figure;
for str_ap = 1:2
    r2 = cell2mat(r(:,str_ap));
    nexttile;
    bar(categorical({'V>A','A>V'}),ap.groupfun(@nanmean,r2,training_order,[]))
    legend({'V+A','V','A'});
    ylabel('Fraction responsive units');
end
linkaxes(gcf().Children.Children);


% Smooth and normalize PSTH data
unit_psth_all_smooth = cellfun(@(x) cellfun(@(x) ...
    smoothdata(x,2,'gaussian',100),x,'uni',false), ...
    unit_psth_all,'uni',false);

softnorm = 1;
unit_psth_all_norm = cellfun(@(x) cellfun(@(x) ...
    (x-nanmean(x(:,1:500,:,:),[2,3]))./(nanmean(x(:,1:500,:,:),[2,3])+softnorm),x,'uni',false), ...
    unit_psth_all_smooth,'uni',false);


% Plot PSTHs for V/A units during V/A stim
str_ap = 2;
use_animals = training_order == 0;

vis_psth_cat = cell2mat(cellfun(@(x) vertcat(x{:,1}),unit_psth_all_norm(use_animals,str_ap),'uni',false));
aud_psth_cat = cell2mat(cellfun(@(x) vertcat(x{:,2}),unit_psth_all_norm(use_animals,str_ap),'uni',false));

v_units_cat = find(cell2mat(cellfun(@(x) vertcat(x{:,1}),responsive_units_all(use_animals,str_ap),'uni',false)));
[~,v_sort_idx] = sort(nanmean(vis_psth_cat(v_units_cat,500:700,3),2));

a_units_cat = find(cell2mat(cellfun(@(x) vertcat(x{:,2}),responsive_units_all(use_animals,str_ap),'uni',false)));
[~,a_sort_idx] = sort(nanmean(aud_psth_cat(a_units_cat,500:700,2),2));

figure; colormap(AP_colormap('BWR'));
h = tiledlayout(2,2);

nexttile;
imagesc(vis_psth_cat(v_units_cat(v_sort_idx),:,3));
clim([-10,10]);
title('Vis (V units)')

nexttile;
imagesc(aud_psth_cat(v_units_cat(v_sort_idx),:,2));
clim([-10,10]);
title('Aud (V units)')

nexttile;
imagesc(vis_psth_cat(a_units_cat(a_sort_idx),:,3));
clim([-10,10]);
title('Vis (A units)')

nexttile;
imagesc(aud_psth_cat(a_units_cat(a_sort_idx),:,2));
clim([-10,10]);
title('Aud (A units)')


% Plot PSTHs for V/A units during V/A stim
str_ap = 2;
use_animals = training_order==1;

vis_psth_cat = cell2mat(cellfun(@(x) vertcat(x{:,1}),unit_psth_all_norm(use_animals,str_ap),'uni',false));
aud_psth_cat = cell2mat(cellfun(@(x) vertcat(x{:,2}),unit_psth_all_norm(use_animals,str_ap),'uni',false));

responsive_units_cat = find(cell2mat(cellfun(@(x) vertcat(x{:,1}),responsive_units_all(use_animals,str_ap),'uni',false)));
[~,sort_idx] = sort(nanmean(vis_psth_cat(responsive_units_cat,500:700,3),2));

figure; colormap(AP_colormap('BWR'));
h = tiledlayout(1,2);

nexttile;
imagesc(vis_psth_cat(responsive_units_cat(sort_idx),:,3));
clim([-10,10]);
title('Vis (V units)')

nexttile;
imagesc(aud_psth_cat(responsive_units_cat(sort_idx),:,2));
clim([-10,10]);
title('Aud (V units)')



% NOTE: organization is {A,B} -> {day,V/A}
str_ap = 2;
responsive_n_cat = cell2mat(vertcat(responsive_units_n{:,str_ap}));
responsive_frac_cat = cell2mat(vertcat(responsive_units_frac{:,str_ap}));

[~,sort_idx] = sort(responsive_n_cat(:,2),'descend');

figure;tiledlayout(2,1);
nexttile;
plot(responsive_n_cat(sort_idx,:))
ylabel('N units')
legend({'V+A','V','A'})

nexttile;
plot(responsive_frac_cat(sort_idx,:));
ylabel('Frac units')
legend({'V','A'})


rec_day = cell2mat(cellfun(@(x) (1:length(x))',responsive_units_n(:,str_ap),'uni',false));
rec_animal = cellfun(@(x,a) repmat({a},length(x),1),responsive_units_n(:,str_ap),animals','uni',false);
rec_animal = vertcat(rec_animal{:});

x = [rec_animal(sort_idx),num2cell(rec_day(sort_idx))];

























