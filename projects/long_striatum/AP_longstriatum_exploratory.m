%% Get learned day

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021'};

% Grab learning day for each mouse
learned_day_all = nan(size(animals));

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};

    use_workflow = {'stim_wheel*'};
    recordings = plab.find_recordings(animal,[],use_workflow);

    surround_time = [-5,5];
    surround_sample_rate = 100;
    surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

    n_trials_water = nan(length(recordings),2);
    frac_move_day = nan(length(recordings),1);
    rxn_med = nan(length(recordings),1);
    frac_move_stimalign = nan(length(recordings),length(surround_time_points));
    rxn_stat_p = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get total trials/water
        n_trials_water(curr_recording,:) = [length(trial_events.timestamps), ...
            sum(([trial_events.values.Outcome] == 1)*6)];

        % Get median stim-outcome time
        n_trials = length([trial_events.timestamps.Outcome]);
        rxn_med(curr_recording) = median(seconds([trial_events.timestamps(1:n_trials).Outcome] - ...
            cellfun(@(x) x(1),{trial_events.timestamps(1:n_trials).StimOn})));

        % Align wheel movement to stim onset
        align_times = stimOn_times;
        pull_times = align_times + surround_time_points;

        frac_move_day(curr_recording) = nanmean(wheel_move);

        event_aligned_wheel_vel = interp1(timelite.timestamps, ...
            wheel_velocity,pull_times);
        event_aligned_wheel_move = interp1(timelite.timestamps, ...
            +wheel_move,pull_times,'previous');

        frac_move_stimalign(curr_recording,:) = nanmean(event_aligned_wheel_move,1);

        % Get association stat
        rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move);

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});
    end

    learned_day = rxn_stat_p < 0.05 & rxn_med < 2;
    learned_day_all(curr_animal_idx) = find(learned_day,1);

    ap.print_progress_fraction(curr_animal_idx,length(animals));

end

%% Batch widefield (master U)

% Loop through animals
animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022'};

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

% Load learned day
am_data_path = 'C:\Users\petersa\Desktop\am_temp';
load(fullfile(am_data_path,'bhv.mat'));

ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day},'uni',false);

a = cat(3,day_V_all{:});

use_animals = ~cellfun(@isempty,day_V_all);
r = grpstats(reshape(a,[],size(a,3))',vertcat(ld{use_animals}));
r2 = plab.wf.svd2px(U_master,reshape(r',[size(a,[1,2]),size(r,1)]));

ld_x = unique(vertcat(ld{use_animals}));

ap.imscroll(r2,ld_x);axis image

% Plot ROI
% (draw ROI first)
figure;
imagesc(t_centers,ld_x,roi.trace)




%% Batch widefield ROI

% Define ROI
[~,roi_mask] = ap.wf_roi([],[],'master');

% Loop through animals
animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022'};

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.01;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

day_roi_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

%     use_workflow = 'lcr_passive';
    use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_roi = nan(length(recordings),length(t_centers));

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.widefield = true;
        ap.load_recording;

        % If widefield isn't aligned, skip
        if ~load_parts.widefield_align
            continue
        end

        % Get roi trace
        roi_trace = ap.wf_roi(wf_U,wf_V,[],[],roi_mask);

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
        
        aligned_trace = reshape(interp1(wf_t,roi_trace',peri_event_t,'previous'), ...
            length(align_times),length(t_centers),[]);

        aligned_trace_baselinesub = aligned_trace - ...
            mean(aligned_trace(:,t_centers < 0),2);

        % Store mean aligned ROI trace
        day_roi(curr_recording,:) = nanmean(aligned_trace_baselinesub,1);

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_roi_all{curr_animal} = day_roi;

    ap.print_progress_fraction(curr_animal,length(animals));

end


load("C:\Users\petersa\Desktop\am_bhv.mat");

use_t = t_centers > 0.1 & t_centers < 0.25;
use_roi = 1;
a = cellfun(@(x) squeeze(nanmean(x(:,use_t,use_roi),2)),day_roi_all,'uni',false);

ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day},'uni',false);

[m,s] = grpstats(cell2mat(a),cell2mat(ld),{'mean','sem'});

ld_x = unique(cell2mat(ld));

figure;errorbar(ld_x,m,s,'k','linewidth',2);
xline(0);

xlabel('Learned day');
ylabel('\DeltaF/F_0')



%% Batch MUA by depth (n depths)

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022'};

plot_depths = true;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set depth groups
n_depths = 3;

day_mua_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
%     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_mua = nan(n_depths,length(t_centers),length(recordings));

    if plot_depths
        figure;
        h = tiledlayout(1,length(recordings));
        title(h,animal);
    end

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Get striatum start = lowest unit density, end = end of probe
        unit_density_bins = 0:100:3840;
        unit_density = histcounts(template_depths,unit_density_bins);
        [~,unit_density_min_bottom_idx] = min(fliplr(unit_density));
        unit_density_min_idx = length(unit_density_bins) - unit_density_min_bottom_idx;
        template_depths_sorted = sort(template_depths);
        str_start =  template_depths_sorted(find(template_depths_sorted >= ...
            unit_density_bins(unit_density_min_idx+1),1));
        str_end = max(channel_positions(:,2));

        % Discretize units by depth groups
        depth_group_edges = round(linspace(str_start,str_end,n_depths+1));
        depth_group = discretize(spike_depths,depth_group_edges);

        % Plot units and depth groupings
        if plot_depths
            nexttile(h);
            set(gca,'YDir','reverse'); hold on
            norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
            unit_dots = scatter(norm_spike_n,template_depths( ...
                unique(spike_templates)),20,'k','filled');
            yline(depth_group_edges,'color','r');
            drawnow;
        end

        % Stim times (quiescent trials only)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        if isfield(trial_events.values,'TrialStimX')
            stim_x = vertcat(trial_events.values.TrialStimX);
            % sometimes not all trials shown?
            stim_x = stim_x(1:length(stimOn_times));
            align_times = {stimOn_times(stim_x == 90 & quiescent_trials)};
        else
            % Stim times (task)
            align_times = {stimOn_times};
        end

        depth_psth = nan(n_depths,length(t_bins)-1,length(align_times));
        for curr_align = 1:length(align_times)
            use_align = align_times{curr_align};
            t_peri_event = bsxfun(@plus,use_align,t_bins);
            for curr_depth = 1:n_depths
                curr_spikes = spike_times_timelite(depth_group == curr_depth);

                curr_spikes_binned = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_peri_event(x,:)), ...
                    [1:size(t_peri_event,1)]','uni',false))./psth_bin_size;
                curr_mean_psth = mean(curr_spikes_binned,1);

                depth_psth(curr_depth,:,curr_align) = curr_mean_psth;
            end
        end

        smooth_size = 100;
        depth_psth_smooth = smoothdata(depth_psth,2,'gaussian',smooth_size);

        depth_psth_smooth_baseline = nanmean(nanmean(depth_psth_smooth(:,t_centers<0,:),2),3);

        softnorm = 150;
        depth_psth_smooth_norm = (depth_psth_smooth - depth_psth_smooth_baseline)./ ...
            (depth_psth_smooth_baseline+softnorm);

        % Store MUA
        day_mua(:,:,curr_recording) = depth_psth_smooth_norm;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_mua_all{curr_animal} = day_mua;

    ap.print_progress_fraction(curr_animal,length(animals));

end

load("C:\Users\petersa\Desktop\am_bhv.mat");

plot_depth = 1; 
use_t = t_centers > 0.05 & t_centers < 0.15;
a = cellfun(@(x) squeeze(mean(x(plot_depth,use_t,:),2)),day_mua_all,'uni',false);
b = cell2mat(cellfun(@(x) permute(x(1,:,:),[3,2,1]),day_mua_all,'uni',false));

ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day}','uni',false);

[m,s,n] = grpstats(cell2mat(a),cell2mat(ld),{'mean','sem','numel'});
r = grpstats(b,cell2mat(ld),'mean');

ld_x = unique(cell2mat(ld));

figure;
subplot(1,2,1);
imagesc(t_centers,unique(vertcat(ld{:})),r);
subplot(1,2,2); hold on;
cellfun(@(x,ld) plot(ld(1:min(length(ld),length(x))),x(1:min(length(ld),length(x)))),a,ld);
errorbar(ld_x,m,s,'k','linewidth',2);
xline(0);
xlabel('Learned day');
ylabel('\DeltaFR/FR_0');

figure;
subplot(1,3,1); hold on;
cellfun(@(rxn,data) plot(...
    rxn(1:min(length(rxn),length(data))), ...
    data(1:min(length(rxn),length(data)))), ...
    {bhv.rxn_med}',a);
xlabel('Reaction time');
ylabel('MUA');
set(gca,'XScale','log');

subplot(1,3,2); hold on;
cellfun(@(rxn,data) plot(...
    rxn(1:min(length(rxn),length(data))), ...
    data(1:min(length(rxn),length(data)))), ...
    {bhv.stim_move_frac_ratio}',a);
xlabel('Move ratio');
ylabel('MUA');

subplot(1,3,3); hold on;
cellfun(@(rxn,ratio,data) plot3(...
    rxn(1:min(length(rxn),length(data))), ...
    ratio(1:min(length(rxn),length(data))), ...
    data(1:min(length(rxn),length(data)))), ...
    {bhv.rxn_med}',{bhv.stim_move_frac_ratio}',a);
xlabel('Reaction time');
ylabel('Move ratio');
zlabel('MUA');
set(gca,'XScale','log');


% Do above with manually excluded days
exclude_days = cellfun(@(x) false(size(x)),ld,'uni',false);
exclude_days{strcmp(animals,'AP008')}([2,3,4]) = true;
exclude_days{strcmp(animals,'AM011')}([1]) = true;
exclude_days{strcmp(animals,'AM012')}([1,4]) = true;
exclude_days{strcmp(animals,'AM016')}([1,2]) = true;
exclude_days{strcmp(animals,'AM017')}([3,4,7,8]) = true;
exclude_days{strcmp(animals,'AM018')}([1,6,7]) = true;
exclude_days{strcmp(animals,'AM019')}([1,2,3]) = true;
exclude_days{strcmp(animals,'AM022')}([2]) = true;

day_mua_subset = day_mua_all;
for curr_animal = 1:length(day_mua_subset)
    day_mua_subset{curr_animal}(:,:,exclude_days{curr_animal}) = NaN;
end

plot_depth = 1; 
use_t = t_centers > 0.05 & t_centers < 0.15;
a = cellfun(@(x) squeeze(mean(x(plot_depth,use_t,:),2)),day_mua_subset,'uni',false);
b = cell2mat(cellfun(@(x) permute(x(1,:,:),[3,2,1]),day_mua_subset,'uni',false));
ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day},'uni',false);

[m,s,n] = grpstats(cell2mat(a),cell2mat(ld),{'mean','sem','numel'});
r = grpstats(b,cell2mat(ld),'mean');

ld_x = unique(cell2mat(ld));

figure;
subplot(1,2,1);
imagesc(t_centers,unique(vertcat(ld{:})),r);
subplot(1,2,2); hold on;
cellfun(@(x,ld) plot(ld,x),a,ld);
errorbar(ld_x,m,s,'k','linewidth',2);
xline(0);
xlabel('Learned day');
ylabel('\DeltaFR/FR_0');

% Plot above only with subset of animals
on_target = [1,2,5,6,7,12];
off_target = [3,4,8,9,10,11,13];
use_animals = on_target;

a = cellfun(@(x) squeeze(mean(x(plot_depth,use_t,:),2)),day_mua_all,'uni',false);
b = cell2mat(cellfun(@(x) permute(x(1,:,:),[3,2,1]),day_mua_all,'uni',false));
ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day},'uni',false);

a2 = a(use_animals);
ld2 = ld(use_animals);
[m,s,n] = grpstats(cell2mat(a2),cell2mat(ld2),{'mean','sem','numel'});

figure;
subplot(1,2,1);
imagesc(t_centers,unique(vertcat(ld2{:})),r);
subplot(1,2,2); hold on;
cellfun(@(x,ld) plot(ld,x),a2,ld2);

ld_x = unique(cell2mat(ld2));
errorbar(ld_x,m,s,'k','linewidth',2);
xline(0);
xlabel('Learned day');
ylabel('\DeltaFR/FR_0');


% Plot responses of all animal separately
figure;
h = tiledlayout(2,length(animals),'TileIndexing','columnmajor');
for curr_animal = 1:length(animals)
    a = permute(day_mua_all{curr_animal}(1,:,:),[3,2,1]);
    h1 = nexttile;
    imagesc(t_centers,ld{curr_animal},a);
    title(h1,animals{curr_animal});
    nexttile;hold on
    set(gca,'ColorOrder',copper(size(a,1)));
    plot(t_centers,a');
end

%% Batch MUA by depth (size depths)

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022'};

plot_depths = false;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set MUA length (microns)
mua_length = 200;

day_mua_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

%     use_workflow = 'lcr_passive';
    use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_mua = cell(length(recordings),1);

    if plot_depths
        figure;
        h = tiledlayout(1,length(recordings));
        title(h,animal);
    end

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Get striatum start = lowest unit density, end = end of probe
        unit_density_bins = 0:100:3840;
        unit_density = histcounts(template_depths,unit_density_bins);
        [~,unit_density_min_bottom_idx] = min(fliplr(unit_density));
        unit_density_min_idx = length(unit_density_bins) - unit_density_min_bottom_idx;
        template_depths_sorted = sort(template_depths);
        str_start =  template_depths_sorted(find(template_depths_sorted >= ...
            unit_density_bins(unit_density_min_idx+1),1));
        str_end = max(channel_positions(:,2));

        % Discretize spikes by depth
        depth_group_edges = str_start:mua_length:str_end;
        if length(depth_group_edges) < 2
            continue    
        end
        depth_group = discretize(spike_depths,depth_group_edges);

        % Plot units and depth groupings
        if plot_depths
            nexttile(h);
            set(gca,'YDir','reverse'); hold on
            norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
            unit_dots = scatter(norm_spike_n,template_depths( ...
                unique(spike_templates)),20,'k','filled');
            yline(depth_group_edges,'color','r');
            drawnow;
        end

        % Stim times (quiescent trials only)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        if isfield(trial_events.values,'TrialStimX')
            stim_x = vertcat(trial_events.values.TrialStimX);
            % sometimes not all trials shown?
            stim_x = stim_x(1:length(stimOn_times));
            align_times = {stimOn_times(stim_x == 90 & quiescent_trials)};
        else
            % Stim times (task)
            align_times = {stimOn_times};
        end

        depth_psth = nan(max(depth_group),length(t_bins)-1,length(align_times));
        for curr_align = 1:length(align_times)
            use_align = align_times{curr_align};
            t_peri_event = bsxfun(@plus,use_align,t_bins);
            for curr_depth = 1:max(depth_group)
                curr_spikes = spike_times_timelite(depth_group == curr_depth);

                curr_spikes_binned = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_peri_event(x,:)), ...
                    [1:size(t_peri_event,1)]','uni',false))./psth_bin_size;
                curr_mean_psth = mean(curr_spikes_binned,1);

                depth_psth(curr_depth,:,curr_align) = curr_mean_psth;
            end
        end

        smooth_size = 100;
        depth_psth_smooth = smoothdata(depth_psth,2,'gaussian',smooth_size);

        depth_psth_smooth_baseline = nanmean(nanmean(depth_psth_smooth(:,t_centers<0,:),2),3);

        softnorm = 150;
        depth_psth_smooth_norm = (depth_psth_smooth - depth_psth_smooth_baseline)./ ...
            (depth_psth_smooth_baseline+softnorm);

        % Store MUA
        day_mua{curr_recording} = depth_psth_smooth_norm;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_mua_all{curr_animal} = day_mua;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% (same division as cortical maps below)



%% Batch cortical maps by striatal depth

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022'};

% Set MUA length (microns)
mua_length = 200;

ctx_map_all = cell(length(animals),1);
ctx_expl_var_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

%     use_workflow = 'lcr_passive';
    use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    ctx_map = cell(length(recordings),1);
    ctx_expl_var = cell(length(recordings),1);
    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts.ephys = true;
        load_parts.widefield = true;
        load_parts.widefield_align = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % If not aligned, skip
        if ~load_parts.widefield_align 
            continue
        end

        % Whole thing in a try/catch
        try

            % Get striatum start = lowest unit density, end = end of probe
            unit_density_bins = 0:100:3840;
            unit_density = histcounts(template_depths,unit_density_bins);
            [~,unit_density_min_bottom_idx] = min(fliplr(unit_density));
            unit_density_min_idx = length(unit_density_bins) - unit_density_min_bottom_idx;
            template_depths_sorted = sort(template_depths);
            str_start =  template_depths_sorted(find(template_depths_sorted >= ...
                unit_density_bins(unit_density_min_idx+1),1));
            str_end = max(channel_positions(:,2));

            % Discretize spikes by depth
            depth_group_edges = str_start:mua_length:str_end;
            depth_group = discretize(spike_depths,depth_group_edges);

            % Get MUA binned by widefield
            sample_rate = (1/mean(diff(wf_t)));
            skip_seconds = 60;
            time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
            time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

            binned_spikes = zeros(max(depth_group),length(time_bins)-1);
            for curr_depth = 1:max(depth_group)
                curr_spike_times = spike_times_timelite(depth_group == curr_depth);
                binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
            end

            binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
            binned_spikes_std(isnan(binned_spikes_std)) = 0;

            use_svs = 1:500;
            kernel_t = [0,0];
            kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
            lambda = 10;
            zs = [false,false];
            cvfold = 5;
            return_constant = false;
            use_constant = true;

            fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

            % Regress cortex to spikes
            [k,predicted_spikes,explained_var] = ...
                AP_regresskernel(fVdf_deconv_resample, ...
                binned_spikes_std,kernel_frames, ...
                lambda,zs,cvfold,return_constant,use_constant);

            ctx_map{curr_recording} = permute(k,[1,3,2]);
            ctx_expl_var{curr_recording} = explained_var.total;

        catch me
            % Any errors: do nothing
        end

        % Prep for next loop
        clearvars('-except',preload_vars{:});

        ap.print_progress_fraction(curr_recording,length(recordings));

    end

    ctx_map_all{curr_animal} = ctx_map;
    ctx_expl_var_all{curr_animal} = ctx_expl_var;

    fprintf('Done %s\n',animal);

end

% Load master U
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);

% Convert V to px
r = cellfun(@(x) cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),x,'uni',false),ctx_map_all,'uni',false);

animal_day_idx_cell = cellfun(@(animal,data) cellfun(@(day,data) ...
    [repmat(animal,size(data,3),1),repmat(day,size(data,3),1),(1:size(data,3))'], ...
    num2cell(1:length(data))',data,'uni',false), ...
    num2cell(1:length(r))',r,'uni',false);
animal_day_idx = cell2mat(vertcat(animal_day_idx{:}));




%%% Try? 
% ROI in mVIS and mPFC to classify?


% Below won't work - changed above to be by 200um rather than depth

% % Plot average cortical map by learned day and depth
% load("C:\Users\petersa\Desktop\am_bhv.mat");
% ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day},'uni',false);
% ld_x = unique(cell2mat(ld));
% 
% day_ctxmap_ld_avg = ap.groupfun(cat(4,ctx_map_all{:}),[],[],[],vertcat(ld{:}));
% 
% figure;
% c = [-1,1].*prctile(day_ctxmap_ld_avg,99,'all');
% h = tiledlayout(size(day_ctxmap_ld_avg,3),size(day_ctxmap_ld_avg,4),'TileSpacing','none');
% for curr_depth = 1:n_depths
%     for curr_ld = 1:length(ld_x)
%         nexttile;
%         imagesc(day_ctxmap_ld_avg(:,:,curr_depth,curr_ld));
%         ap.wf_draw('ccf','k');
%         axis image off;
%         clim(c);
%         colormap(ap.colormap('PWG'));
%         if curr_depth == 1
%             title(ld_x(curr_ld));
%         end
%         drawnow;
%     end
% end



%% Batch cortical regression by striatal depth

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022'};

plot_depths = false;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set depth groups
n_depths = 3;

day_spikes_all = cell(length(animals),1);
day_predicted_spikes_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
%     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_spikes = nan(n_depths,length(t_centers),length(recordings));
    day_predicted_spikes = nan(n_depths,length(t_centers),length(recordings));
  
    if plot_depths
        figure;
        h = tiledlayout(1,length(recordings));
        title(h,animal);
    end

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts.ephys = true;
        load_parts.widefield = true;
        load_parts.widefield_align = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % If not aligned, skip
        if ~load_parts.widefield_align 
            continue
        end

        % Get striatum start = lowest unit density, end = end of probe
        unit_density_bins = 0:100:3840;
        unit_density = histcounts(template_depths,unit_density_bins);
        [~,unit_density_min_bottom_idx] = min(fliplr(unit_density));
        unit_density_min_idx = length(unit_density_bins) - unit_density_min_bottom_idx;
        template_depths_sorted = sort(template_depths);
        str_start =  template_depths_sorted(find(template_depths_sorted >= ...
            unit_density_bins(unit_density_min_idx+1),1));
        str_end = max(channel_positions(:,2));

        % Discretize spikes by depth
        depth_group_edges = round(linspace(str_start,str_end,n_depths+1));
        depth_group = discretize(spike_depths,depth_group_edges);

        % Get MUA binned by widefield
        sample_rate = (1/mean(diff(wf_t)));
        skip_seconds = 60;
        time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

        binned_spikes = zeros(n_depths,length(time_bins)-1);
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timelite(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end

        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;

        use_svs = 1:100;
        kernel_t = [-0.1,0.1];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        lambda = 5;
        zs = [false,false];
        cvfold = 5;
        return_constant = true;
        use_constant = true;

        fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

        % Regress cortex to spikes
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample, ...
            binned_spikes_std,kernel_frames, ...
            lambda,zs,cvfold,return_constant,use_constant);

        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});

        % Stim times (quiescent trials only)
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

        % Align measured and predicted activity to align times    
        peri_event_t = align_times + t_centers;
        
        aligned_spikes = reshape(interp1(time_bin_centers,binned_spikes', ...
            peri_event_t), ...
            length(align_times),length(t_centers),[]);

        aligned_predicted_spikes = reshape(interp1(time_bin_centers, ...
            predicted_spikes',peri_event_t), ...
            length(align_times),length(t_centers),[]);


        % Store mean aligned ROI trace
        day_spikes(:,:,curr_recording) = permute(nanmean(aligned_spikes,1),[3,2,1]);
        day_predicted_spikes(:,:,curr_recording) = permute(nanmean(aligned_predicted_spikes,1),[3,2,1]);

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_spikes_all{curr_animal} = day_spikes;
    day_predicted_spikes_all{curr_animal} = day_predicted_spikes;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Plot average measured and predicted responses
softnorm = 0;
day_spikes_all_norm = cellfun(@(x) (x-nanmean(x(:,t_centers<0,:),2))./(nanmean(x(:,t_centers<0,:),2)+softnorm),day_spikes_all,'uni',false);
day_predicted_spikes_all_norm = cellfun(@(x) (x-nanmean(x(:,t_centers<0,:),2))./(nanmean(x(:,t_centers<0,:),2)+softnorm),day_predicted_spikes_all,'uni',false);

load("C:\Users\petersa\Desktop\am_bhv.mat");

use_depth = 1; 
use_t = t_centers > 0.05 & t_centers < 0.2;
a = cellfun(@(x) squeeze(mean(x(use_depth,use_t,:),2)),day_predicted_spikes_all_norm,'uni',false);
b = cell2mat(cellfun(@(x) permute(x(1,:,:),[3,2,1]),day_predicted_spikes_all_norm,'uni',false));

ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day},'uni',false);

[m,s,n] = grpstats(cell2mat(a),cell2mat(ld),{'mean','sem','numel'});
r = grpstats(b,cell2mat(ld),'mean');

ld_x = unique(cell2mat(ld));

figure;
subplot(1,2,1);
imagesc(t_centers,unique(vertcat(ld{:})),r);
subplot(1,2,2);errorbar(ld_x,m,s,'k','linewidth',2);
xline(0);
xlabel('Learned day');
ylabel('\DeltaFR/FR_0')


use_t = t_centers > 0.05 & t_centers < 0.2;
a = cellfun(@(x) squeeze(mean(x(use_depth,use_t,:),2)),day_spikes_all_norm,'uni',false);
b = cellfun(@(x) squeeze(mean(x(use_depth,use_t,:),2)),day_predicted_spikes_all_norm,'uni',false);

figure; tiledlayout(1,3);
nexttile; hold on;
cellfun(@(x,ld) plot(ld,x),a,ld); title('Measured');
nexttile; hold on;
cellfun(@(x,ld) plot(ld,x),b,ld); title('Predicted');

nexttile
hold on; axis equal
cellfun(@(x,y) plot(x,y),a,b);
cellfun(@(x,y,ld) plot(x(ld==0),y(ld==0),'.k','MarkerSize',20),a,b,ld);
axis tight;
line(xlim,xlim,'color','k');
xlabel('Measured');
ylabel('Predicted');


%% Load MUA / cortex maps and process

am_data_path = 'C:\Users\petersa\Desktop\am_temp';
load(fullfile(am_data_path,'bhv.mat'));
load(fullfile(am_data_path,'mua_passive.mat'));
% load(fullfile(am_data_path,'mua_task.mat'));
load(fullfile(am_data_path,'ctx_maps_task.mat'));
load(fullfile(am_data_path,'wf.mat'));


% Load master U, convert V to px
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);
ctx_map_px = cellfun(@(x) cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),x,'uni',false),ctx_map_all,'uni',false);

use_animals = ~cellfun(@isempty,ctx_map_px);

ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day}','uni',false);

% (learned day)
animal_ld_idx_cell = cellfun(@(animal,ld,data) cellfun(@(ld,data) ...
    [repmat(animal,size(data,3),1),repmat(ld,size(data,3),1),(1:size(data,3))'], ...
    num2cell(ld),data,'uni',false), ...
    num2cell(find(use_animals)),ld(use_animals),ctx_map_px(use_animals),'uni',false);
animal_ld_idx = cell2mat(vertcat(animal_ld_idx_cell{:}));

% (cardinal day)
animal_day_idx_cell = cellfun(@(animal,ld,data) cellfun(@(ld,data) ...
    [repmat(animal,size(data,3),1),repmat(ld,size(data,3),1),(1:size(data,3))'], ...
    num2cell(1:length(ld))',data,'uni',false), ...
    num2cell(find(use_animals)),ld(use_animals),ctx_map_px(use_animals),'uni',false);
animal_day_idx = cell2mat(vertcat(animal_day_idx_cell{:}));

% Concatenate V's
V_cat = cat(3,day_V_all{:});

V_animal_day_idx = cell2mat(cellfun(@(x,animal,ld) ...
    [repmat(animal,size(x,3),1),ld], ...
    day_V_all(use_animals),num2cell(find(use_animals)),ld(use_animals),'uni',false));

% Concatenate maps
c2 = cell2mat(permute(vertcat(ctx_map_px{:}),[2,3,1]));

ap.imscroll(c2); 
axis image off;
clim(max(abs(clim)).*[-1,1].*0.5);
colormap(ap.colormap('BWR',[],1.5));
ap.wf_draw('ccf','k');



% (draw ROI over mVis)
vis_map_weights = roi.trace';

% (draw ROI over mpfc)
mpfc_map_weights = roi.trace';


r = cell2mat(vertcat(day_mua_all{use_animals}));

use_t = [500,700];
r2 = nanmean(r(:,use_t(1):use_t(2)),2);

% Plot vis/mpfc weights and MUA
figure; colormap(jet);
subplot(1,2,1);
scatter3(vis_map_weights,mpfc_map_weights,r2,15,animal_day_idx(:,1),'filled')
xlabel('vis');ylabel('mpfc');zlabel('mua');
axis vis3d; shading interp;

subplot(1,2,2);
f = fit([vis_map_weights,mpfc_map_weights],r2,'lowess');
plot(f,[vis_map_weights,mpfc_map_weights],r2)
xlabel('vis');ylabel('mpfc');zlabel('mua');
axis vis3d; shading interp;

% (above, but by learned day)
figure;tiledlayout('flow','tilespacing','compact'); colormap(jet)
for curr_ld = -3:3
    use_data = animal_ld_idx(:,2) == curr_ld;

    nexttile
    f = fit([vis_map_weights(use_data),mpfc_map_weights(use_data)],r2(use_data),'lowess');
    plot(f,[vis_map_weights(use_data),mpfc_map_weights(use_data)],r2(use_data))
    xlabel('vis');ylabel('mpfc');zlabel('mua');
    axis vis3d; shading interp;
    title(curr_ld);
end

% (above, by learned day on one plot)
figure;hold on
for curr_ld = -3:3
    use_data = animal_ld_idx(:,2) == curr_ld;
    f = fit([vis_map_weights(use_data),mpfc_map_weights(use_data)],r2(use_data),'poly22');
    plot(f);
end
xlabel('vis');ylabel('mpfc');zlabel('mua');
axis vis3d;
shading interp;

% correlation between pixels and MUA?
c3 = 1-reshape(pdist2(reshape(c2,[],size(c2,3)),r2','correlation'),size(c2,[1,2]));
figure;imagesc(c3);
axis image off;
clim([-0.5,0.5]);
colormap(ap.colormap('BWR',[],1.5));
ap.wf_draw('ccf','k');

% Plot maps by threshold cutoffs
ctx_thresh = 0.01;
mua_thresh = 0.1;

curr_use_ctx = vis_map_weights >= ctx_thresh;
curr_use_mua = r2 >= mua_thresh;

figure; 
h = tiledlayout(2,2);
colormap(ap.colormap('BWR',[],1.5));

nexttile([2,1]);
plot(vis_map_weights,r2,'.k');
xlabel('Cortex weight');
ylabel('Striatum vis response');
xline(ctx_thresh,'r');
yline(mua_thresh,'r');

nexttile;
imagesc(nanmean(c2(:,:,curr_use_ctx & curr_use_mua),3));
axis image off; clim([-0.03,0.03]);
ap.wf_draw('ccf','k');

nexttile;
imagesc(nanmean(c2(:,:,curr_use_ctx & ~curr_use_mua),3));
axis image off; clim([-0.03,0.03]);
ap.wf_draw('ccf','k');



% combine used depths
% (needs cardinal day index above)
[grp,~,grp_idx] = unique(animal_day_idx(curr_use_ctx,1:2),'rows');

x1 = ap.groupfun(r(curr_use_ctx,:),grp_idx,[]);
x2 = ap.groupfun(c2(:,:,curr_use_ctx),[],[],grp_idx);

ap.imscroll(x2); 
axis image off;
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('BWR',[],1.5));
ap.wf_draw('ccf','k');

x11 = nanmean(x1(:,use_t(1):use_t(2)),2);

stim_move = nan(size(grp,1),1);
for i = 1:size(stim_move,1)
    stim_move(i) = bhv(grp(i,1)).stim_move_frac_ratio(grp(i,2));
end

figure;plot(stim_move,x11,'.k')
xlabel('stim move ratio');
ylabel('MUA');


% Get vis cortex maps, group striatal responses, plot by LD
ctx_thresh = 0.01;
curr_use_ctx = vis_map_weights >= ctx_thresh;

[grp,~,grp_idx] = unique(animal_ld_idx(curr_use_ctx,1:2),'rows');

x1 = ap.groupfun(r(curr_use_ctx,:),grp_idx,[]);
x11t = ap.groupfun(x1,grp(:,2),[]);
x11 = mean(x1(:,use_t(1):use_t(2)),2);

figure;imagesc([],unique(grp(:,2)),x11t)
ylabel('LD');
title('MUA');

x2 = ap.groupfun(c2(:,:,curr_use_ctx),[],[],grp_idx);
x22 = ap.groupfun(x2,[],[],grp(:,2));
figure; tiledlayout('flow','TileSpacing','none');
c = max(abs(x22),[],'all').*[-1,1].*0.8;
ld_x = unique(grp(:,2));
for i = -3:3
    nexttile; 
    imagesc(x22(:,:,ld_x==i));
    axis image off;
    clim(c);
    ap.wf_draw('ccf','k');
    title(i);
end
colormap(ap.colormap('BWR',[],1.5));

ap.imscroll(x22,ld_x);
axis image off;
clim(c);
ap.wf_draw('ccf','k');
colormap(ap.colormap('BWR',[],1.5));


figure; 
subplot(1,2,1); hold on;
arrayfun(@(x) plot(grp(grp(:,1)==x,2), x11(grp(:,1)==x)),unique(grp(:,1)));
xlabel('Learned day');
ylabel('MUA')

[m,e] = grpstats(x11,grp(:,2),{'mean','sem'});
errorbar(unique(grp(:,2)),m,e,'k');

stim_move = nan(size(grp,1),1);
for i = 1:size(stim_move,1)
    stim_move(i) = bhv(grp(i,1)).stim_move_frac_ratio(ld{grp(i,1)} == grp(i,2));
end

subplot(1,2,2); hold on;
arrayfun(@(x) plot(stim_move(grp(:,1)==x), x11(grp(:,1)==x)),unique(grp(:,1)));
xlabel('Stim move ratio');
ylabel('MUA')

% (combine V's from animals/days used above)
use_v = ismember(V_animal_day_idx,grp,'rows');
V_avg = ap.groupfun(V_cat(:,:,use_v),[],[],V_animal_day_idx(use_v,2));
px_avg = plab.wf.svd2px(U_master,V_avg);
ap.imscroll(px_avg); 
axis image;
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));

% K-means cortex maps
n_k = 4;
[kidx,ck] = kmeans(reshape(c2,[],size(c2,3))',n_k,'Distance','correlation');
ck = reshape(ck',[size(c2,[1,2]),n_k]);

figure;tiledlayout('flow','tilespacing','none');
for i = 1:size(ck,3)
    nexttile;
    imagesc(ck(:,:,i));
    axis image off;
    clim(max(abs(clim)).*[-1,1]);
    colormap(ap.colormap('BWR',[],1.5));
    ap.wf_draw('ccf','k');
end


% Group data as above, but using k-means idx as grouper
curr_use_ctx = kidx == 2;

[grp,~,grp_idx] = unique(animal_ld_idx(curr_use_ctx,1:2),'rows');

x1 = ap.groupfun(r(curr_use_ctx,:),grp_idx,[]);
x11 = mean(x1(:,use_t(1):use_t(2)),2);

x2 = ap.groupfun(c2(:,:,curr_use_ctx),[],[],grp_idx);
x22 = ap.groupfun(x2,[],[],grp(:,2));
figure; tiledlayout('flow','TileSpacing','none');
c = max(abs(x22),[],'all').*[-1,1].*0.8;
ld_x = unique(grp(:,2));
for i = -3:3
    nexttile; 
    imagesc(x22(:,:,ld_x==i));
    axis image off;
    clim(c);
    ap.wf_draw('ccf','k');
    title(i);
end
colormap(ap.colormap('BWR',[],1.5));

figure; 
subplot(1,2,1); hold on;
arrayfun(@(x) plot(grp(grp(:,1)==x,2), x11(grp(:,1)==x)),unique(grp(:,1)));
xlabel('Learned day');
ylabel('MUA')

[m,e] = grpstats(x11,grp(:,2),{'mean','sem'});
errorbar(unique(grp(:,2)),m,e,'k');

stim_move = nan(size(grp,1),1);
for i = 1:size(stim_move,1)
    stim_move(i) = bhv(grp(i,1)).stim_move_frac_ratio(ld{grp(i,1)} == grp(i,2));
end

subplot(1,2,2); hold on;
arrayfun(@(x) plot(stim_move(grp(:,1)==x), x11(grp(:,1)==x)),unique(grp(:,1)));
xlabel('Stim move ratio');
ylabel('MUA')






