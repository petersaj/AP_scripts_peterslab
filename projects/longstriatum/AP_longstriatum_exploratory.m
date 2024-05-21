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
load("C:\Users\petersa\Desktop\am_bhv.mat");
ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day},'uni',false);

a = cat(3,day_V_all{:});

r = grpstats(reshape(a,[],size(a,3))',vertcat(ld{:}));
r2 = plab.wf.svd2px(U_master,reshape(r',[size(a,[1,2]),size(r,1)]));

ap.imscroll(r2);axis image


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



%% Batch MUA by depth

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





%% Batch cortical maps by striatal depth

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
n_depths = 4;

day_ctxmap_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
%     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    U_size = [450,426];
    day_ctxmap = nan(U_size(1),U_size(2),n_depths,length(recordings));

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
        kernel_t = [0,0];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        lambda = 5;
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

        k_px = plab.wf.svd2px(wf_U(:,:,use_svs),k);

        day_ctxmap(:,:,:,curr_recording) = squeeze(k_px);

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_ctxmap_all{curr_animal} = day_ctxmap;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Plot average cortical map by learned day and depth
load("C:\Users\petersa\Desktop\am_bhv.mat");
ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day},'uni',false);
ld_x = unique(cell2mat(ld));

day_ctxmap_ld_avg = ap.groupfun(cat(4,day_ctxmap_all{:}),[],[],[],vertcat(ld{:}));

figure;
c = [-1,1].*prctile(day_ctxmap_ld_avg,99,'all');
h = tiledlayout(size(day_ctxmap_ld_avg,3),size(day_ctxmap_ld_avg,4),'TileSpacing','none');
for curr_depth = 1:n_depths
    for curr_ld = 1:length(ld_x)
        nexttile;
        imagesc(day_ctxmap_ld_avg(:,:,curr_depth,curr_ld));
        ap.wf_draw('ccf','k');
        axis image off;
        clim(c);
        colormap(ap.colormap('PWG'));
        if curr_depth == 1
            title(ld_x(curr_ld));
        end
        drawnow;
    end
end



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

%% Widefield/ephys regression maps (MUA)

% Set upsample value for regression
upsample_factor = 1;
sample_rate = (1/mean(diff(wf_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% (group multiunit by evenly spaced depths)
n_depths = 10;
depth_group_edges = round(linspace(0,4000,n_depths+1));
[depth_group_n,~,depth_group] = histcounts(spike_depths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

% Draw units and borders
figure('units','normalized','position',[0.05,0.2,0.1,0.7]);
unit_axes = axes('YDir','reverse'); hold on; xlabel('Norm spike rate');ylabel('Depth');

if exist('probe_areas','var')
    probe_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
        probe_areas{1}.color_hex_triplet,'uni',false)),[1,3,2]);

    probe_areas_boundaries = probe_areas{1}.probe_depth;
    probe_areas_centers = mean(probe_areas_boundaries,2);

    probe_areas_image_depth = 0:1:max(probe_areas_boundaries,[],'all');
    probe_areas_image_idx = interp1(probe_areas_boundaries(:,1), ...
        1:height(probe_areas{1}),probe_areas_image_depth, ...
        'previous','extrap');
    probe_areas_image = probe_areas_rgb(probe_areas_image_idx,:,:);

    image(unit_axes,[0,1],probe_areas_image_depth,probe_areas_image);
    yline(unique(probe_areas_boundaries(:)),'color','k','linewidth',1);
    set(unit_axes,'YTick',probe_areas_centers,'YTickLabels',probe_areas{1}.acronym);
end

norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
unit_dots = scatter3( ...
    norm_spike_n,template_depths(unique(spike_templates)), ...
    unique(spike_templates),20,'k','filled');
multiunit_lines = arrayfun(@(x) line(xlim,[0,0],'linewidth',2,'visible','off'),1:2);
xlim(unit_axes,[-0.1,1]);
ylim([-50, max(channel_positions(:,2))+50]);
ylabel('Depth (\mum)')
xlabel('Normalized log rate')

yline(depth_group_edges,'linewidth',2,'color','r');
depth_group_centers = movmean(depth_group_edges,2,'endpoints','discard');
text(zeros(length(depth_group_centers),1),depth_group_centers, ...
    num2cell(1:length(depth_group_centers)),'FontSize',20','color','r');
drawnow;

% Bin spikes
n_depths = length(depth_group_edges) - 1;
depth_group = discretize(spike_depths,depth_group_edges);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timelite(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

use_svs = 1:100;
kernel_t = [0,0];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 5;
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

% Convert kernel to pixel space
r_px = plab.wf.svd2px(wf_U(:,:,use_svs),k);

AP_imscroll(r_px,kernel_frames/wf_framerate);
clim([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
colormap(AP_colormap('BWR'));
axis image;

% Get center of mass for each pixel
% (get max r for each pixel, filter out big ones)
r_px_max = squeeze(r_px(:,:,kernel_frames==0,:));% squeeze(max(r_px,[],3));
r_px_max(isnan(r_px_max)) = 0;
r_px_max_norm = bsxfun(@rdivide,r_px_max, ...
    permute(max(reshape(r_px_max,[],n_depths),[],1),[1,3,2]));
r_px_max_norm(isnan(r_px_max_norm)) = 0;
r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depths,[1,3,2])),3)./sum(r_px_max_norm,3);

% Plot map of cortical pixel by preferred depth of probe
r_px_com_col = ind2rgb(round(mat2gray(r_px_com,[1,n_depths])*255),jet(255));
figure;
a1 = axes('YDir','reverse');
imagesc(wf_avg); colormap(gray); clim([0,prctile(wf_avg(:),99.7)]);
axis off; axis image;
a2 = axes('Visible','off');
p = imagesc(r_px_com_col);
axis off; axis image;
set(p,'AlphaData',mat2gray(max(r_px_max_norm,[],3), ...
    [0,double(prctile(reshape(max(r_px_max_norm,[],3),[],1),95))]));
set(gcf,'color','w');

c1 = colorbar('peer',a1,'Visible','off');
c2 = colorbar('peer',a2);
ylabel(c2,'Depth (\mum)');
colormap(c2,jet);
set(c2,'YDir','reverse');
set(c2,'YTick',linspace(0,1,n_depths));
set(c2,'YTickLabel',round(linspace(depth_group_edges(1),depth_group_edges(end),n_depths)));