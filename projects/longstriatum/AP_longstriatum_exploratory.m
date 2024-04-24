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

    AP_print_progress_fraction(curr_animal_idx,length(animals));

end

%% Batch widefield (master U)

% Loop through animals
animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021'};

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

    AP_print_progress_fraction(curr_animal,length(animals));

end

% Load master U
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');

load(master_U_fn);

% Load learned day
load("C:\Users\petersa\Desktop\learned_day_all.mat");
ld = cellfun(@(x,ld) (1:size(x,3))'-ld,day_V_all,num2cell(learned_day_all)','uni',false);

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
    'AM018','AM019','AM021'};

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

    AP_print_progress_fraction(curr_animal,length(animals));

end


load("C:\Users\petersa\Desktop\learned_day_all.mat");

use_t = t_centers > 0.1 & t_centers < 0.25;
use_roi = 1;
a = cellfun(@(x) squeeze(nanmean(x(:,use_t,use_roi),2)),day_roi_all,'uni',false);

ld = cellfun(@(x,ld) (1:length(x))'-ld,a,num2cell(learned_day_all)','uni',false);

[m,s] = grpstats(cell2mat(a),cell2mat(ld),{'mean','sem'});

ld_x = unique(cell2mat(ld));

figure;errorbar(ld_x,m,s,'k','linewidth',2);
xline(0);

xlabel('Learned day');
ylabel('\DeltaF/F_0')



%% Batch MUA by depth

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021'};

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set depth groups
n_depths = 3;
depth_group_edges = round(linspace(0,4000,n_depths+1));

day_mua_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
%     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_mua = nan(n_depths,length(t_centers),length(recordings));

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Make depth groups
        [depth_group_n,~,depth_group] = histcounts(spike_depths,depth_group_edges);
        depth_groups_used = unique(depth_group);
        depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

        n_depths = length(depth_group_edges) - 1;
        depth_group = discretize(spike_depths,depth_group_edges);

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

    AP_print_progress_fraction(curr_animal,length(animals));

end

load("C:\Users\petersa\Desktop\learned_day_all.mat");

use_t = t_centers > 0 & t_centers < 0.2;
a = cellfun(@(x) squeeze(nanmean(x(2,use_t,:),2)),day_mua_all,'uni',false);

ld = cellfun(@(x,ld) (1:length(x))'-ld,a,num2cell(learned_day_all)','uni',false);

[m,s] = grpstats(cell2mat(a),cell2mat(ld),{'mean','sem'});

ld_x = unique(cell2mat(ld));

figure;errorbar(ld_x,m,s,'k','linewidth',2);
xline(0);

xlabel('Learned day');
ylabel('\DeltaFR/FR_0')





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