%% Fig 1B: task wheel velocity

% Average wheel velocity aligned to learning 
animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

wheel_all = cell(length(animals), 1);

warning off;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};

    % Find all task recordings
    workflow = {'stim_wheel_right*'};
    recordings = plab.find_recordings(animal, [], workflow);
    
    for curr_rec = 1:length(recordings)

        rec_day = recordings(curr_rec).day;
        rec_time = recordings(curr_rec).recording{end};
        verbose = false;
        load_parts.behavior = true;
        ap.load_recording

        % Align wheel velocity to stim onset
        align_times = stimOn_times;

        surround_time = [-1,2];
        surround_sample_rate = 100;
        surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);
        pull_times = align_times + surround_time_points;

        [wheel_velocity,wheel_move] = ...
            AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);

        wheel_all{curr_animal}{curr_rec,1} = interp1(timelite.timestamps, ...
            wheel_velocity,pull_times);    
    end
    ap.print_progress_fraction(curr_animal,length(animals));
end

% Load behavior and get learned day
data_path = fullfile(plab.locations.server_path,'Lab','Papers','Marica_2025','data');

% Set stat and p-value to define learning day
use_stat = 'firstmove_mean';
learn_stat_p = 0.05;

% Load behavior
load(fullfile(data_path,'bhv'));

% Set "learned_days" and "days_from_learning"
bhv.learned_days = cellfun(@(x) x < learn_stat_p,bhv.(['stimwheel_pval_',use_stat]));
for curr_animal = unique(bhv.animal)'
    curr_rows = strcmp(bhv.animal,curr_animal);
    bhv.days_from_learning(curr_rows) = ...
        (1:sum(curr_rows))' - ...
        max([NaN,find(bhv.learned_days(curr_rows),1)]);
end

% Plot wheel velocity aligned to stim onset by learned day
wheel_rec = cell2mat(cellfun(@(x) nanmedian(x,1),vertcat(wheel_all{:}),'uni',false));
[wheel_ld,wheel_ld_grp] = ap.groupfun(@nanmean,wheel_rec,bhv.days_from_learning);
wheel_ld_sem = ap.groupfun(@AP_sem,wheel_rec,bhv.days_from_learning);

plot_days = -3:2;
max_ld = max(abs(plot_days));
ld_colors = ap.colormap('BKR',max_ld*2+1);
plot_ld_colors = ld_colors(ismember(-max_ld:max_ld,plot_days),:);

plot_idx = ismember(wheel_ld_grp,plot_days);
figure; set(gca,'ColorOrder',plot_ld_colors); 
ap.errorfill(surround_time_points,wheel_ld(plot_idx,:)', ...
    wheel_ld_sem(plot_idx,:)');
xline(0);
set(gca,'XTick',0:0.5:1);
set(gca,'YTick',[]);
xlabel('Time from stimulus onset');
ap.prettyfig;



%% Fig 3B: Classical condition lick raster

% Classical conditioned animal example: HA020

animal = 'HA020';

rec = plab.find_recordings(animal,[], ...
    'visual_operant_lick_two_stim_static_big_stim');

rec_day = rec(end).day;
rec_time = rec(end).recording{end};
verbose = true;
load_parts.behavior = true;
ap.load_recording;

% Trial times and params
rewarded_x = trial_events.parameters.RewardedX;

n_trials = sum(cellfun(@(x) length(x) == 2,{trial_events.timestamps.StimOn}));

stim_x = vertcat(trial_events.values(1:n_trials).TrialX);
stimOn_times = photodiode_on_times(1:n_trials);
stimOff_times = photodiode_off_times(1:n_trials);

trial_static_stim_time = vertcat(trial_events.values(1:n_trials).TrialStimStaticTime);
trial_quiescence_time = vertcat(trial_events.values(1:n_trials).TrialQuiescence);

% Align licks to stim on times
psth_window = [-0.5,3];

[lick_psth_plus,~,lick_t] = ap.psth(lick_times, ...
    stimOn_times(stim_x == rewarded_x),...
    'window',psth_window,'bin_size',0.01,'smoothing',20);

lick_psth_minus = ap.psth(lick_times, ...
    stimOn_times(stim_x ~= rewarded_x),...
        'window',psth_window,'bin_size',0.01,'smoothing',20);

figure; hold on;
plot(lick_t,lick_psth_plus,'linewidth',2);
plot(lick_t,lick_psth_minus,'linewidth',2);
xline([0,trial_events.parameters.StimStaticMin, ...
    trial_events.parameters.StimStaticMax],'k');
set(gca,'YTick',[0,8]); ylim([0,8]);
set(gca,'XTick',[0:2]);xlim([-0.2,2]);
ap.prettyfig;


%% Fig 3C: Naive/classical conditioned striatum

conditions = ["naive","trained"];

mua = cell(size(conditions));
for condition = conditions

    switch condition
        case "naive"
            animal = 'AP034';
            rec_day = '2025-12-05';
            rec_time = '1129';
            use_depths = [1500,2500];
        case "trained"
            animal = 'HA020';
            rec = plab.find_recordings(animal,[],'lcr_passive_corner*');
            rec_day = '2026-02-10';
            rec_time = '1549';
            use_depths = [750,1500];
    end

    ap.load_recording;

    % Quiescent trials
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        (1:length(stimOn_times))');

    % Stim-aligned PSTH
    stim_x = vertcat(trial_events.values.TrialStimX);

    use_spikes = isbetween(spike_depths,use_depths(1),use_depths(2));

    [mua_psth,~,mua_psth_t] = ap.psth(spike_times_timelite(use_spikes), ...
        stimOn_times(quiescent_trials & stim_x == 90), ...
        'window',[-0.2,1],'smoothing',50);

    mua{ismember(conditions,condition)} = mua_psth;

end

baseline_t = mua_psth_t < 0;
mua_norm = cellfun(@(x) (x-mean(x(baseline_t),2))./mean(x(baseline_t),2),mua,'uni',false);
figure;plot(mua_psth_t,cell2mat(mua_norm')','linewidth',2);
xlim([-0.2,0.7]);
ap.prettyfig;
ap.scalebar(0.2,1);
axis off;
xline([0,0.5]);

%% Fig 3D: Faster learning with classical conditioning

animals = {'HA008','HA009','HA010','HA011','HA012','HA013','HA014','HA015'};

data_all = cell(length(animals), 1);

warning off;
for animal_idx=1:length(animals)
    
    animal = animals{animal_idx};
    disp(animal);

    % Find all task recordings
    workflow = 'stim_wheel_right_stage\d';
    recordings = plab.find_recordings(animal, [], workflow);
    
    data_animal = table;
    for use_rec=1:length(recordings)

        rec_day = recordings(use_rec).day;
        rec_time = recordings(use_rec).recording{end};
        verbose = false;
        load_parts.behavior = true;
        ap.load_recording
    
        % Align wheel velocity to stim onset
        align_times = stimOn_times;

        surround_time = [-1,2];
        surround_sample_rate = 100;
        surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);
        pull_times = align_times + surround_time_points;

        [wheel_velocity,wheel_move] = ...
            AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);

        wheel_aligned = interp1(timelite.timestamps,wheel_velocity,pull_times);    

        % Get association p-value in a few ways
        % (first just grab mean null reaction times for each trial)
        [~,~,~,stim_to_move_nullmean] = ...
            AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_move, 'mean');

        % (mean to firstmove)
        [stimwheel_pval_firstmove_mean,stimwheel_rxn_firstmove_mean,stimwheel_rxn_null_firstmove_mean] = ...
            AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_move, 'mean');

        % Save data in table
        data_animal.animal(use_rec) = {animal};
        data_animal.rec_day(use_rec) = {rec_day};

        % (firstmove mean stats)
        data_animal.stimwheel_pval_firstmove_mean(use_rec) = {stimwheel_pval_firstmove_mean};
        data_animal.stimwheel_rxn_firstmove_mean(use_rec) = {stimwheel_rxn_firstmove_mean};
        data_animal.stimwheel_rxn_null_firstmove_mean(use_rec) = {stimwheel_rxn_null_firstmove_mean};

        % (reaction times)
        data_animal.stim_to_move(use_rec) = {stim_to_move(1:n_trials)};
        data_animal.stim_to_outcome(use_rec) = {stim_to_outcome(1:n_trials)};
        data_animal.trial_outcome(use_rec) = {logical(trial_outcome(1:n_trials))};

        % (wheel velocity)
        data_animal.wheel_aligned(use_rec) = {wheel_aligned};
  
        % Print progress
        ap.print_progress_fraction(use_rec,length(recordings))
    end

    data_all{animal_idx} = data_animal;

end
warning on;

cc_ld = cellfun(@(x) find(cell2mat(x.stimwheel_pval_firstmove_mean)<0.05,1),data_all);

% Get learend day for operant-only behavior
data_path = fullfile(plab.locations.server_path,'Lab','Papers','Marica_2025','data');

% Set stat and p-value to define learning day
use_stat = 'firstmove_mean';
learn_stat_p = 0.05;

% Load behavior
load(fullfile(data_path,'bhv'));

% Set "learned_days" and "days_from_learning"
bhv.learned_days = cellfun(@(x) x < learn_stat_p,bhv.(['stimwheel_pval_',use_stat]));
for curr_animal = unique(bhv.animal)'
    curr_rows = strcmp(bhv.animal,curr_animal);
    bhv.days_from_learning(curr_rows) = ...
        (1:sum(curr_rows))' - ...
        max([NaN,find(bhv.learned_days(curr_rows),1)]);
end

oc_ld = cell2mat(cellfun(@(x) find(bhv.learned_days(strcmp(bhv.animal,x)),1), ...
    unique(bhv.animal,'stable'),'uni',false));

% Plot learning days OC vs CC
figure; hold on;
swarmchart([ones(size(oc_ld));2*ones(size(cc_ld))],[oc_ld;cc_ld],'filled');
errorbar([nanmean(oc_ld),nanmean(cc_ld)],[AP_sem(oc_ld),AP_sem(cc_ld)], ...
    'capsize',0,'linewidth',2,'color','k');
ap.prettyfig;
set(gca,'XTick',1:2,'XTickLabel',["Naive","Post-conditioning"]);
set(gca,'YTick',[1:2:8]);ylim([1,8]);


%% Fig 3B: example 2-stim performance

animal = 'DS022';

recordings = plab.find_recordings(animal,[],'*error_repeat');

curr_rec = 8;

rec_day = recordings(curr_rec).day;
rec_time = recordings(curr_rec).recording{end};
load_parts.bhv = true;
ap.load_recording;

trial_stim = vertcat(trial_events.values(1:n_trials).TaskType);

% Align wheel velocity to stim onset
align_times = stimOn_times(1:n_trials);

surround_time = [-1,2];
surround_sample_rate = 100;
surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);
pull_times = align_times + surround_time_points;

[wheel_velocity,wheel_move] = ...
    AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);
wheel_aligned = interp1(timelite.timestamps,wheel_velocity,pull_times);

wheel_aligned_avg = ap.groupfun(@nanmean,wheel_aligned,trial_stim);

figure('Name',rec_day);
plot(surround_time_points,wheel_aligned_avg'); drawnow;

xlim([-0.2,1]);
ap.scalebar(0.2,200);
axis off
ap.prettyfig;


%% Fig 3D: SNr example neurons

animal = 'DS023';
rec_day = '2025-12-22';

plot_units = [38,25];
stim_colors = {[0.7,0,0],lines(1)};
figure;
h = tiledlayout(3,length(plot_units),'TileSpacing','tight');

for curr_stim = 1:2
    switch curr_stim
        case 1
            rec_time = '1507'; % checkerboard > R
        case 2
            rec_time = '1514'; % grating > L
    end
    ap.load_recording;

    % Quiescent trials
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        (1:length(stimOn_times))');

    % Stim-aligned PSTH
    stim_x = vertcat(trial_events.values.TrialStimX);
    stim_x = stim_x(1:length(stimOn_times));

    [~,spike_group] = ismember(spike_templates,plot_units);

    [psth,raster,psth_t] = ap.psth(spike_times_timelite, ...
        stimOn_times(quiescent_trials & stim_x == 0),spike_group,...
        'window',[-0.3,1],'smoothing',100);

    % Plot
    for curr_unit = 1:length(plot_units)
        % (psth)
        nexttile(h,tilenum(h,1,curr_unit)); hold on;
        plot(psth_t,psth(curr_unit,:),'linewidth',2,'color',stim_colors{curr_stim});

        % (raster)
        nexttile(h,tilenum(h,1+curr_stim,curr_unit));
        [raster_y,raster_x] = find(raster(:,:,curr_unit));
        scatter(psth_t(raster_x),raster_y,10,stim_colors{curr_stim}, ...
            '|','linewidth',1);
        axis off
        set(gca,'YDir','reverse');
    end
end
linkaxes(h.Children,'x');

ap.prettyfig;
xlim([-0.2,1]);


%% Fig 4B: VM stim responses

conditions = ["naive","trained"];

mua = cell(size(conditions));
for condition = conditions

    switch condition
        case "naive"
            animal = 'AM010';
            rec_day = '2023-12-08';
            rec_time = '1517';
            use_depths = [3100,3500];
        case "trained"
            animal = 'AP009';
            rec_day = '2023-07-12';
            rec_time = '1458';
            use_depths = [2700,3200];
    end

    load_parts.ephys = true;
    ap.load_recording;

    % Quiescent trials
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        (1:length(stimOn_times))');

    % Stim-aligned PSTH
    stim_x = vertcat(trial_events.values.TrialStimX);

    use_spikes = isbetween(spike_depths,use_depths(1),use_depths(2));

    [mua_psth,~,mua_psth_t] = ap.psth(spike_times_timelite(use_spikes), ...
        stimOn_times(quiescent_trials & stim_x == 90), ...
        'window',[-0.2,1],'smoothing',50);

    mua{ismember(conditions,condition)} = mua_psth;

end

baseline_t = mua_psth_t < 0;
mua_norm = cellfun(@(x) (x-mean(x(baseline_t),2))./mean(x(baseline_t),2),mua,'uni',false);
figure;plot(mua_psth_t,cell2mat(mua_norm')','linewidth',2);
xlim([-0.2,0.7]);
ap.prettyfig;
ap.scalebar(0.2,1);
axis off;
xline([0,0.5]);


%% Fig 5A-B: Cortex-striatum regression

animal = 'AP025';
rec_day = '2024-09-10';
rec_time = '1802';

ap.load_recording;

% Set timing/sampling
upsample_factor = 1;
sample_rate = (1/mean(diff(wf_t)))*upsample_factor;
skip_seconds = 60;
time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Bin spikes
depth_group_edges = [1100,1600];

n_depths = length(depth_group_edges) - 1;
depth_group = discretize(spike_depths,depth_group_edges);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timelite(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

% Regress cortex to spikes
use_svs = 1:200;
kernel_t = [-0.5,0.5];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 20;
zs = [false,false];
cvfold = 1;
return_constant = false;
use_constant = true;

fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

[k,predicted_spikes,explained_var] = ...
    ap.regresskernel(fVdf_deconv_resample, ...
    binned_spikes_std,kernel_frames, ...
    lambda,zs,cvfold,return_constant,use_constant);

% Convert kernel to pixel space
k_px = plab.wf.svd2px(wf_U(:,:,use_svs),k);

figure;imagesc(nanmean(k_px(:,:,15:21),3))
clim([0,0.002])
colormap(AP_colormap('WG',[],1));
ap.wf_draw('ccf',[0.5,0.5,0.5]);
axis image off;
ap.prettyfig;

% Hand-draw ROIs to plot timecourse
ap.imscroll(k_px); axis image;
figure; hold on

% (run this line for mPFC and PPC ROIs)
plot(-kernel_frames/sample_rate,roi.trace,'linewidth',2);

xline(0);
ap.prettyfig;


%% Fig 5C-D: Cortex-VM regression

animal = 'AM025';
rec_day = '2024-04-25';
rec_time = '1322';

ap.load_recording;

% Set timing/sampling
upsample_factor = 1;
sample_rate = (1/mean(diff(wf_t)))*upsample_factor;
skip_seconds = 60;
time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Bin spikes
depth_group_edges = [2800,3100];

n_depths = length(depth_group_edges) - 1;
depth_group = discretize(spike_depths,depth_group_edges);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timelite(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

% Regress cortex to spikes
use_svs = 1:200;
kernel_t = [-0.5,0.5];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 20;
zs = [false,false];
cvfold = 1;
return_constant = false;
use_constant = true;

fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

[k,predicted_spikes,explained_var] = ...
    ap.regresskernel(fVdf_deconv_resample, ...
    binned_spikes_std,kernel_frames, ...
    lambda,zs,cvfold,return_constant,use_constant);

% Convert kernel to pixel space
k_px = plab.wf.svd2px(wf_U(:,:,use_svs),k);

figure;imagesc(nanmean(k_px(:,:,15:21),3))
clim([0,0.002])
colormap(AP_colormap('WG',[],1));
ap.wf_draw('ccf',[0.5,0.5,0.5]);
axis image off;
ap.prettyfig;

% Hand-draw ROIs to plot timecourse
ap.imscroll(k_px); axis image;

figure; hold on
% (run this line for mPFC and PPC ROIs)
plot(-kernel_frames/sample_rate,roi.trace,'linewidth',2);

xline(0);
ap.prettyfig;



