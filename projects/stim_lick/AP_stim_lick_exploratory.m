
%% Lick raster

[lick_psth,lick_raster,lick_t] = ap.psth(lick_times,stimOn_times,...
    'window',[-4,4],'bin_size',0.01,'smoothing',10);

[lick_trial,lick_t_raster_idx] = find(lick_raster);
lick_t_raster = lick_t(lick_t_raster_idx);

figure;
h = tiledlayout(4,1);

nexttile(1)
plot(lick_t,lick_psth,'k');
xline(0,'r');

nexttile([3,1])
plot(lick_t_raster,lick_trial,'.k');
xline(0,'r');
set(gca,'YDir','reverse');

linkaxes(h.Children,'x');


%% Lick raster (stim move conditioning)

n_trials = length([trial_events.timestamps.Outcome]);

stimOn_times = photodiode_on_times(1:2:end);
stim_move_times = photodiode_off_times(1:2:end);
stim_center_times = photodiode_off_times(2:2:end);

[lick_psth,lick_raster,lick_t] = ap.psth(lick_times,reward_times,...
    'window',[-0.5,1],'bin_size',0.01,'smoothing',10);

figure('name',sprintf('%s %s',animal,rec_day));
h = tiledlayout(4,1);

nexttile(1)
plot(lick_t,lick_psth,'k');
xline(0,'color',[0.7,0.7,0]);

nexttile([3,1]);

% Trial params
trial_static_stim_time = vertcat(trial_events.values(1:n_trials).TrialStimStaticTime);
trial_quiescence_time = vertcat(trial_events.values(1:n_trials).TrialQuiescence);

% Sort by: 
% % trial order 
% sort_idx = 1:n_trials;
% 
% stim static time
[~,sort_idx] = sort(trial_static_stim_time);
% 
% % reward time
% [~,sort_idx] = sort(reward_times - stimOn_times(1:n_trials));
%
% % quiescence time
% [~,sort_idx] = sort(trial_quiescence_time);


[lick_trial,lick_t_raster_idx] = find(lick_raster(sort_idx,:));
lick_t_raster = lick_t(lick_t_raster_idx);

plot(lick_t_raster,lick_trial,'.k');
% hold on;
% plot(-trial_quiescence_time(sort_idx),1:n_trials,'m')
% plot(stim_move_times(sort_idx)-stimOn_times(sort_idx),1:n_trials,'color',[0,0.8,0])
% plot(stim_center_times(sort_idx)-stimOn_times(sort_idx),1:n_trials,'r')
% plot(reward_times(sort_idx)-stimOn_times(sort_idx),1:n_trials,'c')


xline(0,'color',[0.7,0.7,0]);
set(gca,'YDir','reverse');

linkaxes(h.Children,'x');

% figure;plot(reward_times - stim_move_times(1:n_trials) - 1,'.k')
% ylabel('Stim to reward time')

% (lick autocorrelation)
% lick_xcorr_timebin = 0.01;
% [lick_xcorr,lick_lags] = xcorr(histcounts(lick_times,lick_times(1): ...
%     lick_xcorr_timebin:lick_times(end)),10/lick_xcorr_timebin,'coeff');
% figure;plot(lick_lags*lick_xcorr_timebin,lick_xcorr,'k');
% xlabel('Lag');
% ylabel('Lick autocorr');


%% stats (making variant of the stim wheel conditional resampling)

stimOn_times = photodiode_on_times(1:2:end);


% Only look at completed trials
n_trials = length([trial_events.timestamps.Outcome]);

% Get quiescence range from Bonsai parameters
quiescence_range = trial_events.parameters.QuiescenceTimes;

% Loop through trials, get possible valid stim times
stimOn_times_valid = cell(n_trials,1);
for curr_trial = 1:n_trials

    % NOTE: this is only using alternate quiescence periods, not alternate
    % ITIs. No quiescence clock during ITI means no direct Bonsai measure
    % of would-be quiescence resets, and they're not accurately estimatable
    % from NIDAQ because of timing/precision differences. This means that
    % there's fewer "valid" stim times, so less statistical power, and may
    % err on the side of missing learned days.

    % Get time of last quiescence reset (Timelite)
    % (define as stim time - trial quiescence time)
    last_quiescence_reset = stimOn_times(curr_trial) - trial_events.values(curr_trial).TrialQuiescence;

    % Get quiescence reset times (Timelite)
    curr_quiescence_resets_bonsai = vertcat(...
        trial_events.timestamps(curr_trial).QuiescenceStart, ... % start of quiescence period
        trial_events.timestamps(curr_trial).QuiescenceReset);    % all quiescence resets

    curr_quiescence_resets_timelite = last_quiescence_reset + ...
        seconds(curr_quiescence_resets_bonsai - ...
        curr_quiescence_resets_bonsai(end,:));

    % Get quiescence durations between resets, including post-stim lick
    curr_poststim_lick = lick_times(find(lick_times > stimOn_times(curr_trial),1));

    curr_quiescence_durations = ...
        diff(vertcat(curr_quiescence_resets_timelite, ...
        curr_poststim_lick));

    % Get valid quiescence times that would yield same response movement
    % (i.e. last quiescence duration was first over-threshold)
    quiescence_overthresh_grid = curr_quiescence_durations > quiescence_range';
    valid_quiescence_times = quiescence_range( ...
        (sum(quiescence_overthresh_grid,1) == 1) & quiescence_overthresh_grid(end,:));

    %%%%%%% TESTING: THIS IS DIFFERENT NOW

    % For each quiescence threshold, find first instance it was satisfied
    % (this is the time the stim would have appeared, given that q)
    [q_hit,q_first_idx] = max(quiescence_overthresh_grid,[],1);

    % Get quiescence offsets relative to stim time
    curr_quiescence_relative_to_stim = curr_quiescence_resets_timelite - stimOn_times(curr_trial);

    % Get valid stim times: for given quiescence threshold, first instance
    % where that threshold was met
    stimOn_offsets_valid = ...
        reshape(curr_quiescence_relative_to_stim(q_first_idx(q_hit)),[],1) + ...
        quiescence_range(q_hit);
    stimOn_times_valid{curr_trial} = stimOn_times(curr_trial) + stimOn_offsets_valid;

    %%%%%%%%%%%%%

    % (debugging: plot trial)
    if false

        trial_start_t = min(curr_quiescence_resets_timelite) - 2;
        trial_end_t = stimOn_times(curr_trial) + 4;
        curr_t_idx = timelite.timestamps > trial_start_t & timelite.timestamps < trial_end_t;

        figure;
        % (wheel velocity)
        plot(timelite.timestamps(curr_t_idx),wheel_velocity(curr_t_idx));
        % (quiescence resets)
        xline(curr_quiescence_resets_timelite,'c')
        % (stim onset)
        xline(stimOn_times(curr_trial),'g','linewidth',2);
        % (move onset)
        xline(stimOn_times(curr_trial)+stim_to_move(curr_trial),'r','linewidth',2);
        % (possible stim times)
        xline(stimOn_times_valid{curr_trial},'--m');
        % (last quiescence reset)
        xline(stimOn_times(curr_trial)-trial_events.values(curr_trial).TrialQuiescence,'k')

    end
end

% Count windowed licks after real stim
lick_window_t = 0.5;
stim_lick_window = stimOn_times + [0,lick_window_t];

lick_hist_full = histcounts(lick_times,reshape(stim_lick_window',[],1));
lick_stim = mean(lick_hist_full(1:2:end));

% Count windowed licks after alternate stim

n_sample = 1000;
lick_hist_null = nan(n_sample,1);
for curr_sample = 1:n_sample
    curr_shuff_stim_times = cellfun(@(x) datasample(x,1),stimOn_times_valid);
    curr_shuff_stim_lick_window = curr_shuff_stim_times + [0,lick_window_t];

    curr_shuff_lick_hist_full = histcounts(lick_times,reshape(curr_shuff_stim_lick_window',[],1));
    lick_hist_null(curr_sample) = mean(curr_shuff_lick_hist_full(1:2:end));
    ap.print_progress_fraction(curr_sample,n_sample);
end

% Plot null and measured
figure;
histogram(lick_hist_null)
xline(lick_stim,'r');



%%% probably need to change this a bit: at the moment, the averaging time
%%% can run into the actual reward time if there's enough of a difference
%%% (e.g. if actual quiescence was 0.2 and alt was 1.8, +0.5s from actual
%%% might run into real reward time)


%% Lick raster (two-stim move)

% fix issue with DS017 2025-02-18
if ~exist('fix_flag','var');fix_flag = false;end
if strcmp(animal,'DS017') && strcmp(rec_day,'2025-02-18') && ~fix_flag
    photodiode_on_times(64) = [];
    photodiode_off_times(64) = [];
    fix_flag = true;
end

rewarded_x = 90;%trial_events.parameters.RewardedX;

n_trials = sum(cellfun(@(x) length(x) == 2,{trial_events.timestamps.StimOn}));

stim_x = vertcat(trial_events.values(1:n_trials).TrialX);

stim_pd_n = (stim_x==-90)*1 + (stim_x==rewarded_x)*2; % (for center-out): (stim_x==-90)*2 + (stim_x==rewarded_x)*2;
    
stim_pd_on_grouped = mat2cell(photodiode_on_times(1:sum(stim_pd_n)),stim_pd_n);
stim_pd_off_grouped = mat2cell(photodiode_off_times(1:sum(stim_pd_n)),stim_pd_n);

% 2 PD ups for right (on, move), 1 for left (on)
stimOn_times = cellfun(@(x) x(1), stim_pd_on_grouped);

stim_move_times = nan(size(stim_x));
stim_move_times(stim_x == rewarded_x) = cellfun(@(x) x(1), stim_pd_off_grouped(stim_x == rewarded_x));

% (stim center times, or would-be for CS-)
stim_center_times = nan(size(stim_x));
stim_center_times(stim_x == rewarded_x) = cellfun(@(x) x(2), stim_pd_off_grouped(stim_x == rewarded_x));
stim_center_times(stim_x == -90) = cell2mat(stim_pd_off_grouped(stim_x == -90))+1; % (for center-out): cellfun(@(x) x(2), stim_pd_off_grouped(stim_x == -90));

% (reward times, or would-be for CS-)
stim_reward_times = nan(size(stim_x));
% stim_reward_times(stim_x == rewarded_x) = reward_times;
stim_reward_times(stim_x == rewarded_x) = interp1(lick_times,lick_times,stim_center_times(stim_x == rewarded_x),'next');
stim_reward_times(stim_x == -90) = interp1(lick_times,lick_times,stim_center_times(stim_x == -90),'next');

% Trial params
trial_static_stim_time = vertcat(trial_events.values(1:n_trials).TrialStimStaticTime);
trial_quiescence_time = vertcat(trial_events.values(1:n_trials).TrialQuiescence);

% Plot licks aligned to stim onset / (would-be) center
figure('name',sprintf('%s %s',animal,rec_day));
h = tiledlayout(4,3,'TileIndexing','ColumnMajor');

for curr_align = 1:3
    switch curr_align
        case 1
            plot_align = stimOn_times;
        case 2
            plot_align = stim_center_times;
        case 3
            plot_align = stim_reward_times;
    end

    [lick_psth_r,lick_raster_r,lick_t] = ap.psth(lick_times,plot_align(stim_x == rewarded_x & ~isnan(plot_align)),...
        'window',[-7,10],'bin_size',0.01,'smoothing',20);
    [lick_psth_l,lick_raster_l] = ap.psth(lick_times,plot_align(stim_x == -90 & ~isnan(plot_align)),...
        'window',[-7,10],'bin_size',0.01,'smoothing',20);

    nexttile(h,tilenum(h,1,curr_align)); hold on;
    plot(lick_t,lick_psth_r,'r');
    plot(lick_t,lick_psth_l,'b');
    xline(0,'color',[0.7,0.7,0]);

    nexttile(h,tilenum(h,2,curr_align),[3,1]); hold on;

    [lick_trial_r,lick_t_raster_r_idx] = find(lick_raster_r);
    [lick_trial_l,lick_t_raster_l_idx] = find(lick_raster_l);

    lick_t_raster_r = lick_t(lick_t_raster_r_idx);
    lick_t_raster_l = lick_t(lick_t_raster_l_idx);

    plot(lick_t_raster_r,lick_trial_r,'.r');
    plot(lick_t_raster_l,lick_trial_l+size(lick_raster_r,1),'.b');
    hold on; set(gca,'YDir','reverse');
    xline(0,'color',[0.7,0.7,0]);
    axis tight

end

linkaxes(h.Children,'x');

%% Lick raster (two-stim static)

% rewarded_x = trial_events.parameters.RewardedFrequency;
rewarded_x = trial_events.parameters.RewardedX;

n_trials = sum(cellfun(@(x) length(x) == 2,{trial_events.timestamps.StimOn}));

% stim_x = vertcat(trial_events.values(1:n_trials).TrialFrequency);
stim_x = vertcat(trial_events.values(1:n_trials).TrialX);
stimOn_times = photodiode_on_times(1:n_trials);
stimOff_times = photodiode_off_times(1:n_trials);

% (reward times, or would-be for CS-)
stim_reward_times = nan(size(stim_x));
stim_reward_times(stim_x == rewarded_x) = reward_times;
stim_reward_times(stim_x ~= rewarded_x) = interp1(lick_times,lick_times,stimOn_times(stim_x ~= rewarded_x),'next');

% Trial params
trial_static_stim_time = vertcat(trial_events.values(1:n_trials).TrialStimStaticTime);
trial_quiescence_time = vertcat(trial_events.values(1:n_trials).TrialQuiescence);

% Plot licks aligned to stim onset / (would-be) center
% sort_val = 1:n_trials;
sort_val = trial_static_stim_time;
% sort_val = trial_quiescence_time;

figure('name',sprintf('%s %s',animal,rec_day));
h = tiledlayout(4,2,'TileIndexing','ColumnMajor');

for curr_align = 1:2
    switch curr_align
        case 1
            plot_align = stimOn_times;
        case 2
            plot_align = stim_reward_times;
    end

    r_trials = find(stim_x == rewarded_x & ~isnan(plot_align));
    l_trials = find(stim_x ~= rewarded_x & ~isnan(plot_align));

    [lick_psth_r,lick_raster_r,lick_t] = ap.psth(lick_times,plot_align(r_trials),...
        'window',[-7,10],'bin_size',0.01,'smoothing',20);
    [lick_psth_l,lick_raster_l] = ap.psth(lick_times,plot_align(l_trials),...
        'window',[-7,10],'bin_size',0.01,'smoothing',20);

    nexttile(h,tilenum(h,1,curr_align)); hold on;
    plot(lick_t,lick_psth_r,'r');
    plot(lick_t,lick_psth_l,'b');
    xline(0,'color',[0.7,0.7,0]);

    nexttile(h,tilenum(h,2,curr_align),[3,1]); hold on;

    [~,sort_r] = sort(sort_val(r_trials));
    [lick_trial_r,lick_t_raster_r_idx] = find(lick_raster_r(sort_r,:));

    [~,sort_l] = sort(sort_val(l_trials));
    [lick_trial_l,lick_t_raster_l_idx] = find(lick_raster_l(sort_l,:));

    lick_t_raster_r = lick_t(lick_t_raster_r_idx);
    lick_t_raster_l = lick_t(lick_t_raster_l_idx);

    plot(lick_t_raster_r,lick_trial_r,'.r');
    plot(trial_static_stim_time(r_trials(sort_r)),1:length(r_trials),'k','linewidth',2)

    plot(lick_t_raster_l,lick_trial_l+size(lick_raster_r,1),'.b');
    hold on; set(gca,'YDir','reverse');
    xline(0,'color',[0.7,0.7,0]);
    axis tight

end

linkaxes(h.Children,'x');

%% Task kernel (move task)

rewarded_x = trial_events.parameters.RewardedX;

n_trials = sum(cellfun(@(x) length(x) == 2,{trial_events.timestamps.StimOn}));

stim_x = vertcat(trial_events.values(1:n_trials).TrialX);

stim_pd_n = (stim_x==-90)*1 + (stim_x==rewarded_x)*2;
    
stim_pd_on_grouped = mat2cell(photodiode_on_times(1:sum(stim_pd_n)),stim_pd_n);
stim_pd_off_grouped = mat2cell(photodiode_off_times(1:sum(stim_pd_n)),stim_pd_n);

% 2 PD ups for right (on, move), 1 for left (on)
stimOn_times = cellfun(@(x) x(1), stim_pd_on_grouped);

stim_move_times = nan(size(stim_x));
stim_move_times(stim_x == rewarded_x) = cellfun(@(x) x(1), stim_pd_off_grouped(stim_x == rewarded_x));

% (stim center times, or would-be for -90 stim)
stim_center_times = nan(size(stim_x));
stim_center_times(stim_x == rewarded_x) = cellfun(@(x) x(2), stim_pd_off_grouped(stim_x == rewarded_x));
stim_center_times(stim_x == -90) = cell2mat(stim_pd_off_grouped(stim_x == -90))+1;


time_bins = [wf_t;wf_t(end)+1/wf_framerate];

task_events = zeros(0,length(time_bins)-1);
task_events(end+1,:) = histcounts(stimOn_times(stim_x == rewarded_x),time_bins);
task_events(end+1,:) = histcounts(stim_center_times(stim_x == rewarded_x),time_bins);
task_events(end+1,:) = histcounts(stimOn_times(stim_x == -90),time_bins);

n_components = 200;

frame_shifts = -5:20;
lambda = 20;
cv_fold = 3;

skip_t = 60; % seconds start/end to skip for artifacts
skip_frames = round(skip_t*wf_framerate);
[kernels,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
    task_events(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

ap.imscroll(plab.wf.svd2px(wf_U(:,:,1:size(kernels,1)),kernels));
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image;


%% Task kernel (static task)

rewarded_x = trial_events.parameters.RewardedX;
% rewarded_x = trial_events.parameters.RewardedFrequency;

n_trials = sum(cellfun(@(x) length(x) == 2,{trial_events.timestamps.StimOn}));

stim_x = vertcat(trial_events.values(1:n_trials).TrialX);
% stim_x = vertcat(trial_events.values(1:n_trials).TrialFrequency);
stimOn_times = photodiode_on_times(1:n_trials);
stimOff_times = photodiode_off_times(1:n_trials);

time_bins = [wf_t;wf_t(end)+1/wf_framerate];

task_events = zeros(0,length(time_bins)-1);
task_events(end+1,:) = histcounts(stimOn_times(stim_x == rewarded_x),time_bins);
task_events(end+1,:) = histcounts(stimOn_times(stim_x ~= rewarded_x),time_bins);

n_components = 200;

frame_shifts = -5:20;
lambda = 20;
cv_fold = 3;

skip_t = 60; % seconds start/end to skip for artifacts
skip_frames = round(skip_t*wf_framerate);
[kernels,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
    task_events(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

kernels_px = plab.wf.svd2px(wf_U(:,:,1:size(kernels,1)),kernels);
ap.imscroll(kernels_px,frame_shifts);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image;

t = frame_shifts/wf_framerate;
cs_minus_color = ap.colormap('WB');
cs_plus_color = ap.colormap('WR');

stim_t = t > 0 & t < 0.2;
kernels_px_max = squeeze(max(kernels_px(:,:,stim_t,:),[],3));

col_lim = [0,1e-4];

figure;
h = tiledlayout(1,2,'TileSpacing','none');

nexttile; imagesc(kernels_px_max(:,:,1)); 
clim(col_lim); axis image off;
colormap(gca,cs_plus_color);
ap.wf_draw('ccf','k');

nexttile; imagesc(kernels_px_max(:,:,2)); 
clim(col_lim); axis image off;
colormap(gca,cs_minus_color);
ap.wf_draw('ccf','k');

%% Task kernel (wheel task)

time_bins = [wf_t;wf_t(end)+1/wf_framerate];

stim_regressors = histcounts(stimOn_times,time_bins);

n_components = 200;

frame_shifts = -10:20;
lambda = 20;
cv_fold = 3;

skip_t = 60; % seconds start/end to skip for artifacts
skip_frames = round(skip_t*wf_framerate);
[kernels,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
    stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

ap.imscroll(plab.wf.svd2px(wf_U(:,:,1:n_components),kernels))
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image;







