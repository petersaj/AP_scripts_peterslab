
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


%% Lick raster (two-stim move conditioning)

% fix issue with DS017 2025-02-18
if ~exist('fix_flag','var');fix_flag = false;end
if strcmp(animal,'DS017') && strcmp(rec_day,'2025-02-18') && ~fix_flag
    photodiode_on_times(64) = [];
    photodiode_off_times(64) = [];
    fix_flag = true;
end

n_trials = sum(cellfun(@(x) length(x) == 2,{trial_events.timestamps.StimOn}));

stim_x = vertcat(trial_events.values(1:n_trials).TrialX);

use_pd_times = sum((stim_x==-90)*1 + (stim_x==90)*2);
    
stim_pd_on_grouped = mat2cell(photodiode_on_times(1:use_pd_times),(stim_x==-90)*1 + (stim_x==90)*2);
stim_pd_off_grouped = mat2cell(photodiode_off_times(1:use_pd_times),(stim_x==-90)*1 + (stim_x==90)*2);

% 2 PD ups for right (on, move), 1 for left (on)
stimOn_times = cellfun(@(x) x(1), stim_pd_on_grouped);

stim_move_times = nan(size(stim_x));
stim_move_times(stim_x == 90) = cellfun(@(x) x(1), stim_pd_off_grouped(stim_x == 90));

% (stim center times, or would-be for -90 stim)
stim_center_times = nan(size(stim_x));
stim_center_times(stim_x == 90) = cellfun(@(x) x(2), stim_pd_off_grouped(stim_x == 90));
stim_center_times(stim_x == -90) = cell2mat(stim_pd_off_grouped(stim_x == -90))+1;

use_trials = find(stim_x == 90);

[lick_psth,lick_raster,lick_t] = ap.psth(lick_times,stim_center_times(use_trials),...
    'window',[-7,10],'bin_size',0.01,'smoothing',20);

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
% trial order 
sort_idx = 1:length(use_trials);
%
% stim staix] = sort(trial_static_stim_time(use_trials));
% 
% % reward time
% [~,sort_idx] = sort(reward_times(use_trials) - stimOn_times(use_trials));
%
% % quiescence time
% [~,sort_idx] = sort(trial_quiescence_time(use_trials));


[lick_trial,lick_t_raster_idx] = find(lick_raster(sort_idx,:));
lick_t_raster = lick_t(lick_t_raster_idx);

plot(lick_t_raster,lick_trial,'.k');
hold on; set(gca,'YDir','reverse');
plot(-trial_quiescence_time(use_trials(sort_idx)),1:length(use_trials),'m')
xline(0,'color',[0.7,0.7,0]);

% plot(stim_move_times(use_trials(sort_idx))-stimOn_times(use_trials(sort_idx)),1:length(use_trials),'color',[0,0.8,0])
% plot(stim_center_times(use_trials(sort_idx))-stimOn_times(use_trials(sort_idx)),1:length(use_trials),'r')
% plot(reward_times(sort_idx)-stimOn_times(use_trials(sort_idx)),1:length(use_trials),'c')

linkaxes(h.Children,'x');





