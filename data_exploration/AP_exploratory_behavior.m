%% Exploratory behavior analysis


%% Align mousecam to event

use_cam = mousecam_fn;
use_t = mousecam_times;

align_times = stimOn_times;
align_category = vertcat(trial_events.values.TrialStimX);

% (get only quiescent trials)
[wheel_velocity,wheel_move] = AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);
framerate = 30;
wheel_window = [0,0.5];
wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
wheel_window_t_peri_event = align_times_all + wheel_window_t;
event_aligned_move = interp1(timelite.timestamps, ...
    +wheel_move,wheel_window_t_peri_event,'previous');
quiescent_trials = ~any(event_aligned_move,2);

use_align = align_times(align_category == 90 & quiescent_trials);

surround_frames = 30;

% Initialize video reader, get average and average difference
vr = VideoReader(use_cam);
cam_im1 = read(vr,1);

cam_align_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    surround_frames*2+1);
cam_align_diff_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    surround_frames*2);

frame_t_offset = nan(size(use_align));
for curr_align = 1:length(use_align)

    % Find closest camera frame to timepoint
    [frame_t_offset(curr_align),curr_frame] = ...
        min(abs(use_align(curr_align) - use_t));

    % Pull surrounding frames
    curr_surround_frames = curr_frame + [-surround_frames,surround_frames];
    curr_clip = double(squeeze(read(vr,curr_surround_frames)));
    curr_clip_diff = abs(diff(curr_clip,[],3));

    cam_align_avg = cam_align_avg + curr_clip./length(use_align);
    cam_align_diff_avg = cam_align_diff_avg + curr_clip_diff./length(use_align);

    AP_print_progress_fraction(curr_align,length(use_align));
end

surround_t = [-surround_frames:surround_frames]./vr.FrameRate;
AP_imscroll(cam_align_avg,surround_t)
axis image;

surround_t = [-surround_frames:surround_frames]./vr.FrameRate;
AP_imscroll(cam_align_diff_avg,surround_t(2:end))
axis image;

% Plot difference within window
use_t = [0,0.2];
use_t_idx = surround_t >= use_t(1) & surround_t <= use_t(2);
figure;
imagesc(nanmean(cam_align_diff_avg(:,:,use_t_idx(2:end)),3));
axis image off;


%% Align wheel to event

align_times = photodiode_times(1:2:end);

surround_time = [-10,10];
surround_sample_rate = 100;
surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);
pull_times = align_times + surround_time_points;

[wheel_velocity,wheel_move] = ...
    AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);

event_aligned_wheel_vel = interp1(timelite.timestamps, ...
    wheel_velocity,pull_times);
event_aligned_wheel_move = interp1(timelite.timestamps, ...
    +wheel_move,pull_times,'previous');

figure;
subplot(2,2,1);
imagesc(surround_time_points,[],event_aligned_wheel_vel)
xline(0,'color','r');
clim(max(abs(clim)).*[-1,1])
colormap(gca,AP_colormap('BWR'));

subplot(2,2,3);
plot(surround_time_points,nanmean(event_aligned_wheel_vel,1));
xline(0,'color','r');

subplot(2,2,2);
imagesc(surround_time_points,[],event_aligned_wheel_move)
xline(0,'color','r');
ylabel('Velocity');
xlabel('Time from event');

subplot(2,2,4);
plot(surround_time_points,nanmean(event_aligned_wheel_move,1));
xline(0,'color','r');
ylabel('Move prob.');
xlabel('Time from event');


[~,sort_idx] = sort([trial_events.values.TrialQuiescence]);
figure;
imagesc(surround_time_points,[],event_aligned_wheel_move(sort_idx,:));
xline(0,'color','r','linewidth',2);
hold on;
plot(-[trial_events.values(sort_idx).TrialQuiescence],1:length(trial_events.values),'b','linewidth',2);


%% Behavior across days

animal = 'AP009';
use_workflow = {'stim_wheel_right_stage1','stim_wheel_right_stage2'};
recordings = ap.find_recordings(animal,[],use_workflow);

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
    rec_time = recordings(curr_recording).protocol{end};
    load_parts.widefield = false;
    ap.load_experiment;

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
    AP_print_progress_fraction(curr_recording,length(recordings));

end

learned_day = rxn_stat_p < 0.05;

relative_day = days(datetime({recordings.day}) - datetime({recordings(1).day}))+1;
nonrecorded_day = setdiff(1:length(recordings),relative_day);

figure('Name',animal);
t = tiledlayout('flow');
title(t,animal);

nexttile;
yyaxis left; plot(relative_day,n_trials_water(:,1));
ylabel('# trials');
yyaxis right; plot(relative_day,frac_move_day);
ylabel('Fraction time moving');
xlabel('Day');
if any(nonrecorded_day)
    xline(nonrecorded_day,'--k');
end
xline(relative_day(learned_day),'g')

nexttile;
yyaxis left
plot(relative_day,rxn_med)
set(gca,'YScale','log');
ylabel('Med. rxn');
xlabel('Day');
if any(nonrecorded_day)
    xline(nonrecorded_day,'--k');
end
xline(relative_day(learned_day),'g')

yyaxis right
prestim_max = max(frac_move_stimalign(:,surround_time_points < 0),[],2);
poststim_max = max(frac_move_stimalign(:,surround_time_points > 0),[],2);
plot(relative_day,(poststim_max-prestim_max)./(poststim_max+prestim_max));
yline(0);
ylabel('pre/post move idx');
xlabel('Day');

nexttile;
imagesc(surround_time_points,[],frac_move_stimalign); hold on;
clim([0,1]);
colormap(gca,AP_colormap('WK'));
set(gca,'YTick',1:length(recordings),'YTickLabel',string(datetime({recordings.day},'format','MM-dd')));
xlabel('Time from stim');
plot(0,find(learned_day),'.g')

nexttile; hold on
AP_errorfill(surround_time_points,frac_move_stimalign(learned_day,:)', ...
    0.02,[0,1,0],1,false);
set(gca,'ColorOrder',copper(length(recordings)));
plot(surround_time_points,frac_move_stimalign','linewidth',2);
xline(0,'color','k');
ylabel('Fraction moving');
xlabel('Time from stim');





















