%% Exploratory behavior analysis


%% Align mousecam to event

use_cam = mousecam_fn;
use_t = mousecam_times;

% (passive)
stim_window = [0,0.5];
quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
    timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
    timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
    1:length(stimOn_times))';

stim_x = vertcat(trial_events.values.StimFrequence);
use_align = stimOn_times(stim_x == 8000 & quiescent_trials);

% % (task)
% use_align = stimOn_times;


surround_frames = 60;

% Initialize video reader, get average and average difference
vr = VideoReader(use_cam);
cam_im1 = read(vr,1);

cam_align_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    surround_frames*2+1);
cam_align_diff_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    surround_frames*2);

frame_t_offset = nan(size(use_align));
for curr_align = 2:length(use_align)

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

%% Align mousecam ROI to event

use_cam = mousecam_fn;
use_t = mousecam_times;

% (passive)
stim_window = [0,0.5];
quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
    timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
    timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
    1:length(stimOn_times))';

stim_x = vertcat(trial_events.values.TrialStimX);
use_align = stimOn_times(stim_x == 90 & quiescent_trials);

% (task)
% use_align = stimOn_times;

surround_frames = 60;

% Initialize video reader, get average and average difference
vr = VideoReader(use_cam);
cam_im1 = read(vr,1);

% Draw ROI
h = figure;imagesc(cam_im1);axis image; 
roi_mask = roipoly;
close(h);

cam_roi_diff_align = nan(length(use_align),surround_frames*2);

% (would probably be way faster and reasonable to just load in the entire
% movie?)
for curr_align = 1:length(use_align)

    % Find closest camera frame to timepoint
    curr_frame = interp1(mousecam_times,1:length(mousecam_times), ...
        use_align(curr_align),'nearest');

    % Pull surrounding frames
    curr_surround_frames = curr_frame + [-surround_frames,surround_frames];
    if any(curr_surround_frames < 0) || any(curr_surround_frames > vr.NumFrames)
        continue
    end

    curr_clip_diff_flat = reshape(abs(diff(double(squeeze( ...
        read(vr,curr_surround_frames))),[],3)),[],surround_frames*2);

    cam_roi_diff_align(curr_align,:) = ...
        ((roi_mask(:))'*curr_clip_diff_flat)./sum(roi_mask,'all');

    AP_print_progress_fraction(curr_align,length(use_align));
end


surround_t = [-surround_frames:surround_frames]./vr.FrameRate;

figure;imagesc(surround_t(2:end),[],cam_roi_diff_align);
figure; hold on;
plot(surround_t(2:end),nanmean(cam_roi_diff_align,1));
plot(surround_t(2:end),nanmedian(cam_roi_diff_align,1));


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

animals = {'AP020','AP021','AP022'};

% Create master tiled layout
figure;
t = tiledlayout(1,length(animals),'TileSpacing','tight');

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
        AP_print_progress_fraction(curr_recording,length(recordings));

    end

    % Define learned day from reaction stat p-value and reaction time
    learned_day = rxn_stat_p < 0.05 & rxn_med < 2;

    relative_day = days(datetime({recordings.day}) - datetime({recordings(1).day}))+1;
    nonrecorded_day = setdiff(1:length(recordings),relative_day);

    % Draw in tiled layout nested in master
    t_animal = tiledlayout(t,4,1);
    t_animal.Layout.Tile = curr_animal_idx;
    title(t_animal,animal);

    nexttile(t_animal);
    yyaxis left; plot(relative_day,n_trials_water(:,1));
    ylabel('# trials');
    yyaxis right; plot(relative_day,frac_move_day);
    ylabel('Fraction time moving');
    xlabel('Day');
    if any(nonrecorded_day)
        xline(nonrecorded_day,'--k');
    end
    if any(learned_day)
        xline(relative_day(learned_day),'g');
    end

    nexttile(t_animal);
    yyaxis left
    plot(relative_day,rxn_med)
    set(gca,'YScale','log');
    ylabel('Med. rxn');
    xlabel('Day');
    if any(nonrecorded_day)
        xline(nonrecorded_day,'--k');
    end
    if any(learned_day)
        xline(relative_day(learned_day),'g');
    end

    yyaxis right
    prestim_max = max(frac_move_stimalign(:,surround_time_points < 0),[],2);
    poststim_max = max(frac_move_stimalign(:,surround_time_points > 0),[],2);
    plot(relative_day,(poststim_max-prestim_max)./(poststim_max+prestim_max));
    yline(0);
    ylabel('pre/post move idx');
    xlabel('Day');

    nexttile(t_animal);
    imagesc(surround_time_points,[],frac_move_stimalign); hold on;
    clim([0,1]);
    colormap(gca,AP_colormap('WK'));
    set(gca,'YTick',1:length(recordings),'YTickLabel', ...
        cellfun(@(day,num) sprintf('%d (%s)',num,day(6:end)), ...
        {recordings.day},num2cell(1:length(recordings)),'uni',false));
    xlabel('Time from stim');
    if any(learned_day)
        plot(0,find(learned_day),'.g')
    end

    nexttile(t_animal); hold on
    set(gca,'ColorOrder',copper(length(recordings)));
    plot(surround_time_points,frac_move_stimalign','linewidth',2);
    xline(0,'color','k');
    ylabel('Fraction moving');
    xlabel('Time from stim');
    if any(learned_day)
        AP_errorfill(surround_time_points,frac_move_stimalign(learned_day,:)', ...
            0.02,[0,1,0],0.1,false);
        
        % Store learned day across animals
        learned_day_all(curr_animal_idx) = find(learned_day,1);
    end

    drawnow;

end



















