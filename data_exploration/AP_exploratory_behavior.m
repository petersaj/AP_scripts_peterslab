%% Exploratory behavior analysis


%% Align mousecam to event

% (passive)
stim_window = [0,0.2];
quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
    timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
    timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
    1:length(stimOn_times))';

% use_align = stimOn_times(quiescent_trials);

% stim_x = vertcat(trial_events.values.StimFrequence);
% use_align = stimOn_times(stim_x == 8000  & quiescent_trials);

% stim_x = vertcat(trial_events.values.TrialStimX);
% use_align = stimOn_times(stim_x == 90);

% stim_x = vertcat(trial_events.values.TrialX);
% use_align = stimOn_times(stim_x(1:n_trials) == 90);

% modality = vertcat(trial_events.values(1:n_trials).TaskType);
% use_align = stimOn_times(modality == 1);

% use_align = reward_times;

use_align = stimOn_times;

% use_align = stimOff_times(trial_opacity == 1);

% stim_x = vertcat(trial_events.values.PictureID);
% use_align = stimOn_times(stim_x == 2 & quiescent_trials);

% use_align = stim_center_times(stim_x==90);

% % (task)
% use_align = stimOn_times;
% use_align = stim_move_time;

% Initialize video reader, get average and average difference
vr = VideoReader(mousecam_fn);
cam_im1 = read(vr,1);

surround_times = [-0.2,1];

mousecam_framerate = vr.FrameRate;
surround_frames = round(surround_times*mousecam_framerate);
grab_frames = interp1(mousecam_times,1:length(mousecam_times), ...
    use_align,'previous') + surround_frames;

grab_frames_use = find(~any(isnan(grab_frames),2) & all(grab_frames>0,2));

cam_align_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    diff(surround_frames)+1);
for curr_align = grab_frames_use'
    curr_clip = double(squeeze(read(vr,grab_frames(curr_align,:))));
    cam_align_avg = cam_align_avg + curr_clip./length(grab_frames_use);
    AP_print_progress_fraction(curr_align,length(use_align));
end

surround_t = (surround_frames(1):surround_frames(2))./vr.FrameRate;
AP_imscroll(cam_align_avg,surround_t)
axis image;

surround_t_diff = surround_t(2:end) + diff(surround_t)/2;
AP_imscroll(abs(diff(cam_align_avg,[],3)),surround_t_diff)
axis image;



%% Align mousecam ROI to event

use_cam = mousecam_fn;
use_t = mousecam_times;

% (passive)
% stim_type = vertcat(trial_events.values.TrialStimX);
stim_type = vertcat(trial_events.values.StimFrequence);
use_align = stimOn_times(stim_type == 8000);

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

align_times = stimOn_times;
% align_times = photodiode_times(1:2:end);
% align_times = stimOn_times(align_category_all == 90);
% align_times = stimOn_times(stim_x == 90);
% align_times = stim_move_time;
% align_times = iti_fastmove_times;

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


% [~,sort_idx] = sort([trial_events.values.TrialQuiescence]);
% figure;
% imagesc(surround_time_points,[],event_aligned_wheel_move(sort_idx,:));
% xline(0,'color','r','linewidth',2);
% hold on;
% plot(-[trial_events.values(sort_idx).TrialQuiescence],1:length(trial_events.values),'b','linewidth',2);


%% Behavior across days

% animals = { ...
%     'AM011','AM012','AM014','AM015','AM016','AM017', ...
%     'AM018','AM019','AM021','AM022','AM026','AM029', ...
%     'AP023','AP025'};

animals = {'DS000','DS004','DS014','DS015','DS016'};
animals = {'DS014'};
% Set reaction statistic to use
use_stat = 'mean';
learn_p = 0.05;

% Create master tiled layout
figure('name',sprintf('%s, p < %.2g',use_stat,learn_p));
t = tiledlayout(1,length(animals),'TileSpacing','tight');

% Grab learning day for each mouse
bhv = struct;

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};

    % use_workflow = 'stim_wheel*';
    % use_workflow = 'stim_wheel_right_stage\d';
    use_workflow = 'stim_wheel_right_stage\d_audio_volume';
    % use_workflow = '*audio_volume*';
%     use_workflow = '*audio_frequency*';
%     use_workflow = '*no_change*';
%     use_workflow = '*size*';
%     use_workflow = '*opacity*';
%     use_workflow = '*angle';
%     use_workflow = '*angle_size60';

    recordings = plab.find_recordings(animal,[],use_workflow);

    surround_time = [-5,5];
    surround_sample_rate = 100;
    surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

    n_trials_success = nan(length(recordings),2);
    frac_move_day = nan(length(recordings),1);
    frac_move_stimalign = nan(length(recordings),length(surround_time_points));

    rxn_stat_p = nan(length(recordings),1);
    rxn_stat = nan(length(recordings),1);
    rxn_null_stat = nan(length(recordings),1);

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
        n_trials_success(curr_recording,:) = ...
            [length([trial_events.values.Outcome]), ...
            sum([trial_events.values.Outcome])];

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
        % (skip if only a few trials)
        if n_trials < 10
            continue
        end
       
        [rxn_stat_p(curr_recording), ...
            rxn_stat(curr_recording),rxn_null_stat(curr_recording)] = ...
            AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

        %%%% CAN SUBSTITUTE: 
        % stim_to_move OR stim_to_lastmove
        %%%%

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});
        AP_print_progress_fraction(curr_recording,length(recordings));

    end

    % Define learned day from reaction stat p-value and reaction time
    learned_day = rxn_stat_p < learn_p;

    relative_day = days(datetime({recordings.day}) - datetime({recordings(1).day}))+1;
    nonrecorded_day = setdiff(1:length(recordings),relative_day);

    % Draw in tiled layout nested in master
    t_animal = tiledlayout(t,4,1);
    t_animal.Layout.Tile = curr_animal_idx;
    title(t_animal,animal);

    nexttile(t_animal);
    yyaxis left; plot(relative_day,n_trials_success);
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
    plot(relative_day,rxn_stat)
    set(gca,'YScale','log');
    ylabel(sprintf('Rxn stat: %s',use_stat));
    xlabel('Day');
    if any(nonrecorded_day)
        xline(nonrecorded_day,'--k');
    end
    if any(learned_day)
        xline(relative_day(learned_day),'g');
    end

%     yyaxis right
%     prestim_max = max(frac_move_stimalign(:,surround_time_points < 0),[],2);
%     poststim_max = max(frac_move_stimalign(:,surround_time_points > 0),[],2);
%     stim_move_frac_ratio = (poststim_max-prestim_max)./(poststim_max+prestim_max);
%     plot(relative_day,stim_move_frac_ratio);
%     yline(0);
%     ylabel('pre/post move idx');
%     xlabel('Day');

    yyaxis right
    plot(relative_day,(rxn_stat-rxn_null_stat)./(rxn_stat+rxn_null_stat));
    yline(0);
    ylabel(sprintf('Rxn stat idx: %s',use_stat));
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
    end

    drawnow;

    % Store behavior across animals
    bhv(curr_animal_idx).rxn_stat = rxn_stat;
    bhv(curr_animal_idx).rxn_stat = rxn_null_stat;
    bhv(curr_animal_idx).learned_day = learned_day;

end



















