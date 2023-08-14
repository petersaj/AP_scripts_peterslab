%% Exploratory widefield analysis

%% Experiment scroller

% ap.expscroll(wf_U_raw{1},wf_V_raw{1},wf_t_all{1})
% ap.expscroll(wf_U_raw{2},wf_V_raw{2},wf_t_all{2})
% ap.expscroll(wf_U,wf_Vdf,wf_times)
% ap.expscroll(wf_U,wf_V,wf_times)

ap.expscroll(wf_U,wf_V,wf_times,mousecam_fn,mousecam_times)


%% Align widefield to event

% % (passive)
% align_times_all = stimOn_times;
% align_category_all = vertcat(trial_events.values.TrialStimX);
% % (get only quiescent trials)
% stim_window = [0,0.5];
% quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
%     timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
%     timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
%     1:length(stimOn_times))';
% align_times = align_times_all(quiescent_trials);
% align_category = align_category_all(quiescent_trials);

% (task stim)
align_times = stimOn_times;
align_category = ones(size(align_times));

% % (task move)
% align_times = stim_move_time;
% align_category = ones(size(align_times));

% % (task rewards)
% align_times = reward_times;
% align_category = ones(size(align_times));

% % (sparse noise)
% px_x = 23;
% px_y = 4;
% align_times = ...
%     stim_times(find( ...
%     (noise_locations(px_y,px_x,1:end-1) == 128 & ...
%     noise_locations(px_y,px_x,2:end) == 255) | ...
%     (noise_locations(px_y,px_x,1:end-1) == 128 & ...
%     noise_locations(px_y,px_x,2:end) == 0))+1);
% align_category = ones(size(align_times));

surround_window = [-1,2];

surround_samplerate = 35;
t = surround_window(1):1/surround_samplerate:surround_window(2);
peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

use_U = wf_U;
use_V = wf_V;
use_wf_t = wf_times;

aligned_v = reshape(interp1(use_wf_t,use_V',peri_event_t,'previous'), ...
    length(align_times),length(t),[]);

align_id = findgroups(align_category);
aligned_v_avg = permute(splitapply(@nanmean,aligned_v,align_id),[3,2,1]);
aligned_v_avg_baselined = aligned_v_avg - nanmean(aligned_v_avg(:,t < 0,:),2);

aligned_px_avg = plab.wf.svd2px(use_U,aligned_v_avg_baselined);

AP_imscroll(aligned_px_avg,t);
colormap(AP_colormap('PWG'));
clim(prctile(abs(aligned_px_avg(:)),99.5).*[-1,1]);
axis image;

%% Sparse noise retinotopy (single day)

ap.wf_retinotopy

%% Create animal alignment (sparse noise retinotopy: batch, save, align)

animal = 'AP010';

workflow = 'sparse_noise';
recordings = plab.find_recordings(animal,[],workflow);

vfs_all = cell(length(recordings),1);        
for curr_day = 1:length(recordings)
    rec_day = recordings(curr_day).day;
    rec_time = recordings(curr_day).recording{end};
    verbose = true;
    ap.load_recording;

    ap.wf_retinotopy;
    vfs_all{curr_day} = vfs;
end

% Save retinotopy from all days
retinotopy = struct('animal',animal,'day',reshape({recordings.day},[],1),'vfs',vfs);
alignment_path = fullfile(plab.locations.server_path, ...
    'Users','Andy_Peters','widefield_alignment');
retinotopy_path = fullfile(alignment_path,'retinotopy');
retinotopy_fn = fullfile(retinotopy_path,sprintf('%s_retinotopy.mat',animal));
save(retinotopy_fn,'retinotopy');
disp(['Saved ' retinotopy_fn]);

% Align VFS across days
aligned_vfs = cell(length(retinotopy),1);
for curr_day = 1:length(retinotopy)
        aligned_vfs{curr_day} = ap.align_wf(retinotopy(curr_day).vfs, ...
            retinotopy.animal,retinotopy(curr_day).day,'day_only');
end

% Plot day-aligned VFS
AP_imscroll(cat(3,aligned_vfs{:}));
colormap(AP_colormap('BWR'));clim([-1,1]);
axis image off;
title('Day-aligned VFS')

% Align across-day mean VFS to master animal
% (TEMPORARY: flip VFS sign, flipped from previous)
vfs_mean = nanmean(cat(3,aligned_vfs{:}),3);
ap.align_widefield(-vfs_mean,animal,[],'new_animal');

% Plot master-aligned average VFS with CCF overlay
master_aligned_vfs = cell(length(retinotopy),1);
for curr_day = 1:length(retinotopy)
        master_aligned_vfs{curr_day} = ...
            ap.align_widefield(retinotopy(curr_day).vfs, ...
            animal,retinotopy(curr_day).day);
end
figure;imagesc(nanmean(cat(3,master_aligned_vfs{:}),3));
colormap(AP_colormap('BWR'));clim([-1,1]);
axis image off;
ap.draw_wf_ccf('ccf_aligned',[0.5,0.5,0.5]);
title('Master-aligned average VFS');


%% TESTING BATCH PASSIVE WIDEFIELD

animal = 'AP009';
use_workflow = 'lcr_passive';
recordings = ap.find_recordings(animal,[],use_workflow);

% temp: cut out bad days
% recordings(ismember({recordings.day},{'2023-05-04','2023-05-05'})) = [];
% recordings(ismember({recordings.day},{'2023-05-08'})) = [];

wf_px = cell(size(recordings));

for curr_recording = 1:length(recordings)
    
    % Grab pre-load vars
    preload_vars = who;

    % Load data
    rec_day = recordings(curr_recording).day;
    rec_time = recordings(curr_recording).recording{end};
    if ~recordings(curr_recording).widefield(end)
        continue
    end
    ap.load_recording;

    % Align to stim and store
    % (passive)
    align_times_all = photodiode_times(1:2:end);
    align_category_all = vertcat(trial_events.values.TrialStimX);
    % (get only quiescent trials)
    framerate = 30;
    wheel_window = [0,0.5];
    wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
    wheel_window_t_peri_event = align_times_all + wheel_window_t;
    event_aligned_move = interp1(timelite.timestamps, ...
        +wheel_move,wheel_window_t_peri_event,'previous');
    quiescent_trials = ~any(event_aligned_move,2);
    align_times = align_times_all(quiescent_trials);
    align_category = align_category_all(quiescent_trials);

    surround_window = [-0.5,1];
    surround_samplerate = 35;
    t = surround_window(1):1/surround_samplerate:surround_window(2);
    peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

    aligned_v = reshape(interp1(wf_times,wf_V',peri_event_t,'previous'), ...
        length(align_times),length(t),[]);

    align_id = findgroups(align_category);
    aligned_v_avg = permute(splitapply(@nanmean,aligned_v,align_id),[3,2,1]);
    aligned_v_avg_baselined = aligned_v_avg - nanmean(aligned_v_avg(:,t < 0,:),2);

    aligned_px_avg = plab.wf.svd2px(wf_U,aligned_v_avg_baselined);

    wf_px{curr_recording} = aligned_px_avg;

    AP_print_progress_fraction(curr_recording,length(recordings));

    % Clear vars except pre-load for next loop
    clearvars('-except',preload_vars{:});

end

surround_window = [-0.5,1];
surround_samplerate = 35;
t = surround_window(1):1/surround_samplerate:surround_window(2);
a = cellfun(@(x) max(x(:,:,t > 0 & t < 0.2,3),[],3), ...
    wf_px(cellfun(@(x) ~isempty(x),wf_px)),'uni',false);
c = [0,max(cellfun(@(x) max(x(:)),a))];

figure('Name',animal');
tiledlayout('flow')
for i = 1:length(a)
nexttile;imagesc(a{i});clim(c);axis image off;
end




