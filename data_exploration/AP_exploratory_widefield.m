%% Exploratory widefield analysis

%% Experiment scroller

% ap.expscroll(wf_U_raw{1},wf_V_raw{1},wf_t_all{1})
% ap.expscroll(wf_U_raw{2},wf_V_raw{2},wf_t_all{2})
% ap.expscroll(wf_U,wf_Vdf,wf_times)
% ap.expscroll(wf_U,wf_V,wf_times)

ap.expscroll(wf_U,wf_V,wf_times,mousecam_fn,mousecam_times)


%% Align widefield to event

if contains(bonsai_workflow,'lcr')

    % LCR passive: align to quiescent stim onset
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        1:length(stimOn_times))';

    align_times = stimOn_times(quiescent_trials);
    align_category_all = vertcat(trial_events.values.TrialStimX);
    align_category = align_category_all(quiescent_trials);

    baseline_times = stimOn_times(quiescent_trials);

elseif contains(bonsai_workflow,'stim_wheel')
    
    % Task: align to stim/move/reward
    rewarded_trials = logical([trial_events.values.Outcome]');

    use_trials = rewarded_trials(1:n_trials);
    align_times = reshape([ ...
        stimOn_times(use_trials), ...
        stim_move_time(use_trials), ...
        reward_times(use_trials)],[],1);
    align_category = reshape(ones(sum(use_trials),3).*[1,2,3],[],1);
    baseline_times = repmat(stimOn_times(use_trials),3,1);

%     use_trials = rewarded_trials(1:n_trials);
%     align_times = reshape([ ...
%         stimOn_times(use_trials), ...
%         stim_move_time(use_trials), ...
%         reward_times(1:sum(rewarded_trials))],[],1);
%     align_category = reshape(ones(sum(use_trials),3).*[1,2,3],[],1);
%     baseline_times = repmat(stimOn_times(use_trials),3,1);

elseif contains(bonsai_workflow,'sparse_noise')

    % (sparse noise)
    px_x = 23;
    px_y = 4;
    align_times = ...
        stim_times(find( ...
        (noise_locations(px_y,px_x,1:end-1) == 128 & ...
        noise_locations(px_y,px_x,2:end) == 255) | ...
        (noise_locations(px_y,px_x,1:end-1) == 128 & ...
        noise_locations(px_y,px_x,2:end) == 0))+1);
    align_category = ones(size(align_times));

end

surround_window = [-1,2];
baseline_window = [-0.5,-0.1];

surround_samplerate = 60;

t = surround_window(1):1/surround_samplerate:surround_window(2);
baseline_t = baseline_window(1):1/surround_samplerate:baseline_window(2);

peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);
baseline_event_t = reshape(baseline_times,[],1) + reshape(baseline_t,1,[]);

use_U = wf_U;
use_V = wf_V;
use_wf_t = wf_times;

aligned_v = reshape(interp1(use_wf_t,use_V',peri_event_t,'previous'), ...
    length(align_times),length(t),[]);
aligned_baseline_v = nanmean(reshape(interp1(use_wf_t,use_V',baseline_event_t,'previous'), ...
    length(baseline_times),length(baseline_t),[]),2);

aligned_v_baselinesub = aligned_v - aligned_baseline_v;

align_id = findgroups(reshape(align_category,[],1));

aligned_v_avg = permute(splitapply(@nanmean,aligned_v_baselinesub,align_id),[3,2,1]);
aligned_px_avg = plab.wf.svd2px(use_U,aligned_v_avg);

AP_imscroll(aligned_px_avg,t);
colormap(AP_colormap('PWG'));
clim(prctile(abs(aligned_px_avg(:)),100).*[-1,1]);
axis image;


%% Sparse noise retinotopy (single day)

ap.widefield_retinotopy


%% ~~~~~~~~ ALIGN WIDEFIELD

%% Create day alignment

animal = 'AP010';

ap.align_widefield([],animal,[],'new_days');

%% Check day alignment

animal = 'AP010';

recordings = plab.find_recordings(animal);
wf_days_idx = cellfun(@(x) any(x),{recordings.widefield});
wf_recordings = recordings(wf_days_idx);

avg_im_aligned = cell(size(wf_recordings));
for curr_day = 1:length(wf_recordings)
    day = wf_recordings(curr_day).day;

    img_path = plab.locations.make_server_filename( ...
        animal,day,[],'widefield');

    avg_im_h = readNPY([img_path filesep 'meanImage_violet.npy']);
    avg_im_n = readNPY([img_path filesep 'meanImage_blue.npy']);
    avg_im_aligned{curr_day} = [ap.align_widefield(avg_im_n,animal,day), ...
        ap.align_widefield(avg_im_h,animal,day)];
end

% Plot average
c = prctile(reshape([avg_im_aligned{:}],[],1),[0,99.9]);
AP_imscroll(cat(3,avg_im_aligned{:}),{recordings.day});
caxis(c);
axis image;
set(gcf,'Name',animal);

%% Create animal alignment (sparse noise retinotopy: batch, save, align)
% NOTE: need to do day alignment first

animal = 'AP010';
overwrite_retinotopy = false;

% Check if aligned mean retinotopy is saved, create if not
alignment_path = fullfile(plab.locations.server_path, ...
    'Users','Andy_Peters','widefield_alignment');
retinotopy_path = fullfile(alignment_path,'retinotopy');
retinotopy_fn = fullfile(retinotopy_path,sprintf('%s_retinotopy.mat',animal));

if overwrite_retinotopy || ~exist(retinotopy_fn,'file')
    workflow = 'sparse_noise';
    recordings = plab.find_recordings(animal,[],workflow);

    vfs_all = cell(length(recordings),1);
    disp('Creating retinotopy...');
    for curr_day = 1:length(recordings)
        rec_day = recordings(curr_day).day;
        rec_time = recordings(curr_day).recording{end};
        verbose = true;
        ap.load_recording;

        ap.widefield_retinotopy;
        vfs_all{curr_day} = vfs;
    end

    % Save retinotopy from all days
    retinotopy = struct('animal',animal,'day',reshape({recordings.day},[],1),'vfs',vfs);
    save(retinotopy_fn,'retinotopy');
    fprintf('Saved %s\n',retinotopy_fn);
end

% Load pre-saved retinotopy
load(retinotopy_fn);

% Align VFS across days
aligned_vfs = cell(length(retinotopy),1);
for curr_day = 1:length(retinotopy)
        aligned_vfs{curr_day} = ap.align_widefield(retinotopy(curr_day).vfs, ...
            retinotopy.animal,retinotopy(curr_day).day,'day_only');
end

% Plot day-aligned VFS
AP_imscroll(cat(3,aligned_vfs{:}));
colormap(AP_colormap('BWR'));clim([-1,1]);
axis image off;
title('Day-aligned VFS')

% Align across-day mean VFS to master animal
vfs_mean = nanmean(cat(3,aligned_vfs{:}),3);
ap.align_widefield(vfs_mean,animal,[],'new_animal');

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

%% ~~~~~~~~ BATCH


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




