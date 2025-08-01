%% Exploratory widefield analysis


%% Align widefield to event

if contains(bonsai_workflow,'passive')

    % LCR passive: align to quiescent stim onset
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        1:length(stimOn_times))';

    align_times = stimOn_times(quiescent_trials);
    if isfield(trial_events.values,'TrialStimX')
        align_category_all = vertcat(trial_events.values.TrialStimX);
    elseif isfield(trial_events.values,'StimFrequence')
        align_category_all = vertcat(trial_events.values.StimFrequence);
    end
    align_category = align_category_all(quiescent_trials);

    baseline_times = stimOn_times(quiescent_trials);

elseif contains(bonsai_workflow,'stim_wheel')

    % Task: align to stim/move/reward
    rewarded_trials = logical([trial_events.values.Outcome]');

    if isfield(trial_events.values,'TaskType')
        use_trials = rewarded_trials(1:n_trials) & ...
            vertcat(trial_events.values(1:n_trials).TaskType) == 0;
    else
        use_trials = true(size(stim_to_move));
    end

    align_times = [ ...
        stimOn_times(use_trials); ...
        stim_move_time(use_trials); ...
        reward_times];
    align_category = vertcat( ...
        1*ones(sum(use_trials),1), ...
        2*ones(sum(use_trials),1), ...
        3*ones(length(reward_times),1));
    baseline_times = vertcat(...
        stimOn_times(use_trials), ...
        stimOn_times(use_trials), ...
        reward_times);

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

    baseline_times = align_times;

elseif contains(bonsai_workflow,'visual')

    % Visual conditioning: stim and reward
    if isfield(trial_events.values,'TrialX')
        stim_x = vertcat(trial_events.values.TrialX);
    else
        stim_x = repmat(90,length(stimOn_times),1);
    end
    align_times = stimOn_times(1:n_trials);
    align_category = stim_x(1:n_trials);
    baseline_times = align_times;

%     % (new: just one stim)
%     align_times = stimOn_times;
%     align_category = ones(size(align_times));
%     baseline_times = stimOn_times;

end

surround_window = [-0.5,1];
baseline_window = [-0.5,-0.1];

surround_samplerate = 35;

t = surround_window(1):1/surround_samplerate:surround_window(2);
baseline_t = baseline_window(1):1/surround_samplerate:baseline_window(2);

peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);
baseline_event_t = reshape(baseline_times,[],1) + reshape(baseline_t,1,[]);

use_U = wf_U;
use_V = wf_V;
use_wf_t = wf_t;

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
set(gcf,'name',sprintf('%s %s %s',animal,rec_day,bonsai_workflow));


%% Align ROI to event

[roi_trace,roi_mask] = ap.wf_roi(wf_U,wf_V,wf_avg);

% % LCR passive: align to quiescent stim onset
% stim_window = [0,0.5];
% quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
%     timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
%     timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
%     1:length(stimOn_times))';
% 
% align_times = stimOn_times(quiescent_trials);
% if isfield(trial_events.values,'TrialStimX')
%     align_category_all = vertcat(trial_events.values.TrialStimX);
% elseif isfield(trial_events.values,'StimFrequence')
%     align_category_all = vertcat(trial_events.values.StimFrequence);
% end
% align_category = align_category_all(quiescent_trials);
% baseline_times = stimOn_times(quiescent_trials);

% (task)
align_times = stimOn_times(1:n_trials);
baseline_times = stimOn_times(1:n_trials);
align_category = ones(size(align_times));

% Align ROI trace to align times
surround_window = [-1,2];
baseline_window = [-0.5,-0.1];

surround_samplerate = 35;

t = surround_window(1):1/surround_samplerate:surround_window(2);
peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

aligned_trace = reshape(interp1(wf_t,roi_trace',peri_event_t,'previous'), ...
    length(align_times),length(t),[]);

aligned_trace_baselinesub = aligned_trace - ...
    mean(aligned_trace(:,t >= baseline_window(1) & t <= baseline_window(2)),2);

align_id = findgroups(reshape(align_category,[],1));

figure;
h = tiledlayout(1,max(align_id));
for curr_id = 1:max(align_id)

    curr_trials = find(align_id == curr_id);
    [~,sort_idx] = sort(stim_to_move(curr_trials));

    nexttile;
    imagesc(t,[],aligned_trace_baselinesub(curr_trials(sort_idx),:))
    clim(0.8*max(max(aligned_trace_baselinesub,[],'all')).*[-1,1]);
    xline(0);
end
colormap(AP_colormap('PWG'));



%% Passive kernel

switch bonsai_workflow
    case 'lcr_passive'
        stim_type = vertcat(trial_events.values.TrialStimX);
    case 'hml_passive_audio'
        stim_type = vertcat(trial_events.values.StimFrequence);
end

time_bins = [wf_t;wf_t(end)+1/wf_framerate];

stim_regressors = cell2mat(arrayfun(@(x) ...
    histcounts(stimOn_times(stim_type == x),time_bins), ...
    unique(stim_type),'uni',false));

n_components = 500;

frame_shifts = -5:20;
lambda = 10;
cv = 1;

skip_t = 10; % seconds start/end to skip for artifacts
skip_frames = round(skip_t*wf_framerate);
[kernels,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
    stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv);

kernels_px = plab.wf.svd2px(wf_U(:,:,1:size(kernels,1)),kernels);
ap.imscroll(kernels_px,frame_shifts);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image

% trying SVM/SVR
% mdl = fitcsvm(wf_V',stim_regressors(1,:)');
% cvmodel = crossval(mdl);
% classLoss = kfoldLoss(cvmodel);
% 
% r = mdl.predict(wf_V');

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

nexttile; imagesc(kernels_px_max(:,:,3)); 
clim(col_lim); axis image off;
colormap(gca,cs_minus_color);
ap.wf_draw('ccf','k');



%% Sparse noise retinotopy (single day)

% Load data
animal = 'AP026';
workflow = 'sparse_noise';

% % Specific day
% rec_day = '2023-12-06';
% rec_time = plab.find_recordings(animal,rec_day,workflow).recording{end};

% Relative day
recordings = plab.find_recordings(animal,[],workflow);
use_day = 1;
% use_day = length(recordings);
rec_day = recordings(use_day).day;
rec_time = recordings(use_day).recording{end};

load_parts.widefield = true;

verbose = true;
ap.load_recording;

% Get retinotopy
ap.wf_retinotopy


%% ~~~~~~~~ ALIGN WIDEFIELD

%% Create alignments

animal = 'AP005';

% Get and save VFS maps for animal
plab.wf.retinotopy_vfs_batch(animal);

% Create across-day alignments
plab.wf.wf_align([],animal,[],'new_days');

% Create across-animal alignments
plab.wf.wf_align([],animal,[],'new_animal');


%% View aligned days

animal = 'HA009';

recordings = plab.find_recordings(animal);
wf_days_idx = cellfun(@(x) any(x),{recordings.widefield});
wf_recordings = recordings(wf_days_idx);

avg_im_aligned = cell(size(wf_recordings));
for curr_day = 1:length(wf_recordings)
    day = wf_recordings(curr_day).day;

    img_path = plab.locations.filename('server', ...
        animal,day,[],'widefield');

    avg_im_n = readNPY([img_path filesep 'meanImage_blue.npy']);
    avg_im_h = readNPY([img_path filesep 'meanImage_violet.npy']);

%         % (to concatenate)
%         avg_im_aligned{curr_day} = [plab.wf.wf_align(avg_im_n,animal,day), ...
%             plab.wf.wf_align(avg_im_h,animal,day)];

    % (blue only)
    avg_im_aligned{curr_day} = plab.wf.wf_align(avg_im_n,animal,day);
end

% Plot average
c = prctile(reshape([avg_im_aligned{:}],[],1),[0,99.9]);
AP_imscroll(cat(3,avg_im_aligned{:}),{wf_recordings.day});
caxis(c);
axis image;
set(gcf,'Name',animal);

%% Regression task > widefield

% Parameters for regression
regression_params.use_svs = 1:100;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Get time points to bin
time_bin_centers = wf_t;
time_bins = [wf_t;wf_t(end)+1/wf_framerate];

% Regressors
stim_regressors = histcounts(stimOn_times,time_bins);
reward_regressors = histcounts(reward_times,time_bins);

stim_move_regressors = histcounts(stim_move_time,time_bins);
nonstim_move_times = ...
    setdiff(timelite.timestamps(find(diff(wheel_move) == 1)+1), ...
    stim_move_regressors);
nonstim_move_regressors = histcounts(nonstim_move_times,time_bins);

% Concatenate selected regressors, set parameters
task_regressors = {stim_regressors;reward_regressors;stim_move_regressors;nonstim_move_regressors};
task_regressor_labels = {'Stim','Reward','Stim move','Nonstim move'};

task_t_shifts = { ...
    [-0.2,2]; ... % stim
    [-0.2,2];  ... % outcome
    [-0.2,2];  ... % nonstim move
    [-0.2,2]};    % stim move

% task_regressors = {stim_regressors;reward_regressors};
% task_regressor_labels = {'Stim','Reward'};
% 
% task_t_shifts = { ...
%     [-0.2,2]; ... % stim
%     [-0.2,2]};    % reward

task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(wf_framerate)): ...
    round(x(2)*(wf_framerate)),task_t_shifts,'uni',false);
lambda = 0;
zs = [false,false];
cvfold = 5;
use_constant = true;
return_constant = false;
use_svs = 100;

[fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
    AP_regresskernel(task_regressors,wf_V(1:use_svs,:),task_regressor_sample_shifts, ...
    lambda,zs,cvfold,return_constant,use_constant);
% 
% fluor_taskpred = ...
%     interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
% 
% fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
%     interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
%     t_peri_event),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));



a = cellfun(@(x) plab.wf.svd2px(wf_U(:,:,1:use_svs), ...
    permute(x,[3,2,1])),fluor_taskpred_k,'uni',false);

AP_imscroll(cat(4,a{:}));
axis image;
clim(max(abs(clim)).*[-1,1]); 
colormap(AP_colormap('PWG'));


%% ~~~~~~~~ BATCH


%% TESTING BATCH PASSIVE WIDEFIELD

animal = 'HA002';
passive_workflow = 'lcr_passive';
% passive_workflow = 'hml_passive_audio';
recordings_passive = plab.find_recordings(animal,[],passive_workflow);

% training_workflow = 'stim_wheel*';
training_workflow = 'visual*';
% training_workflow = '*audio_volume*';
recordings_training = plab.find_recordings(animal,[],training_workflow);

% (use recordings on training days)
recordings = recordings_passive( ...
    cellfun(@any,{recordings_passive.widefield}) & ...
    ismember({recordings_passive.day},{recordings_training.day}));

% % (use recordings on or before last training day)
% recordings = recordings_passive( ...
%     cellfun(@any,{recordings_passive.widefield}) & ...
%     ~[recordings_passive.ephys] & ...
%     days(datetime({recordings_passive.day}) - datetime(recordings_training(end).day)) <= 0);

% (use all passive recordings)
% recordings = recordings_passive( ...
%     cellfun(@any,{recordings_passive.widefield}) & ...
%     ~[recordings_passive.ephys]);

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

    try
    load_parts.widefield = true;
    ap.load_recording;
    catch me
        warning('%s %s %s: load error, skipping \n >> %s', ...
            animal,rec_day,rec_time,me.message)
        continue
    end

    % Get quiescent trials and stim onsets/ids
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        1:length(stimOn_times))';

    align_times = stimOn_times(quiescent_trials);
    align_category_all = vertcat(trial_events.values.TrialStimX);
%     align_category_all = vertcat(trial_events.values.StimFrequence);
    align_category = align_category_all(quiescent_trials);

    % Align to stim onset
    surround_window = [-0.5,1];
    surround_samplerate = 35;
    t = surround_window(1):1/surround_samplerate:surround_window(2);
    peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

    aligned_v = reshape(interp1(wf_t,wf_V',peri_event_t,'previous'), ...
        length(align_times),length(t),[]);

    align_id = findgroups(align_category);
    aligned_v_avg = permute(splitapply(@nanmean,aligned_v,align_id),[3,2,1]);
    aligned_v_avg_baselined = aligned_v_avg - nanmean(aligned_v_avg(:,t < 0,:),2);

    % Convert to pixels and package
    aligned_px_avg = plab.wf.svd2px(wf_U,aligned_v_avg_baselined);
    wf_px{curr_recording} = aligned_px_avg;

    % Prep for next loop
    AP_print_progress_fraction(curr_recording,length(recordings));
    clearvars('-except',preload_vars{:});

end

surround_window = [-0.5,1];
surround_samplerate = 35;
t = surround_window(1):1/surround_samplerate:surround_window(2);
a = cellfun(@(x) max(x(:,:,t > 0 & t < 0.2,3),[],3), ...
    wf_px(cellfun(@(x) ~isempty(x),wf_px)),'uni',false);
c = (max(cellfun(@(x) max(x(:)),a)).*[-1,1])/2;

figure('Name',animal');
tiledlayout('flow','TileSpacing','none')
for i = 1:length(a)
    nexttile;imagesc(a{i}); axis image off;
    clim(c); colormap(AP_colormap('PWG'));
end
ap.imscroll(cat(3,a{:}));
axis image;
clim(max(abs(clim)).*[-1,1]); colormap(AP_colormap('PWG'));

a = cat(5,wf_px{:});
b = squeeze(a(:,:,:,3,:));
ap.imscroll(b);
axis image;
clim(max(abs(clim)).*[-1,1]); colormap(AP_colormap('PWG'));

% % (reflect widefield - taken out for now)
% a = cellfun(@(x) x-ap.wf_reflect(x),wf_px,'uni',false);
% b = cellfun(@(x) mean(x(:,:,t > 0.05 & t < 0.15,3),3),a,'uni',false);
% c = (max(cellfun(@(x) max(x(:)),a)).*[-1,1])/2;
% figure('Name',animal');
% tiledlayout('flow')
% for i = 1:length(a)
%     nexttile; imagesc(b{i}); axis image off;
% end
% AP_imscroll(cat(3,b{:}));
% axis image;
% clim(c); colormap(AP_colormap('PWG'));

a = cat(5,wf_px{:});
b = squeeze(a(:,:,:,3,:));
b2 = squeeze(mean(b(:,:,t > 0 & t < 0.2,:),3));









