%% Exploratory widefield analysis


%% Align widefield to event

if contains(bonsai_workflow,'lcr')

    % LCR passive: align to quiescent stim onset
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        1:length(stimOn_times))';

    quiescent_trials(70:end) = false;

    align_times = stimOn_times(quiescent_trials);
    align_category_all = vertcat(trial_events.values.TrialStimX);
    align_category = align_category_all(quiescent_trials);

    baseline_times = stimOn_times(quiescent_trials);

elseif contains(bonsai_workflow,'stim_wheel')

    % Task: align to stim/move/reward
    rewarded_trials = logical([trial_events.values.Outcome]');

    use_trials = rewarded_trials(1:n_trials);
    %     align_times = [ ...
    %         stimOn_times(use_trials); ...
    %         stim_move_time(use_trials); ...
    %         reward_times_task(use_trials)];
    align_times = [ ...
        stimOn_times(use_trials); ...
        stim_move_time(use_trials); ...
        reward_times];
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

elseif contains(bonsai_workflow,'visual_conditioning')

    % Only use completed (rewarded) trials
    n_trials = length(reward_times);

    % Visual conditioning: stim and reward
    align_times = [ ...
        stimOn_times(1:n_trials); ...
        reward_times(1:n_trials)];
    align_category = reshape(ones(n_trials,2).*[1,2],[],1);
    baseline_times = repmat(stimOn_times(1:n_trials),2,1);

end

surround_window = [-1,4];
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

%% Sparse noise retinotopy (single day)

% Load data
animal = 'AP014';
workflow = 'sparse_noise';

% % Specific day
% rec_day = '2023-12-06';
% rec_time = plab.find_recordings(animal,rec_day,workflow).recording{end};

% Relative day
recordings = plab.find_recordings(animal,[],workflow);
% use_day = 1;
use_day = length(recordings);
rec_day = recordings(use_day).day;
rec_time = recordings(use_day).recording{end};

load_parts.widefield = true;

verbose = true;
ap.load_recording;

% Get retinotopy
ap.wf_retinotopy


%% ~~~~~~~~ ALIGN WIDEFIELD

%% Create day alignment

animal = 'AP016';

ap.wf_align([],animal,[],'new_days');

%% Batch sparse noise retinotopy

% animals = {'AP014','AP015'};

animals = {'AP016'};

for curr_animal = 1:length(animals)

    preload_vars = who;

    animal = animals{curr_animal};

    workflow = 'sparse_noise';
    recordings = plab.find_recordings(animal,[],workflow);

    vfs_all = cell(length(recordings),1);
    disp('Creating retinotopy...');
    for curr_day = 1:length(recordings)
        rec_day = recordings(curr_day).day;
        rec_time = recordings(curr_day).recording{end};

        load_parts.widefield = true;
        load_parts.widefield_align = false;
        verbose = true;
        try
            ap.load_recording;
            ap.wf_retinotopy;
            vfs_all{curr_day} = vfs;
        catch me
            % If there's an error, skip to next day
            continue
        end
    end

    % Save retinotopy from all days which have VFS
    retinotopy_fn = fullfile(plab.locations.server_path, ...
        'Users','Andy_Peters','widefield_alignment','retinotopy', ...
        sprintf('retinotopy_%s.mat',animal));

    use_recordings = cellfun(@(x) ~isempty(x),vfs_all);

    retinotopy = struct;
    retinotopy.animal = animal;
    retinotopy.day = {recordings(use_recordings).day};
    retinotopy.vfs = vfs_all;

    save(retinotopy_fn,'retinotopy');
    fprintf('Saved %s\n',retinotopy_fn);

    % Clear variables for next loop
    clearvars -except preload_vars

end


%% Create animal alignment (sparse noise retinotopy: batch, save, align)
% NOTE: need to do day alignment first

animal = 'AP016';

% Load pre-saved retinotopy
retinotopy_fn = fullfile(plab.locations.server_path, ...
    'Users','Andy_Peters','widefield_alignment','retinotopy', ...
    sprintf('retinotopy_%s.mat',animal));

if exist(retinotopy_fn,'file')
    load(retinotopy_fn);
else
    error('No retinotopy saved: %s',animal);
end

% Align VFS across days
aligned_vfs = cell(length(retinotopy),1);
for curr_day = 1:length(retinotopy)
    aligned_vfs{curr_day} = ap.wf_align(retinotopy(curr_day).vfs, ...
        retinotopy(1).animal,retinotopy(curr_day).day,'day_only');
end

% Plot day-aligned VFS
AP_imscroll(cat(3,aligned_vfs{:}));
colormap(AP_colormap('BWR'));clim([-1,1]);
axis image off;
title('Day-aligned VFS')

% Align across-day mean VFS to master animal
vfs_mean = nanmean(cat(3,aligned_vfs{:}),3);
ap.wf_align(vfs_mean,animal,[],'new_animal');

% Plot master-aligned average VFS with CCF overlay
master_aligned_vfs = cell(length(retinotopy),1);
for curr_day = 1:length(retinotopy)
    master_aligned_vfs{curr_day} = ...
        ap.wf_align(retinotopy(curr_day).vfs, ...
        animal,retinotopy(curr_day).day);
end
figure;imagesc(nanmean(cat(3,master_aligned_vfs{:}),3));
colormap(AP_colormap('BWR'));clim([-1,1]);
axis image off;
ap.wf_draw('ccf',[0.5,0.5,0.5]);
title('Master-aligned average VFS');

%% View aligned days

animal = 'AP016';

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

    %     % (to concatenate)
    %     avg_im_aligned{curr_day} = [ap.wf_align(avg_im_n,animal,day), ...
    %         ap.wf_align(avg_im_h,animal,day)];

    % (blue only)
    avg_im_aligned{curr_day} = ap.wf_align(avg_im_n,animal,day);
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

% stim_move_regressors = histcounts(stim_move_time,time_bins);
% nonstim_move_times = ...
%     setdiff(timelite.timestamps(find(diff(wheel_move) == 1)+1), ...
%     stim_move_regressors);
% nonstim_move_regressors = histcounts(nonstim_move_times,time_bins);


% Concatenate selected regressors, set parameters
% task_regressors = {stim_regressors;reward_regressors;stim_move_regressors;nonstim_move_regressors};
% task_regressor_labels = {'Stim','Reward','Stim move','Nonstim move'};
% 
% task_t_shifts = { ...
%     [-0.2,2]; ... % stim
%     [-0.2,2];  ... % outcome
%     [-0.2,2];  ... % nonstim move
%     [-0.2,2]};    % stim move

task_regressors = {stim_regressors;reward_regressors};
task_regressor_labels = {'Stim','Reward'};

task_t_shifts = { ...
    [-0.2,2]; ... % stim
    [-0.2,2]};    % reward

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

animal = 'AP010';
use_workflow = 'lcr_passive';
recordings = plab.find_recordings(animal,[],use_workflow);

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

    load_parts.widefield = true;
    ap.load_recording;

    % Get quiescent trials and stim onsets/ids
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        1:length(stimOn_times))';

    align_times = stimOn_times(quiescent_trials);
    align_category_all = vertcat(trial_events.values.TrialStimX);
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
AP_imscroll(cat(3,a{:}));
axis image;
clim(max(abs(clim)).*[-1,1]); colormap(AP_colormap('PWG'));

a = cat(5,wf_px{:});
b = squeeze(a(:,:,:,3,:));
AP_imscroll(b);
axis image;
clim(max(abs(clim)).*[-1,1]); colormap(AP_colormap('PWG'));

% (reflect widefield - taken out for now)
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












