%% Da's secondary-ish project with multiple visual stimuli

%% TESTING: all animals that have learned the stim task

% Find all animals from task
% task_dir = dir(fullfile(plab.locations.server_data_path, ...
%     '**','bonsai','**','stim_wheel_right_stage2.bonsai'));

task_dir_topfolders = unique(extractBetween({task_dir.folder}, ...
    [plab.locations.server_data_path,filesep], ...
    filesep));

animals = task_dir_topfolders(matches(task_dir_topfolders,lettersPattern(2)+digitsPattern(3)));

% Get learned days
learned_days = struct('animal',animals,'days',cell(size(animals)),'learned',cell(size(animals)));
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    task_recordings = plab.find_recordings(animal,[],{'stim_wheel_right_stage1','stim_wheel_right_stage2'});

    for curr_recording = 1:length(task_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = task_recordings(curr_recording).day;
        rec_time = task_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % (only use if there were at least 10 trials)
        if n_trials < 10
            continue
        end

        % Get reaction stat
        use_stat = 'mean';
        rxn_stat_p = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);
        learned_day = rxn_stat_p < 0.05;

        % Package data
        learned_days(curr_animal).days{curr_recording} = rec_day;
        learned_days(curr_animal).learned(curr_recording) = learned_day;

        ap.print_progress_fraction(curr_recording,length(task_recordings));
    end
end

% Save
save_path = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data';
save_fn = fullfile(save_path,'position_learned_days');
save(save_fn,'learned_days');


%% Task-specific: AVERAGE passive stim

% task_workflow = 'stim_wheel_right_stage\d_size_up';
% task_workflow = 'stim_wheel_right_stage\d_angle_size60';
task_workflow = 'stim_wheel_right_stage\d_audio_volume';

% passive_workflow = 'lcr_passive';
passive_workflow = 'hml_passive_audio';

% animals = {'DS019','DS020','DS021','HA003','HA004'};

animals = {'AP021','AP022','DS001','DS007','DS010','DS011'}; % V>A
animals = {'DS000','DS003','DS004','DS013','DS014','DS015','DS016'}; % A>V


% Grab V's from passive stimuli on post-learning days
stim_v = cell(size(animals));
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    % Find learned days
    task_recordings = plab.find_recordings(animal,[],task_workflow);
    curr_rxn_stat_p = nan(size(task_recordings));

    for curr_recording = 1:length(task_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = task_recordings(curr_recording).day;
        rec_time = task_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get reaction stat
        use_stat = 'mean';
        curr_rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    % Grab passive LCR response from learned days
    learned_passive_recordings = plab.find_recordings(animal, ...
        {task_recordings(curr_rxn_stat_p < 0.05).day},passive_workflow);

    for curr_recording = 1:length(learned_passive_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = learned_passive_recordings(curr_recording).day;
        rec_time = learned_passive_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % Stim-align data (quiescent trials)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        trial_t_window = [-0.3,1];
        trial_t_rate = 35;
        trial_t = trial_t_window(1):1/trial_t_rate:trial_t_window(2);
        interp_t = reshape(stimOn_times(quiescent_trials),[],1) + trial_t;

        n_wf_components = 200;
        wf_V_stimalign = interp1(wf_t,wf_V(1:n_wf_components,:)',interp_t,'previous');

        switch passive_workflow
            case 'lcr_passive'
                stim_type = vertcat(trial_events.values.TrialStimX);
            case 'hml_passive_audio'
                stim_type = vertcat(trial_events.values.StimFrequence);
        end

        wf_V_stimalign_avg = nan(n_wf_components,length(trial_t),length(unique(stim_type)),class(wf_V_stimalign));
        stim_types_used_idx = ismember(unique(stim_type),unique(stim_type(quiescent_trials)));

        wf_V_stimalign_avg(:,:,stim_types_used_idx) = ...
            permute(ap.groupfun(@mean,wf_V_stimalign,stim_type(quiescent_trials),[],[]),[3,2,1]);

        stim_v{curr_animal}{curr_recording} = wf_V_stimalign_avg;

        ap.print_progress_fraction(curr_recording,length(learned_passive_recordings));
    end

end
disp('Done.');

U_master = plab.wf.load_master_U;

stim_v_avg = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(cat(4,x{:}),4),stim_v,'uni',false),[1,3,4,2])),4);

stim_v_avg = stim_v_avg - nanmean(stim_v_avg(:,trial_t<0,:),2);

stim_v_px = plab.wf.svd2px(U_master(:,:,1:n_wf_components),stim_v_avg);
ap.imscroll(stim_v_px,trial_t);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(AP_colormap('PWG'));
ap.wf_draw('ccf','k');


%% Task-specific: KERNEL passive stim (linear)

% task_workflow = 'stim_wheel_right_stage\d_size_up';
% task_workflow = 'stim_wheel_right_stage\d_angle_size60';
task_workflow = 'stim_wheel_right_stage\d_audio_volume';
% task_workflow = 'stim_wheel_right_stage\d';

% passive_workflow = 'lcr_passive';
passive_workflow = 'hml_passive_audio';

% animals = {'DS019','DS020','DS021','HA003','HA004'};

animals = {'AP021','AP022','DS001','DS007','DS010','DS011'}; % V>A
% animals = {'DS000','DS003','DS004','DS013','DS014','DS015','DS016'}; % A>V


% Grab V's from passive stimuli on post-learning days
stim_v = cell(size(animals));
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    % Find learned days
    task_recordings = plab.find_recordings(animal,[],task_workflow);
    curr_rxn_stat_p = nan(size(task_recordings));

    for curr_recording = 1:length(task_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = task_recordings(curr_recording).day;
        rec_time = task_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get reaction stat
        use_stat = 'mean';
        curr_rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    % Grab passive LCR response from learned days
    learned_passive_recordings = plab.find_recordings(animal, ...
        {task_recordings(curr_rxn_stat_p < 0.05).day},passive_workflow);

    for curr_recording = 1:length(learned_passive_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = learned_passive_recordings(curr_recording).day;
        rec_time = learned_passive_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

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

        frame_shifts = 0:30;
        lambda = 20;

        [kernels,predicted_signals,explained_var] = ...
            ap.regresskernel(wf_V(1:n_components,:),stim_regressors,-frame_shifts,lambda);

        stim_v{curr_animal}{curr_recording} = kernels;

        ap.print_progress_fraction(curr_recording,length(learned_passive_recordings));
    end

end
disp('Done.');

U_master = plab.wf.load_master_U;

stim_v_avg = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(cat(4,x{:}),4),stim_v,'uni',false),[1,3,4,2])),4);

stim_v_px = plab.wf.svd2px(U_master(:,:,1:size(stim_v_avg,1)),stim_v_avg);
ap.imscroll(stim_v_px);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(AP_colormap('PWG'));
ap.wf_draw('ccf','k');



% x = cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:200),x),cat(2,stim_v{:}),'uni',false);
% x2 = nanmedian(cell2mat(permute(x,[1,3,4,5,2])),5);
% ap.imscroll(x2);
% axis image;
% clim(max(abs(clim)).*[-1,1]);
% colormap(AP_colormap('PWG'));
% ap.wf_draw('ccf','k');


%% Task-specific: KERNEL passive stim (logistic)

% task_workflow = 'stim_wheel_right_stage\d_size_up';
% task_workflow = 'stim_wheel_right_stage\d_angle_size60';
task_workflow = 'stim_wheel_right_stage\d_audio_volume';
% task_workflow = 'stim_wheel_right_stage\d';

% passive_workflow = 'lcr_passive';
passive_workflow = 'hml_passive_audio';

% animals = {'DS019','DS020','DS021','HA003','HA004'};

animals = {'AP021','AP022','DS001','DS007','DS010','DS011'}; % V>A
% animals = {'DS000','DS003','DS004','DS013','DS014','DS015','DS016'}; % A>V


% Grab V's from passive stimuli on post-learning days
stim_v = cell(size(animals));
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    % Find learned days
    task_recordings = plab.find_recordings(animal,[],task_workflow);
    curr_rxn_stat_p = nan(size(task_recordings));

    for curr_recording = 1:length(task_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = task_recordings(curr_recording).day;
        rec_time = task_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % (temp: skip if <10 trials - usually this is because short bad one
        % and long good one, not always last - find the good one)
        if n_trials < 10
            continue
        end

        % Get reaction stat
        use_stat = 'mean';
        curr_rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    % Grab passive LCR response from learned days
    learned_passive_recordings = plab.find_recordings(animal, ...
        {task_recordings(curr_rxn_stat_p < 0.05).day},passive_workflow);

    for curr_recording = 1:length(learned_passive_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = learned_passive_recordings(curr_recording).day;
        rec_time = learned_passive_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        switch bonsai_workflow
            case 'lcr_passive'
                stim_type = vertcat(trial_events.values.TrialStimX);
            case 'hml_passive_audio'
                stim_type = vertcat(trial_events.values.StimFrequence);
        end


        switch bonsai_workflow
            case {'lcr_passive','lcr_passive_size60'}
                stim_type = vertcat(trial_events.values.TrialStimX);
                use_stim = 3;
            case 'hml_passive_audio'
                stim_type = vertcat(trial_events.values.StimFrequence);
                use_stim = 2;
        end

        time_bins = [wf_t;wf_t(end)+1/wf_framerate];

        stim_regressors = cell2mat(arrayfun(@(x) ...
            histcounts(stimOn_times(stim_type == x),time_bins), ...
            unique(stim_type),'uni',false));

        n_components = 100;

        timesteps = -10:30;
        regressors = repmat(wf_V(1:n_components,:)',1,1,length(timesteps));
        for i = 1:length(timesteps)
            regressors(:,:,i) = circshift(regressors(:,:,i),-timesteps(i),1);
        end

        regressors_flat = reshape(regressors,[],prod(size(regressors,[2,3])));

        % Model - elastic net (can't get pure ridge alpha=0)
        [B, FitInfo] = lassoglm(regressors_flat, stim_regressors(use_stim,:)','binomial','Alpha',1e-5,'Lambda',1e-2);

        kernels = reshape(B,n_components,[]);

        stim_v{curr_animal}{curr_recording} = kernels;

        ap.print_progress_fraction(curr_recording,length(learned_passive_recordings));
    end

end
disp('Done.');

U_master = plab.wf.load_master_U;

stim_v_avg = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(cat(4,x{:}),4),stim_v,'uni',false),[1,3,4,2])),4);

stim_v_px = plab.wf.svd2px(U_master(:,:,1:size(stim_v_avg,1)),stim_v_avg);
ap.imscroll(stim_v_px);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(AP_colormap('PWG'));
ap.wf_draw('ccf','k');


%% Task-specific: KERNEL passive stim (logistic classification)

% task_workflow = 'stim_wheel_right_stage\d_size_up';
% task_workflow = 'stim_wheel_right_stage\d_angle_size60';
task_workflow = 'stim_wheel_right_stage\d_audio_volume';
% task_workflow = 'stim_wheel_right_stage\d';

% passive_workflow = 'lcr_passive';
passive_workflow = 'hml_passive_audio';

% animals = {'DS019','DS020','DS021','HA003','HA004'};

% animals = {'AP021','AP022','DS001','DS007','DS010','DS011'}; % V>A
animals = {'DS000','DS003','DS004','DS013','DS014','DS015','DS016'}; % A>V


% Grab V's from passive stimuli on post-learning days
stim_v = cell(size(animals));
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    % Find learned days
    task_recordings = plab.find_recordings(animal,[],task_workflow);
    curr_rxn_stat_p = nan(size(task_recordings));

    for curr_recording = 1:length(task_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = task_recordings(curr_recording).day;
        rec_time = task_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % (temp: skip if <10 trials - usually this is because short bad one
        % and long good one, not always last - find the good one)
        if n_trials < 10
            continue
        end

        % Get reaction stat
        use_stat = 'mean';
        curr_rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    % Grab passive LCR response from learned days
    learned_passive_recordings = plab.find_recordings(animal, ...
        {task_recordings(curr_rxn_stat_p < 0.05).day},passive_workflow);

    for curr_recording = 1:length(learned_passive_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = learned_passive_recordings(curr_recording).day;
        rec_time = learned_passive_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        switch bonsai_workflow
            case 'lcr_passive'
                stim_type = vertcat(trial_events.values.TrialStimX);
            case 'hml_passive_audio'
                stim_type = vertcat(trial_events.values.StimFrequence);
        end


        switch bonsai_workflow
            case {'lcr_passive','lcr_passive_size60'}
                stim_type = vertcat(trial_events.values.TrialStimX);
                use_stim = 3;
            case 'hml_passive_audio'
                stim_type = vertcat(trial_events.values.StimFrequence);
                use_stim = 2;
        end


        % max snapshot, lassoglm (doesn't look good)
        [~,~,stim_type_unique] = unique(stim_type);

        stim_window = stimOn_times(stim_type_unique == use_stim) + [0,0.2];
        stim_max_px = cell2mat(arrayfun(@(x) max(plab.wf.svd2px(wf_U, ...
            wf_V(:,wf_t >= stim_window(x,1) ...
            & wf_t <= stim_window(x,2))),[],3), ...
            permute(1:size(stim_window,1),[1,3,2]),'uni',false));

        baseline_window = stimOn_times + [-0.3,0];
        baseline_max_px = cell2mat(arrayfun(@(x) max(plab.wf.svd2px(wf_U, ...
            wf_V(:,wf_t >= baseline_window(x,1) ...
            & wf_t <= baseline_window(x,2))),[],3), ...
            permute(1:size(baseline_window,1),[1,3,2]),'uni',false));

        wf_snapshots = reshape(wf_U,[],size(wf_U,3))'* ...
            cat(2,reshape(stim_max_px,[],size(stim_max_px,3)), ...
            reshape(baseline_max_px,[],size(baseline_max_px,3)));
        stim_shapshot = [true(size(stim_window,1),1);false(size(baseline_window,1),1)];


        n_components = 2000;

        model = fitclinear(wf_snapshots(1:n_components,:)', stim_shapshot, ...
            'Learner', 'logistic', 'Regularization', 'ridge', 'Lambda', 0);

        kernels = model.Beta;

        stim_v{curr_animal}{curr_recording} = kernels;

        ap.print_progress_fraction(curr_recording,length(learned_passive_recordings));
    end

end
disp('Done.');

U_master = plab.wf.load_master_U;

stim_v_avg = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(cat(4,x{:}),4),stim_v,'uni',false),[1,3,4,2])),4);

stim_v_px = plab.wf.svd2px(U_master(:,:,1:size(stim_v_avg,1)),stim_v_avg);
ap.imscroll(stim_v_px);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(AP_colormap('PWG'));
ap.wf_draw('ccf','k');


%% Testing kernel parameters

% % VA V-learn
% animal = 'DS011';
% rec_day = '2024-07-20';
% rec_time = '1755';
% verbose = true;
% ap.load_recording;

% VA A-learn
animal = 'DS011';
rec_day = '2024-07-28';
rec_time = '1042';
verbose = true;
ap.load_recording;


switch bonsai_workflow
    case {'lcr_passive','lcr_passive_size60'}
        stim_type = vertcat(trial_events.values.TrialStimX);
    case 'hml_passive_audio'
        stim_type = vertcat(trial_events.values.StimFrequence);
end

time_bins = [wf_t;wf_t(end)+1/wf_framerate];

stim_regressors = cell2mat(arrayfun(@(x) ...
    histcounts(stimOn_times(stim_type == x),time_bins), ...
    unique(stim_type),'uni',false));

% Regular linear regression (to binary variable, so bad practice)
n_components = 200;

frame_shifts = round(-wf_framerate*0:wf_framerate*1);
lambda = 20;

[kernels,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,:),stim_regressors,-frame_shifts,lambda);

ap.imscroll(plab.wf.svd2px(wf_U(:,:,1:size(kernels,1)),kernels));
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));

%% Logistic regression


n_components = 100;

timesteps = -10:30;
regressors = repmat(wf_V(1:n_components,:)',1,1,length(timesteps));
for i = 1:length(timesteps)
    regressors(:,:,i) = circshift(regressors(:,:,i),-timesteps(i),1);
end

regressors_flat = reshape(regressors,[],prod(size(regressors,[2,3])));

use_stim = 2;

% GLM - not regularized
model = fitglm(regressors_flat, stim_regressors(use_stim,:)','linear','Distribution','binomial','Link','logit');

px = plab.wf.svd2px(wf_U(:,:,1:n_components),reshape(model.Coefficients.Estimate(2:end),n_components,[]));
ap.imscroll(px)


% Model - elastic net (can't get pure ridge alpha=0)
[B, FitInfo] = lassoglm(regressors_flat, stim_regressors(use_stim,:)','binomial','Alpha',1e-5,'Lambda',1e-2);

px = plab.wf.svd2px(wf_U(:,:,1:n_components),reshape(B,n_components,[]));
ap.imscroll(px); axis image;
colormap(ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);


% Model - elastic net (can't get pure ridge alpha=0) (do in chunks?)
n_split = 4;
split_idx = ceil((1:length(wf_t))/(length(wf_t)/n_split));

B = nan(size(regressors_flat,2),n_split,'single');
for curr_split = 1:n_split
    B(:,curr_split) = lassoglm(regressors_flat(split_idx==curr_split,:), ...
        stim_regressors(use_stim,split_idx==curr_split)','binomial','Alpha',1e-5,'Lambda',1e-2);
    ap.print_progress_fraction(curr_split,n_split)
end

px = plab.wf.svd2px(wf_U(:,:,1:n_components),reshape(B,n_components,length(timesteps),[]));
ap.imscroll(median(px,4)); axis image;
colormap(ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);


% Exclude other stim times (didn't change anything)
alt_stim_times = stimOn_times(~ismember(stim_type,90)) + [0,1];
use_t = mod(discretize(wf_t,reshape(alt_stim_times',[],1)),2) ~= 1;

[B, FitInfo] = lassoglm(regressors_flat(use_t,:), stim_regressors(use_stim,use_t)','binomial','Alpha',1e-7,'Lambda',1e-1,'CV',5);

px = plab.wf.svd2px(wf_U(:,:,1:n_components),reshape(B,n_components,[]));
ap.imscroll(px); axis image;


% fitclinear (this looks just exactly like the average?)
model = fitclinear(regressors_flat, stim_regressors(use_stim,:)', 'Learner', 'logistic', 'Regularization', 'ridge', 'Lambda', 1);
px = plab.wf.svd2px(wf_U(:,:,1:n_components),reshape(model.Beta,n_components,[]));
ap.imscroll(px); axis image;
colormap(ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);


% max snapshot, lassoglm (doesn't look good)
[~,~,stim_type_unique] = unique(stim_type);

stim_window = stimOn_times(stim_type_unique == use_stim) + [0,0.2];
stim_max_px = cell2mat(arrayfun(@(x) max(plab.wf.svd2px(wf_U, ...
    wf_V(:,wf_t >= stim_window(x,1) ...
    & wf_t <= stim_window(x,2))),[],3), ...
    permute(1:size(stim_window,1),[1,3,2]),'uni',false));

baseline_window = stimOn_times + [-0.3,0];
baseline_max_px = cell2mat(arrayfun(@(x) max(plab.wf.svd2px(wf_U, ...
    wf_V(:,wf_t >= baseline_window(x,1) ...
    & wf_t <= baseline_window(x,2))),[],3), ...
    permute(1:size(baseline_window,1),[1,3,2]),'uni',false));

wf_snapshots = reshape(wf_U,[],size(wf_U,3))'* ...
    cat(2,reshape(stim_max_px,[],size(stim_max_px,3)), ...
    reshape(baseline_max_px,[],size(baseline_max_px,3)));
stim_shapshot = [true(size(stim_window,1),1);false(size(baseline_window,1),1)];

n_components = 100;
[B, FitInfo] = lassoglm(wf_snapshots(1:n_components,:)', stim_shapshot,'binomial','Alpha',0.1,'Lambda',1);
px = plab.wf.svd2px(wf_U(:,:,1:n_components),B);

figure;
imagesc(px); 
axis image;
colormap(ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);


% max snapshot, fitclinear (doesn't look good) 
n_components = 2000;

model = fitclinear(wf_snapshots(1:n_components,:)', stim_shapshot, ...
    'Learner', 'logistic', 'Regularization', 'ridge', 'Lambda', 0);
px = plab.wf.svd2px(wf_U(:,:,1:n_components),model.Beta);
figure;imagesc(px); axis image;
colormap(ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);













