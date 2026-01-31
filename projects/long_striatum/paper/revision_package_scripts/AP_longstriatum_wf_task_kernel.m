%% Load and package task/passive widefield kernels

save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Marica_2025\data';

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

data_all = cell(length(animals),1);

for animal_idx=1:length(animals)

    animal = animals{animal_idx};

    % Find all task recordings
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);

    data_animal = table;

    for use_rec=1:length(recordings_task)

        rec_day = recordings_task(use_rec).day;
        load_parts.behavior = true;
        load_parts.widefield = true;
        load_parts.widefield_master = true;

        % Set regression parameters
        n_components = 200;
        frame_shifts = -10:30;
        lambda = 20;
        cv_fold = 1;
        skip_frames = 300; % frames start/end to skip for artifacts

        % Load task, regress widefield to stim, store
        rec_time = recordings_task(use_rec).recording{end};
        ap.load_recording

        time_bins = [wf_t;wf_t(end)+1/wf_framerate];

        stim_regressors = histcounts(stimOn_times,time_bins);

        [task_kernels,predicted_signals,explained_var] = ...
            ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
            stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

        % Load passive, regress widefield to stim, store
        workflow_passive = {'lcr_passive'};
        rec_time = plab.find_recordings(animal, rec_day, workflow_passive).recording{end};
        ap.load_recording

        time_bins = [wf_t;wf_t(end)+1/wf_framerate];

        stim_x = vertcat(trial_events.values.TrialStimX);

        stim_regressors = cell2mat(arrayfun(@(x) ...
            histcounts(stimOn_times(stim_x(1:length(stimOn_times)) == x),time_bins), ...
            unique(stim_x),'uni',false));

        [passive_kernels,predicted_signals,explained_var] = ...
            ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
            stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

        % Save data in table
        data_animal.animal(use_rec) = {animal};
        data_animal.rec_day(use_rec) = {rec_day};

        data_animal.kernel_frames(use_rec) = {frame_shifts};
        data_animal.task_kernels(use_rec) = {task_kernels};
        data_animal.passive_kernels(use_rec) = {passive_kernels};

        disp(['Done day ' num2str(use_rec)])
        
    end
    
    % Add current animal to full dataset
    data_all{animal_idx} = data_animal;
    disp(['Done ' animal])

end

% Concatenate data into one table and save
wf_kernels = vertcat(data_all{:});
save_name = fullfile(save_path, 'wf_kernels');
save(save_name, "wf_kernels", "-v7.3");

