% Batch package/analyze HA stim lick data

%% Set animals

animals = [ ...
    "HA007", "HA008", "HA009","HA010", "HA011", "HA012", "HA013", "HA014", "HA015", ... % first dataset
    "AP030","AP031","AP032", ... % never learned static
    "HA016", "HA017", "HA018", ... % small to big
    "HA019", "HA020", ... % start big
    "AP036", "AP037"];

data_path = "C:\Users\petersa\Documents\PetersLab\analysis\stim_lick\data";
save_filename = fullfile(data_path,'animals');
save(save_filename,'animals');
fprintf('Saved: %s\n',save_filename);


%% Package: behavior

disp('Packaging: behavior');

% Get animals
data_path = "C:\Users\petersa\Documents\PetersLab\analysis\stim_lick\data";
load(fullfile(data_path,'animals.mat'));

% Loop through animals and package
bhv = struct;
for curr_animal = 1:length(animals)

    animal = animals(curr_animal);
    fprintf('Grabbing: %s\n',animal)

    animal_recordings = plab.find_recordings(animal,[],'*lick_two_stim*');
    wf_recordings = cellfun(@(x) x(end),{animal_recordings.widefield});
    use_recordings = animal_recordings(wf_recordings);

    for curr_recording = 1:length(use_recordings)

        % Grab pre-load variables
        preload_vars = who;

        rec_day = use_recordings(curr_recording).day;
        rec_time = use_recordings(curr_recording).recording{end};

        load_parts.behavior = true;
        try
            ap.load_recording;
        catch me
            warning('LOAD ERROR %s %s: %s',animal,rec_day,me.message)
            continue
        end

        n_trials = sum(cellfun(@(x) length(x) == 2,{trial_events.timestamps.StimOn}));

        % Get trial parameters
        trial_stim_x = vertcat(trial_events.values(1:n_trials).TrialX);
        trial_static_stim_time = vertcat(trial_events.values(1:n_trials).TrialStimStaticTime);
        trial_quiescence_time = vertcat(trial_events.values(1:n_trials).TrialQuiescence);

        % Lick PSTHs
        lick_window = [-7,10];
        lick_binsize = 0.01;

        if contains(bonsai_workflow,'move')
            % Moving stim

            % Get stim on/move times
            % (photodiode CS- = on/off for CS+ = on,pulse on move,off)
            stim_pd_n = (trial_stim_x == -90)*1 + (trial_stim_x == 90)*2;
            stim_pd_on_grouped = mat2cell(photodiode_on_times(1:sum(stim_pd_n)),stim_pd_n);
            stim_pd_off_grouped = mat2cell(photodiode_off_times(1:sum(stim_pd_n)),stim_pd_n);

            stimOn_times = cellfun(@(x) x(1), stim_pd_on_grouped);
            stim_move_times = cellfun(@(x) x(1), stim_pd_off_grouped);

            % PSTHs for stim on/move (CS- would-be)
            [~,lick_raster_stim_on,lick_t] = ...
                ap.psth(lick_times,stimOn_times, ...
                'window',lick_window,'bin_size',lick_binsize);
            [~,lick_raster_stim_move] = ...
                ap.psth(lick_times,stim_move_times, ...
                'window',lick_window,'bin_size',lick_binsize);

        elseif contains(bonsai_workflow,'static')
            % Static stim

            % Get stim on/reward available times
            stimOn_times = photodiode_on_times(1:n_trials);          

            % PSTHs for stim on
            [~,lick_raster_stim_on,lick_t] = ...
                ap.psth(lick_times,stimOn_times, ...
                'window',lick_window,'bin_size',lick_binsize);
        end

        % Store variables
        bhv(curr_animal,curr_recording).animal = animal;
        bhv(curr_animal,curr_recording).rec_day = rec_day;
        bhv(curr_animal,curr_recording).rec_time = rec_time;
        bhv(curr_animal,curr_recording).workflow = bonsai_workflow;

        bhv(curr_animal,curr_recording).trial_stim_x = trial_stim_x;
        bhv(curr_animal,curr_recording).trial_static_stim_time = trial_static_stim_time;
        bhv(curr_animal,curr_recording).trial_quiescence_time = trial_quiescence_time;

        % (convert lick raster to binary for compression)
        bhv(curr_animal,curr_recording).lick_t = lick_t;
        bhv(curr_animal,curr_recording).lick_raster_stim_on = sparse(logical(lick_raster_stim_on));
        if exist('lick_raster_stim_move','var')
            bhv(curr_animal,curr_recording).lick_raster_stim_move = lick_raster_stim_move;
        end

        % Clear load variables
        clearvars('-except',preload_vars{:});

        ap.print_progress_fraction(curr_recording,length(use_recordings));
    end
    fprintf('Done: %s\n',animal)
end

% Save
data_path = "C:\Users\petersa\Documents\PetersLab\analysis\stim_lick\data";
save_filename = fullfile(data_path,'bhv');
save(save_filename,'bhv');


%% Package: task stim kernels

disp('Packaging: task stim kernels');

% Get animals
data_path = "C:\Users\petersa\Documents\PetersLab\analysis\stim_lick\data";
load(fullfile(data_path,'animals.mat'));

% Loop through animals and package
wf_task = struct;
for curr_animal = 1:length(animals)

    animal = animals(curr_animal);
    fprintf('Grabbing: %s\n',animal)

    animal_recordings = plab.find_recordings(animal,[],'*lick_two_stim*');
    wf_recordings = cellfun(@(x) x(end),{animal_recordings.widefield});
    use_recordings = animal_recordings(wf_recordings);

    for curr_recording = 1:length(use_recordings)

        % Grab pre-load variables
        preload_vars = who;

        rec_day = use_recordings(curr_recording).day;
        rec_time = use_recordings(curr_recording).recording{end};

        load_parts.widefield = true;
        load_parts.widefield_master = true;
        try
            ap.load_recording;
        catch me
            warning('LOAD ERROR %s %s: %s',animal,rec_day,me.message)
            continue
        end

        n_trials = sum(cellfun(@(x) length(x) == 2,{trial_events.timestamps.StimOn}));

        % Set parameters for regression
        time_bins = [wf_t;wf_t(end)+1/wf_framerate];
        n_components = 200;
        frame_shifts = -5:40;
        lambda = 5;
        cv_fold = 5;

        skip_t = 60; % seconds start/end to skip for artifacts
        skip_frames = round(skip_t*wf_framerate);
     
        % Get stim on/move times  
        trial_stim_x = vertcat(trial_events.values(1:n_trials).TrialX);

        if contains(bonsai_workflow,'move')
            % Moving stim
            % (photodiode CS- = on/off for CS+ = on,pulse on move,off)
            stim_pd_n = (trial_stim_x == -90)*1 + (trial_stim_x == 90)*2;
            stim_pd_on_grouped = mat2cell(photodiode_on_times(1:sum(stim_pd_n)),stim_pd_n);
            stim_pd_off_grouped = mat2cell(photodiode_off_times(1:sum(stim_pd_n)),stim_pd_n);

            stimOn_times = cellfun(@(x) x(1), stim_pd_on_grouped);
            stim_move_times = cellfun(@(x) x(1), stim_pd_off_grouped);

        elseif contains(bonsai_workflow,'static')
            % Static stim
            stimOn_times = photodiode_on_times(1:n_trials);          
        end

        % Set regressors
        stim_regressors = zeros(0,length(time_bins)-1);
        stim_regressors(1,:) = histcounts(stimOn_times(trial_stim_x == 90),time_bins);
        stim_regressors(2,:) = histcounts(stimOn_times(trial_stim_x == -90),time_bins);
        if contains(bonsai_workflow,'move')
            stim_regressors(3,:) = histcounts(stim_move_times(trial_stim_x == 90),time_bins);
        end
      
        % Get stim kernel
        [kernels,predicted_signals,explained_var] = ...
            ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
            stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

        % Store variables
        wf_task(curr_animal,curr_recording).animal = animal;
        wf_task(curr_animal,curr_recording).rec_day = rec_day;
        wf_task(curr_animal,curr_recording).rec_time = rec_time;
        wf_task(curr_animal,curr_recording).workflow = bonsai_workflow;

        wf_task(curr_animal,curr_recording).frame_shifts = frame_shifts;
        wf_task(curr_animal,curr_recording).stim_kernel = kernels;

        % Clear load variables
        clearvars('-except',preload_vars{:});

        ap.print_progress_fraction(curr_recording,length(use_recordings));
    end
    fprintf('Done: %s\n',animal)
end

% Save
save_filename = fullfile(data_path,'wf_task');
save(save_filename,'wf_task');


%% Package: passive stim kernels

disp('Packaging: passive stim kernels');

% Get animals
data_path = "C:\Users\petersa\Documents\PetersLab\analysis\stim_lick\data";
load(fullfile(data_path,'animals.mat'));

% Loop through animals and package
wf_passive = struct;
for curr_animal = 1:length(animals)

    animal = animals(curr_animal);
    fprintf('Grabbing: %s\n',animal)

    animal_recordings = plab.find_recordings(animal,[],'lcr_passive*');
    wf_recordings = cellfun(@(x) x(end),{animal_recordings.widefield});
    use_recordings = animal_recordings(wf_recordings);

    for curr_recording = 1:length(use_recordings)

        % Grab pre-load variables
        preload_vars = who;

        rec_day = use_recordings(curr_recording).day;
        rec_time = use_recordings(curr_recording).recording{end};

        load_parts.widefield = true;
        load_parts.widefield_master = true;
        try
            ap.load_recording;
        catch me
            warning('LOAD ERROR %s %s: %s',animal,rec_day,me.message)
            continue
        end

        n_trials = sum(cellfun(@(x) length(x) == 2,{trial_events.timestamps.StimOn}));

        % Set parameters for regression
        time_bins = [wf_t;wf_t(end)+1/wf_framerate];
        n_components = 200;
        frame_shifts = -5:40;
        lambda = 5;
        cv_fold = 5;

        skip_t = 60; % seconds start/end to skip for artifacts
        skip_frames = round(skip_t*wf_framerate);
     
        % Get stim kernels
        stim_x = vertcat(trial_events.values.TrialStimX);
        stim_x_unique = unique(stim_x);

        stim_regressors = cell2mat(arrayfun(@(x) ...
            histcounts(stimOn_times(stim_x == x),time_bins), ...
            stim_x_unique,'uni',false));
      
        [kernels,predicted_signals,explained_var] = ...
            ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
            stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

        % Store variables
        wf_passive(curr_animal,curr_recording).animal = animal;
        wf_passive(curr_animal,curr_recording).rec_day = rec_day;
        wf_passive(curr_animal,curr_recording).rec_time = rec_time;
        wf_passive(curr_animal,curr_recording).workflow = bonsai_workflow;

        wf_passive(curr_animal,curr_recording).stim_x_unique = stim_x_unique;
        wf_passive(curr_animal,curr_recording).frame_shifts = frame_shifts;
        wf_passive(curr_animal,curr_recording).stim_kernel = kernels;

        % Clear load variables
        clearvars('-except',preload_vars{:});

        ap.print_progress_fraction(curr_recording,length(use_recordings));
    end
    fprintf('Done: %s\n',animal)
end

% Save
save_filename = fullfile(data_path,'wf_passive');
save(save_filename,'wf_passive');



