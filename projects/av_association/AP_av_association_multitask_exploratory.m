%% Da's secondary-ish project with multiple visual stimuli



%%

% % Get all animals used in the size-up task
% size_up_dir = dir(fullfile(plab.locations.server_data_path,'**','bonsai','**','*size_up*.bonsai'));
% animals = unique(extractBetween({x.folder}, ...
%     [plab.locations.server_data_path,filesep], ...
%     filesep));

animals = {'DS019','DS020','DS021','HA003','HA004'};

% Grab V's from passive stimuli on post-learning days
avg_v = cell(size(animals));
for curr_animal = 1:length(animals)

    % Find learned days
    curr_recordings = plab.find_recordings(animals{curr_animal},[],'*size*');
    curr_rxn_stat_p = nan(size(curr_recordings));

    for curr_recording = 1:length(recordings)
        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get reaction stat
        use_stat = 'mean';
        curr_rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

        ap.print_progress_fraction(curr_recording,length(recordings));
    end

    % Grab passive LCR response from learned days
    learned_passive_recordings = plab.find_recordings(animal, ...
        {curr_recordings(curr_rxn_stat_p < 0.05).day},'lcr_passive');
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
        interp_t = reshape(stimOn_times,[],1) + trial_t;

        n_wf_components = 200;
        wf_V_stimalign = interp1(wf_t,wf_V(1:n_wf_components,:)',interp_t,'previous');

    end

end














