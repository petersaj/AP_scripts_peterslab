% Load Bonsai data

if verbose; disp('Loading Bonsai...'); end

%% Load general Bonsai events

% Get Bonsai workflow
bonsai_dir = dir(plab.locations.filename('server',animal,rec_day,rec_time,'bonsai'));
bonsai_workflow = bonsai_dir([bonsai_dir.isdir] & ~contains({bonsai_dir.name},'.')).name;

% Load Bonsai events (should be included in every workflow)
bonsai_events_fn = plab.locations.filename('server', ...
    animal,rec_day,rec_time,'bonsai','bonsai_events.csv');

if exist(bonsai_events_fn,'file')
    % Set Bonsai timestamp format
    bonsai_table_opts = detectImportOptions(bonsai_events_fn);
    bonsai_table_opts = setvaropts(bonsai_table_opts,'Timestamp','Type','datetime', ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSZ','TimeZone','local', ...
        'DatetimeFormat','yyyy-MM-dd HH:mm:ss.SSS');

    % Load Bonsai CSV file
    bonsai_events_raw = readtable(bonsai_events_fn,bonsai_table_opts);

    % Check for NaT timestamps, throw warning and flag if so
    if any(isnat(bonsai_events_raw.Timestamp))
        bad_bonsai_csv = true;
        warning('Bonsai file ends improperly: %s',bonsai_events_fn);
    else
        bad_bonsai_csv = false;
    end

    % Create nested structure for trial events
    trial_events = struct('parameters',cell(1),'values',cell(1),'timestamps',cell(1));

    % Save anything in "Trial 0" as a parameter
    parameter_idx = bonsai_events_raw.Trial == 0;
    unique_parameters = unique(bonsai_events_raw.Event(parameter_idx));
    for curr_parameter = unique_parameters'
        curr_parameter_idx = parameter_idx & strcmp(bonsai_events_raw.Event,curr_parameter);
        trial_events.parameters.(cell2mat(curr_parameter)) = bonsai_events_raw.Value(curr_parameter_idx);
    end

    % Loop through trials (excluding 0), save values and timestamps for all events
    % (exclude entries with empty events - happens sometimes on last entry)
    empty_events = cellfun(@isempty,bonsai_events_raw.Event);
    n_trials = max(bonsai_events_raw.Trial);
    unique_events = unique(bonsai_events_raw.Event(~parameter_idx & ~empty_events));
    for curr_trial = 1:n_trials
        for curr_event = unique_events'
            curr_event_idx = bonsai_events_raw.Trial == curr_trial & strcmp(bonsai_events_raw.Event,curr_event);
            trial_events.values(curr_trial).(cell2mat(curr_event)) = bonsai_events_raw.Value(curr_event_idx);
            trial_events.timestamps(curr_trial).(cell2mat(curr_event)) = bonsai_events_raw.Timestamp(curr_event_idx);
        end
    end
end


%% Workflow-specific loading

if contains(bonsai_workflow,'stim_wheel')
    % Task: stim and response times

    % Stim times: when photodiode flips to 1
    stimOn_times = photodiode_times(photodiode_values == 1);

    % Use only trials with outcome
    n_trials = length([trial_events.timestamps.Outcome]);

    % Find the last move stop before stim on
    % (sometimes this isn't after the stimulus: Bonsai's quiescence clock
    % isn't very accurate)
    last_prestim_move_stop = cellfun(@(x) ...
        find(~wheel_move(timelite.timestamps <= x),1,'last') + 1, ...
        num2cell(stimOn_times(1:n_trials)));

    % Find the first move start after the pre-stim move stop
    % (usually this is the first post-stim move, sometimes this is before
    % stim if quiescence clock didn't work)
    stim_move_time = cellfun(@(x) ...
        timelite.timestamps((x-1) + find(wheel_move(x:end),1,'first')), ...
        num2cell(last_prestim_move_stop));

    stim_to_move = stim_move_time - stimOn_times(1:n_trials);

elseif contains(bonsai_workflow,'lcr_passive')
    % Passive protocol stim on times

    % Photodiode bug (old, now fixed): screen could flick to black briefly
    % when clicking on another window. This are always brief, and no way to
    % tell when it happened, so compensate by removing all flips that
    % happen with short duration
    photodiode_flicker = find(diff(photodiode_times) < 0.1);
    if any(photodiode_flicker)
        warning('Photodiode flicker? removing')
        photodiode_times(photodiode_flicker+[0,1]) = [];
        photodiode_values(photodiode_flicker+[0,1]) = [];
    end

    stimOn_times = photodiode_times(photodiode_values == 1);

    % If bad Bonsai CSV: truncate stim times to what was recorded
    if bad_bonsai_csv
        n_bonsai_stim = vertcat(trial_events.values.TrialStimX);
        stimOn_times = stimOn_times(1:n_bonsai_stim);
    end

elseif strcmp(bonsai_workflow,'sparse_noise')
    % Sparse noise: get noise locations and times

    bonsai_noise_fn = plab.locations.filename('server', ...
        animal,rec_day,rec_time,'bonsai','NoiseLocations.bin');
    fid = fopen(bonsai_noise_fn);

    n_x_squares = trial_events.parameters.ScreenExtentX./trial_events.parameters.StimSize;
    n_y_squares = trial_events.parameters.ScreenExtentY./trial_events.parameters.StimSize;

    noise_locations = reshape(fread(fid),n_y_squares,n_x_squares,[]);
    fclose(fid);

    % Get stim times from photodiode (extrapolate: sparse noise photodiode
    % flips every N stim to give a more robust signal)
    if ~isfield(trial_events.parameters,'NthPhotodiodeFlip')
        % (if it wasn't defined, default to flipping on every stim)
        trial_events.parameters.NthPhotodiodeFlip = 1;
    end
    photodiode_stim_idx = 1:trial_events.parameters.NthPhotodiodeFlip:size(noise_locations,3);
    % (check that the number of photodiode flips is expected)
    if length(photodiode_stim_idx) ~= length(photodiode_times)
        if length(photodiode_stim_idx) < length(photodiode_times)
            % (rarely: Bonsai square to black temporarily, don't know why)
            % (fix? find when time differences on either side are less than a
            % threshold and remove those flips)
            photodiode_diff = diff(photodiode_times);
            photodiode_diff_thresh = mean(photodiode_diff)/2;
            bad_photodiode_idx = ...
                find(photodiode_diff(1:end-1) < photodiode_diff_thresh & ...
                photodiode_diff(2:end) < photodiode_diff_thresh) + 1;

            if (length(photodiode_times) - length(photodiode_stim_idx)) == ...
                    length(bad_photodiode_idx)
                % (if detected bad flips even out the numbers, apply fix)
                photodiode_times(bad_photodiode_idx) = [];
            else
                % (otherwise, error)
                error('Sparse noise: photodiode > stim, unfixable')
            end
        else
            error('Sparse noise: photodiode < stim')
        end
    end

    stim_times = interp1(photodiode_stim_idx,photodiode_times, ...
        1:size(noise_locations,3),'linear','extrap')';

end


%% For testing: convert bonsai times to timelite times

% %%%% testing: bonsai times into timelite times
% % % (reward)
% % timelite2bonsai = reward_times;
% % bonsai2timelite_datetime = vertcat(trial_events.timestamps([trial_events.values.Outcome] == 1).Outcome);
% % (stim)
% timelite2bonsai = stimOn_times;
% bonsai2timelite_datetime = cellfun(@(x) x(1),{trial_events.timestamps.StimOn}');
% bonsai_relative_t = bonsai2timelite_datetime(1);
% bonsai2timelite = seconds(bonsai2timelite_datetime - bonsai_relative_t);
%
% % event_datetime = cellfun(@(x) x(1),{trial_events.timestamps.StimOn}');
% event_datetime = vertcat(trial_events.timestamps.QuiescenceStart);
% % event_datetime = vertcat(trial_events.timestamps.QuiescenceReset);
%
% event_t_bonsai = seconds(event_datetime - bonsai_relative_t);
%
% event_t_tl = interp1(bonsai2timelite,timelite2bonsai,event_t_bonsai,'linear','extrap');









