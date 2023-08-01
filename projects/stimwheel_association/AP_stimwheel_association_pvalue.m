function p = AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_move)
% p = AP_stimwheel_association_pvalue(trial_events)
% Get p-value for whether reaction times are faster than chance (the animal
% has a stim-wheel association)
%
% Does this by conditional resampling: determines all valid stim on times
% for each trial if other quiescence periods had been selected, gets what
% the reaction times would have been for each of those valid stim times,
% compares to actual reactin times.
%
% Inputs:
% stimOn_times: stim times from Timelite
% trial_events: Bonsai trial event structure from ap.load_bonsai
% stim_to_move: measured reaction times
%
% Output:
% p: p-value for median reaction time

% Only look at completed trials
n_trials = length([trial_events.timestamps.Outcome]);

% Get quiescence range from Bonsai parameters
quiescence_range = trial_events.parameters.QuiescenceTimes;

% Loop through trials, get possible valid stim times
stimOn_times_valid = cell(n_trials,1);
for curr_trial = 1:n_trials

    % Old bug: before trial 1 delay, sometimes first trial quiescence
    % wasn't saved properly. If no trial quiescence, skip trial.
    if isempty(trial_events.values(curr_trial).TrialQuiescence)
        continue
    end

    % NOTE: this is only using alternate quiescence periods, not alternate
    % ITIs. No quiescence clock during ITI means no direct Bonsai measure
    % of would-be quiescence resets, and they're not accurately estimatable
    % from NIDAQ because of timing/precision differences. This means that
    % there's fewer "valid" stim times, so less statistical power, and may
    % err on the side of missing learned days.

    % Get quiescence durations
    response_move_timestamp = trial_events.timestamps(curr_trial).StimOn(1) + ...
        seconds(stim_to_move(curr_trial));

    curr_quiescence_resets = vertcat(...
        trial_events.timestamps(curr_trial).QuiescenceStart, ... % start of quiescence period
        trial_events.timestamps(curr_trial).QuiescenceReset, ... % all quiescence resets
        response_move_timestamp);                                % first post-stim movement

    curr_quiescence_durations = seconds(diff(curr_quiescence_resets));

    % Get valid quiescence times that would yield same response movement
    % (i.e. last quiescence duration was first over-threshold)
    quiescence_overthresh_grid = curr_quiescence_durations > quiescence_range';
    valid_quiescence_times = quiescence_range( ...
        (sum(quiescence_overthresh_grid,1) == 1) & quiescence_overthresh_grid(end,:));

    % Get valid stim offsets (timelite/phodiode)
    % (get offsets between actual and valid quiescence times, apply to
    % actual stim times to get all valid stim times)
    valid_quiescence_offsets =  ...
        valid_quiescence_times - trial_events.values(curr_trial).TrialQuiescence;

    stimOn_times_valid{curr_trial} = stimOn_times(curr_trial) + valid_quiescence_offsets;

end

move_times = stimOn_times(1:n_trials) + stim_to_move;
stim_to_move_valid = cellfun(@(stim_time,move_time) move_time-stim_time, ...
    stimOn_times_valid,num2cell(move_times),'uni',false);

% Create null reaction distribution (only from trials with reaction
% times > 0ms: these are bad quiescence clock)
null_use_trials = (stim_to_move > 0) & ~cellfun(@isempty,stim_to_move_valid);

n_samples = 10000;
stim_to_move_null = nan(length(stim_to_move),n_samples);
stim_to_move_null(null_use_trials,:) = ...
    cell2mat(cellfun(@(x) datasample(x,n_samples)', ...
    stim_to_move_valid(null_use_trials),'uni',false));

% Get reaction statistic (median)
rxn_stat = nanmedian(stim_to_move(null_use_trials),1);
rxn_null_stat = nanmedian(stim_to_move_null(null_use_trials,:),1);

rxn_stat_rank = tiedrank(horzcat(rxn_stat,rxn_null_stat));
p = rxn_stat_rank(1)./(n_samples+1);




