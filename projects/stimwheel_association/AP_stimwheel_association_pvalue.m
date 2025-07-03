function [p,rxn_stat,rxn_null_stat,stim_to_move_nullmean] = AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_move,use_stat)
% [p,rxn_stat,rxn_null_stat,stim_to_move_null_mean] = AP_stimwheel_association_pvalue(stimOn_times,trial_events,stim_to_move)
% Get p-value for whether reaction times are faster than chance (the animal
% has a stim-wheel association)
%
% Does this by conditional resampling: determines all valid stim onset
% times for each trial if other quiescence periods had been selected, gets
% what the reaction times would have been for each of those valid stim
% times, compares to actual reaction times.
%
% Inputs:
% stimOn_times: stim times from Timelite
% trial_events: Bonsai trial event structure from ap.load_bonsai
% stim_to_move: measured reaction times
% use_stat (optional): reaction statistic to use, 'mad' (default), 'mean', 'median'
%
% Output:
% p: p-value for median reaction time
% rxn_stat: statistic for reaction times
% rxn_null_stat: statistic for null reaction times
% stim_to_move_nullmean: mean null reaction time for each trial

% Only look at completed trials (all relevant events occured)
n_trials = min([length(stimOn_times), ...
    length(stim_to_move), ...
    length([trial_events.timestamps.Outcome])]);

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

    % Get time of last quiescence reset (Timelite)
    % (define as stim time - trial quiescence time)
    last_quiescence_reset = stimOn_times(curr_trial) - trial_events.values(curr_trial).TrialQuiescence;

    % Get quiescence reset times (Timelite)
    curr_quiescence_resets_bonsai = vertcat(...
        trial_events.timestamps(curr_trial).QuiescenceStart, ... % start of quiescence period
        trial_events.timestamps(curr_trial).QuiescenceReset);    % all quiescence resets

    curr_quiescence_resets_timelite = last_quiescence_reset + ...
        seconds(curr_quiescence_resets_bonsai - ...
        curr_quiescence_resets_bonsai(end,:));

    % Get quiescence durations between resets, including post-stim move
    curr_poststim_move_timelite = stimOn_times(curr_trial) + stim_to_move(curr_trial);

    curr_quiescence_durations = ...
        diff(vertcat(curr_quiescence_resets_timelite, ...
        curr_poststim_move_timelite));

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

    % (debugging: plot trial)
    if false

        trial_start_t = stimOff_times(curr_trial-1);
        trial_end_t = stimOn_times(curr_trial) + 4;
        curr_t_idx = timelite.timestamps > trial_start_t & timelite.timestamps < trial_end_t;

        figure; hold on;
        % (wheel velocity)
        plot(timelite.timestamps(curr_t_idx),wheel_velocity(curr_t_idx));
        % (quiescence resets)
        plot(curr_quiescence_resets_timelite,0,'|k')
        % (stim onset)
        xline(stimOn_times(curr_trial),'g','linewidth',2);
        % (move onset)
        xline(stimOn_times(curr_trial)+stim_to_move(curr_trial),'r','linewidth',2);
        % (possible stim times)
        xline(stimOn_times_valid{curr_trial},'--m');
        % (last quiescence reset)
        xline(stimOn_times(curr_trial)-trial_events.values(curr_trial).TrialQuiescence,'k')

    end

end

move_times = stimOn_times(1:n_trials) + stim_to_move;
stim_to_move_valid = cellfun(@(stim_time,move_time) move_time-stim_time, ...
    stimOn_times_valid,num2cell(move_times),'uni',false);

% Create null reaction distribution 
% (only from trials with reaction times > 0ms: negative reaction times mean
% the quiescence time wasn't working exactly)
stim_rxn_threshold = -Inf;
null_use_trials = (stim_to_move > stim_rxn_threshold) & ~cellfun(@isempty,stim_to_move_valid);

n_samples = 10000;
stim_to_move_null = nan(length(stim_to_move),n_samples);
stim_to_move_null(null_use_trials,:) = ...
    cell2mat(cellfun(@(x) datasample(x,n_samples)', ...
    stim_to_move_valid(null_use_trials),'uni',false));

stim_to_move_nullmean = mean(stim_to_move_null,2);

% Get reaction statistic
% (default is mad)
if nargin < 4 || isempty(use_stat)
    use_stat = 'mad';
end

% %%%%%% TESTING: NAN-OUT SHORT TRIALS
% rxn_cutoff = 0.1;
% stim_to_move = stim_to_move.*AP_nanout(stim_to_move < rxn_cutoff);
% stim_to_move_null = stim_to_move_null.*AP_nanout(stim_to_move_null < rxn_cutoff);
% %%%%%%%%%%%%%%

switch use_stat
    case 'mad'
        rxn_stat = mad(stim_to_move(null_use_trials),1,1);
        rxn_null_stat_distribution = mad(stim_to_move_null(null_use_trials,:),1,1);

    case 'std'
        rxn_stat = std(stim_to_move(null_use_trials),[],1,'omitmissing');
        rxn_null_stat_distribution = std(stim_to_move_null(null_use_trials,:),[],1,'omitmissing');

    case 'mean'
        rxn_stat = mean(stim_to_move(null_use_trials),1,'omitmissing');
        rxn_null_stat_distribution = mean(stim_to_move_null(null_use_trials,:),1,'omitmissing');

    case 'median'
        rxn_stat = median(stim_to_move(null_use_trials),1,'omitmissing');
        rxn_null_stat_distribution = median(stim_to_move_null(null_use_trials,:),1,'omitmissing');

end

rxn_stat_rank = tiedrank(horzcat(rxn_stat,rxn_null_stat_distribution));
p = rxn_stat_rank(1)./(n_samples+1);

% Output mean statistic of null reaction times
rxn_null_stat = mean(rxn_null_stat_distribution);






