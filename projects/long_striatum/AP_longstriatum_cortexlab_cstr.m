%% Corticostriatal imaging during operant learning from cortex lab
% Preprocess, load, plot data
%
% Only 3 animals learned operant, and their expression patterns are very
% variable, so probably not enough to say anything

%% TESTING: get null times by ITI, rather than quiescence

% (this might not work, because you can only use smaller ITIs, and if
% there's a change in spontaneous move rate with time that could give
% difference)

%%% !! might not need to do this - ITIs taken into account anyway. So just
%%% need to set "possible quiescence times" as fixed value

alt_stimOn_times = cell(n_trials,1);
alt_stimOn_trialparams = cell(n_trials,1);
% (skip trial 1: no ITI and bad quiescence watch)
for curr_trial = 2:n_trials

    %%%%%%%%% STIM RESPONSE VS. NULL

    t = Timeline.rawDAQTimestamps';

    % Get quiescence reset threshold
    % (only stored in script! whoops, this was dumb)
    expDef_fn = [fileparts(block_filename) filesep day '_' ...
        num2str(experiment) '_' animal '_expDef.m'];
    if ~exist(expDef_fn,'file')
        error('%s %s: no expDef.m',animal,day)
    end
    expDef_text = fileread(expDef_fn);
    [~,quiescThreshold_txt] = regexp(expDef_text, ...
        'quiescThreshold = (\d*)','match','tokens');
    quiescThreshold = str2num(quiescThreshold_txt{1}{1});

    % Get quiescence reset times
    %             % (real: from new trial onset)
    %             quiescence_reset_t = AP_find_stimWheel_quiescence;
    % (extrapolated: from response to response)
    quiescence_reset_t_extrap = AP_extrap_stimWheel_quiescence;

    alt_stimOn_times = cell(n_trials,1);
    alt_stimOn_trialparams = cell(n_trials,1);
    % (skip trial 1: no ITI and bad quiescence watch)
    for curr_trial = 2:n_trials

        % Pull out current trial times (last to next response)
        curr_trial_t_idx = t >= signals_events.responseTimes(curr_trial-1) & ...
            t <= signals_events.responseTimes(curr_trial);
        curr_trial_t = t(curr_trial_t_idx);

        % Get quiescence reset times for different ITIs
        param_timestep = 0.1; % (hardcoded in expDef)
        possible_iti = max([block.paramsValues.itiMin]):param_timestep:max([block.paramsValues.itiMax]);
        possible_quiescence = max([block.paramsValues.quiescenceMin]):param_timestep:max([block.paramsValues.quiescenceMax]);

        t_from_quiescence_reset_trialitis = nan(length(curr_trial_t),length(possible_iti));
        for curr_possible_iti = 1:length(possible_iti)
            curr_possible_itiend = signals_events.responseTimes(curr_trial-1) + ...
                possible_iti(curr_possible_iti);
            curr_quiescence_resets = sort([quiescence_reset_t_extrap;curr_possible_itiend]);

            t_from_quiescence_reset_trialitis(:,curr_possible_iti) = ...
                curr_trial_t - interp1(curr_quiescence_resets, ...
                curr_quiescence_resets,curr_trial_t,'previous','extrap');
        end

        % Find alternate stim times which would have given same first move

        % (getting possible iti + quiescence crosses)
        alt_iti_reached = ((t(curr_trial_t_idx) - curr_trial_t(1)) > possible_iti);
        alt_quiescence_reached = ...
            t_from_quiescence_reset_trialitis > permute(possible_quiescence,[1,3,2]);

        % (get possible stim times as iti x quiescence grid)
        [alt_stim_value,alt_stim_idx] = max( ...
            permute(alt_iti_reached & alt_quiescence_reached,[2,3,1]),[],3);
        alt_stim_t = curr_trial_t(alt_stim_idx);
        alt_stimOn_times_all = alt_stim_t(alt_stim_value);

        % (get alt stim times that would have resulted in the same
        % first movement since that's the measured value)
        stim_leeway = 0.1;
        curr_wheel_move_alt_stim_idx = ...
            arrayfun(@(stim) find(wheel_starts > stim-stim_leeway,1,'first'), ...
            alt_stimOn_times_all);
        use_alt_stimOn_times = ...
            curr_wheel_move_alt_stim_idx == wheel_move_stim_idx(curr_trial);

        % (make sure that real parameters = real stim time: missed
        % wheel clicks sometimes give non-reproducible traces, in
        % which case the trial shouldn't be used for stats)
        curr_quiescence_idx = find(possible_quiescence == quiescence_t(curr_trial));
        curr_iti_idx = find(possible_iti == iti_t(curr_trial-1));
        curr_block_stimOn = signals_events.stimOnTimes(curr_trial);
        curr_alt_stim_offset = curr_block_stimOn - ...
            alt_stim_t(curr_iti_idx,curr_quiescence_idx);
        if curr_alt_stim_offset > 0.01
            continue
        end

        % (apply the block vs actual stim on time delay to all
        % times - note this is regularly backwards in time??)
        curr_stim_pd_offset = stimOn_times(curr_trial) - curr_block_stimOn;
        alt_stimOn_times_all_pd = alt_stimOn_times_all + curr_stim_pd_offset;

        % (store alternate stim times)
        alt_stimOn_times{curr_trial} = alt_stimOn_times_all_pd(use_alt_stimOn_times);

        % (trial plot)
        figure; hold on;
        t_plot_scale = 0.1;
        plot(t(curr_trial_t_idx),wheel_velocity(curr_trial_t_idx),'k')
        plot(t(curr_trial_t_idx),[0;diff(wheel_position(curr_trial_t_idx))]*0.1,'g')
        plot(alt_stimOn_times_all,0,'ob');
        plot(alt_stimOn_times{curr_trial},0,'.r');
        line(repmat(curr_trial_t(1)+signals_events.trialITIValues(curr_trial-1),2,1),ylim);
        line(xlim,repmat(signals_events.trialQuiescenceValues(curr_trial),2,1)*t_plot_scale,'color','m');
        line(repmat(curr_block_stimOn,1,2),ylim,'color','r','linestyle','--');
        line(repmat(stimOn_times(curr_trial),1,2),ylim,'color','k','linestyle','--');
        drawnow;

    end
end





%% Grab and save behavior (operant)

animal_group = 'cstr_fixedq';
% animals = {'AP092','AP093','AP094','AP095','AP096','AP097'}; % standard task, only 3 learned well?]
% animals = {'AP092','AP094','AP096'}; % 3 that learned above
% animals = {'AP093','AP095','AP097'}; % 3 that didn't learn above
animals = {'AP089','AP090','AP091'}; % fixed-quiescence

protocol = 'AP_stimWheelRight';
flexible_name = false;
bhv = struct;

% Flags
plot_flag = true;
save_bhv = false;

for curr_animal = 1:length(animals)
    
    preload_vars = who;
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments_externalhd(animal,protocol,flexible_name);
    
    if isempty(experiments)
        disp(['No behavior data: ' animal]);
        continue
    end
    
    disp([animal ', day:'])
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        
        % If multiple experiments, only use the last one
        % (usually multiple happens if mess ups and final one is good)
        for curr_experiment = length(experiment_num)
            
            experiment = experiment_num(curr_experiment);
            
            % Load (behavior/timeline only)
            load_parts.imaging = false;
            AP_load_experiment_externalhd;
            
            % Get protocol name
            [~,curr_protocol] = fileparts(block.expDef);
            
            % Time of session (in minutes)curr_day
            session_duration = block.duration/60;
            
            % Trial counts (use only full trials with a response)
            total_water = sum(block.outputs.rewardValues);
            
            % Hit fraction
            frac_hit = nanmean(trial_outcome == 1);
            
            % (probability of any movement around stim)
            stim_surround_t_centers = -10:0.1:10;
            stim_surround_times = stimOn_times + stim_surround_t_centers;
            stim_surround_move = interp1(Timeline.rawDAQTimestamps, ...
                +wheel_move,stim_surround_times,'previous');
            
            % Resetting quiescence period for each trial
            if isfield(block.events,'trialQuiescenceValues')
                quiescence_t = block.events.trialQuiescenceValues(1:n_trials);
            else
                quiescence_t = 0.5*ones(size(block.events.stimOnTimes));
            end

            % Get parameters only stored in script
            % (didn't save as param in signals)
            expDef_fn = [fileparts(block_filename) filesep day '_' ...
                num2str(experiment) '_' animal '_expDef.m'];
            if ~exist(expDef_fn,'file')
                error('%s %s: no expDef.m',animal,day)
            end
            expDef_text = fileread(expDef_fn);
            % - quiescence reset threshold
            [~,quiescThreshold_txt] = regexp(expDef_text, ...
                'quiescThreshold = (\d*)','match','tokens');
            quiescThreshold = str2num(quiescThreshold_txt{1}{1});
            % - quiscence step size
            [~,quiescenceStepSize_txt] = regexp(expDef_text, ...
                'quiescence_interval = (\d*)','match','tokens');
            if ~isempty(quiescenceStepSize_txt)
                quiescence_step_size = str2num(quiescenceStepSize_txt{1}{1});
            else
                % (if not included, it was fixed-quiescence)
                quiescence_step_size = [];
            end
            % - ITI step size
            [~,itiStepSize_txt] = regexp(expDef_text, ...
                'iti_interval = (\d*)','match','tokens');
            if ~isempty(itiStepSize_txt)
                iti_step_size = str2num(itiStepSize_txt{1}{1});
            else
                % (if not included, it was 1)
                iti_step_size = 1;
            end
            % - possible ITI times
            if isfield(block.paramsValues,itiMin)
                possible_iti = max([block.paramsValues.itiMin]):iti_step_size:max([block.paramsValues.itiMax]);
            else
                % (if not parameterized, was hard-coded)
                [~,itiRange_txt] = regexp(expDef_text, ...
                    'itiRange = \[(\d*):(\d*)\];','match','tokens');
                possible_iti = str2num(itiRange_txt{1}{1}):str2num(itiRange_txt{1}{2});
            end
            % -possible quiescence times
            if isfield(block.paramsValues,quiescenceMin) 
                possible_quiescence = max([block.paramsValues.quiescenceMin]):quiescence_step_size:max([block.paramsValues.quiescenceMax]);
            else
                % (if not parameterized, was fixed and hard-coded)
                [~,quiescenceTimeTxt_txt] = regexp(expDef_text, ...
                    'prestimQuiescentTime = (\d+(\.\d+)?)','match','tokens');
                possible_quiescence = str2num(quiescenceTimeTxt_txt{1}{1});
            end

            % ITI/quiescence time for each trial and parameters
            if isfield(block.events,'trialITIValues')
                % (for variable-quiescence operant: stored by trial)
                iti_t = block.events.trialITIValues;
            else
                % (for fixed-quiescence operant: defined as end-reponse)
                iti_t = block.events.endTrialTimes - block.events.responseTimes(1:n_trials);
                % (set signals events - used in
                % AP_extrap_stimWheel_quiescence)
                signals_events.trialITIValues = iti_t;
                signals_events.trialQuiescenceValues = quiescence_t;                
            end

            % Wheel movements/biases
            % mm in clicks from +hw.DaqRotaryEncoder, lilrig encoder = 100
            wheel_click2mm = 0.4869;
            wheel_mm = sum(abs(diff(wheel_position)))*wheel_click2mm;

            left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
            right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
            wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
                (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));

            % Time between prior move stop and post-stim move start
            prestim_move_pause = ...
                wheel_starts(wheel_move_stim_idx(2:end)) - ...
                wheel_stops(wheel_move_stim_idx(2:end)-1);

            % Time between all movements
            all_move_pause = wheel_starts(2:end)-wheel_stops(1:end-1);

            % Non-stim movements per second
            wheel_iti_move_rate = arrayfun(@(x) ...
                sum(wheel_starts > signals_events.responseTimes(x) & ...
                wheel_starts < stimOn_times(x+1))./ ...
                (stimOn_times(x+1) - signals_events.responseTimes(x)), ...
                1:n_trials-1);


            wheel_iti_move_rate = arrayfun(@(x) ...
                sum(wheel_starts > signals_events.responseTimes(x) & ...
                wheel_starts < stimOn_times(x+1))./ ...
                (stimOn_times(x+1) - signals_events.responseTimes(x)), ...
                1:n_trials-1);

            wheel_iti_move_rate = arrayfun(@(x) ...
                sum(wheel_starts > signals_events.responseTimes(x) & ...
                wheel_starts < stimOn_times(x+1))./ ...
                (stimOn_times(x+1) - signals_events.responseTimes(x)), ...
                1:n_trials-1);

            nonresponse_move_rate = ...
                (length(wheel_starts)-length(wheel_move_response_idx))./ ...
                block.duration;

            %%%%%%%%% STIM RESPONSE VS. NULL
            
            t = Timeline.rawDAQTimestamps';
            
            % Get quiescence reset times
            %             % (real: from new trial onset)
            %             quiescence_reset_t = AP_find_stimWheel_quiescence;
            % (extrapolated: from response to response)
            quiescence_reset_t_extrap = AP_extrap_stimWheel_quiescence;
            
            alt_stimOn_times = cell(n_trials,1);
            alt_stimOn_trialparams = cell(n_trials,1);
            % (skip trial 1: no ITI and bad quiescence watch)
            for curr_trial = 2:n_trials
                
                % Pull out current trial times (last to next response)
                curr_trial_t_idx = t >= signals_events.responseTimes(curr_trial-1) & ...
                    t <= signals_events.responseTimes(curr_trial);
                curr_trial_t = t(curr_trial_t_idx);
                
                t_from_quiescence_reset_trialitis = nan(length(curr_trial_t),length(possible_iti));
                for curr_possible_iti = 1:length(possible_iti)
                    curr_possible_itiend = signals_events.responseTimes(curr_trial-1) + ...
                        possible_iti(curr_possible_iti);
                    curr_quiescence_resets = sort([quiescence_reset_t_extrap;curr_possible_itiend]);
                    
                    t_from_quiescence_reset_trialitis(:,curr_possible_iti) = ...
                        curr_trial_t - interp1(curr_quiescence_resets, ...
                        curr_quiescence_resets,curr_trial_t,'previous','extrap');
                end
                
                % Find alternate stim times which would have given same first move
                
                % (getting possible iti + quiescence crosses)
                alt_iti_reached = ((t(curr_trial_t_idx) - curr_trial_t(1)) > possible_iti);
                alt_quiescence_reached = ...
                    t_from_quiescence_reset_trialitis > permute(possible_quiescence,[1,3,2]);
                
                % (get possible stim times as iti x quiescence grid)
                [alt_stim_value,alt_stim_idx] = max( ...
                    permute(alt_iti_reached & alt_quiescence_reached,[2,3,1]),[],3);
                alt_stim_t = curr_trial_t(alt_stim_idx);
                alt_stimOn_times_all = alt_stim_t(alt_stim_value);
                
                % (get alt stim times that would have resulted in the same
                % first movement since that's the measured value)
                stim_leeway = 0.1;
                curr_wheel_move_alt_stim_idx = ...
                    arrayfun(@(stim) find(wheel_starts > stim-stim_leeway,1,'first'), ...
                    alt_stimOn_times_all);
                use_alt_stimOn_times = ...
                    curr_wheel_move_alt_stim_idx == wheel_move_stim_idx(curr_trial);
                
                % (make sure that real parameters = real stim time: missed
                % wheel clicks sometimes give non-reproducible traces, in
                % which case the trial shouldn't be used for stats)
                curr_quiescence_idx = find(possible_quiescence == quiescence_t(curr_trial));
                curr_iti_idx = find(possible_iti == iti_t(curr_trial-1));
                curr_block_stimOn = signals_events.stimOnTimes(curr_trial);
                curr_alt_stim_offset = curr_block_stimOn - ...
                    alt_stim_t(curr_iti_idx,curr_quiescence_idx);
                if curr_alt_stim_offset > 0.01
                    continue
                end
                
                % (apply the block vs actual stim on time delay to all
                % times - note this is regularly backwards in time??)
                curr_stim_pd_offset = stimOn_times(curr_trial) - curr_block_stimOn;
                alt_stimOn_times_all_pd = alt_stimOn_times_all + curr_stim_pd_offset;
                
                % (store alternate stim times)
                alt_stimOn_times{curr_trial} = alt_stimOn_times_all_pd(use_alt_stimOn_times);
                
                % (trial plot)
                figure; hold on;
                t_plot_scale = 0.1;
                plot(t(curr_trial_t_idx),wheel_velocity(curr_trial_t_idx),'k')
                plot(t(curr_trial_t_idx),[0;diff(wheel_position(curr_trial_t_idx))]*0.1,'g')
                line(repmat(curr_trial_t(1)+signals_events.trialITIValues(curr_trial-1),2,1),ylim);
                line(xlim,repmat(signals_events.trialQuiescenceValues(curr_trial),2,1)*t_plot_scale,'color','m');
                line(repmat(curr_block_stimOn,1,2),ylim,'color','r','linestyle','--');
                line(repmat(stimOn_times(curr_trial),1,2),ylim,'color','k','linestyle','--');
                
                plot(alt_stimOn_times_all,0,'ob');
                plot(alt_stimOn_times{curr_trial},0,'.r');

                drawnow;
                
            end
            
            % Get would-be reaction time after alt stim times
            stim_leeway = 0.1;
            wheel_move_alt_stim_idx = ...
                arrayfun(@(stim) find(wheel_starts > stim-stim_leeway,1,'first'), ...
                cell2mat(alt_stimOn_times));
            
            alt_stim_to_move = ...
                mat2cell(wheel_starts(wheel_move_alt_stim_idx) - cell2mat(alt_stimOn_times), ...
                cellfun(@length,alt_stimOn_times));
            
            % Compare same stim/alt-stim trials (sub-sample alt-stim)
            use_alt_trials = cellfun(@(x) ~isempty(x),alt_stimOn_times);
            surround_t_centers = 0:0.01:2;
            
            % (day total)
            curr_trials = use_alt_trials;
            curr_stim_surround_times = stimOn_times(curr_trials) + surround_t_centers;
            curr_stim_surround_move = interp1(t,+wheel_move,curr_stim_surround_times,'previous');
            
            n_alt = 1000; % random subset alt trials for same n
            curr_alt_subset = cell2mat(cellfun(@(x) randsample(x,n_alt,true)', ...
                alt_stimOn_times(curr_trials),'uni',false));
            
            curr_alt_stim_surround_times = permute(curr_alt_subset,[1,3,2]) + surround_t_centers;
            curr_alt_stim_surround_move = interp1(t,+wheel_move,curr_alt_stim_surround_times,'previous');
            
            stim_move_max = max(nanmean(curr_stim_surround_move,1));
            alt_stim_move_max = nanmean(max(nanmean(curr_alt_stim_surround_move,1),[],2));
            
            stim_response_idx = (stim_move_max-alt_stim_move_max)./ ...
                (stim_move_max+alt_stim_move_max);
            
            % (daysplit)
            n_daysplit = 4;
            day_split_idx = min(floor(linspace(1,n_daysplit+1,n_trials)),n_daysplit)';
            
            stim_move_max_daysplit = nan(1,n_daysplit);
            alt_stim_move_max_daysplit = nan(1,n_daysplit);
            for curr_daysplit = 1:n_daysplit
                
                curr_trials = day_split_idx == curr_daysplit & use_alt_trials;
                
                curr_stim_surround_times = stimOn_times(curr_trials) + surround_t_centers;
                curr_stim_surround_move = interp1(t,+wheel_move,curr_stim_surround_times,'previous');
                
                n_alt = 1000; % random subset alt trials for same n
                curr_alt_subset = cell2mat(cellfun(@(x) randsample(x,n_alt,true)', ...
                    alt_stimOn_times(curr_trials),'uni',false));
                
                curr_alt_stim_surround_times = permute(curr_alt_subset,[1,3,2]) + surround_t_centers;
                curr_alt_stim_surround_move = interp1(t,+wheel_move,curr_alt_stim_surround_times,'previous');
                
                stim_move_max_daysplit(curr_daysplit) = max(nanmean(curr_stim_surround_move,1));
                alt_stim_move_max_daysplit(curr_daysplit) = nanmean(max(nanmean(curr_alt_stim_surround_move,1),[],2));
                
            end
            
            stim_response_idx_daysplit = (stim_move_max_daysplit-alt_stim_move_max_daysplit)./ ...
                (stim_move_max_daysplit+alt_stim_move_max_daysplit);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Store in behavior structure
            % (general timings)
            bhv(curr_animal).animal = animal;
            bhv(curr_animal).day{curr_day,1} = day;
            bhv(curr_animal).protocol{curr_day,1} = curr_protocol;
            bhv(curr_animal).session_duration(curr_day,1) = session_duration;
            bhv(curr_animal).n_trials(curr_day,1) = n_trials;
            bhv(curr_animal).frac_hit(curr_day,1) = frac_hit;
            bhv(curr_animal).total_water(curr_day,1) = total_water;
            bhv(curr_animal).wheel_mm(curr_day,1) = wheel_mm;
            bhv(curr_animal).wheel_velocity(curr_day,1) = nansum(abs(wheel_velocity));
            bhv(curr_animal).stim_move_t{curr_day,1} = stim_to_move;
            bhv(curr_animal).stim_feedback_t{curr_day,1} = stim_to_feedback;
            bhv(curr_animal).quiescence_t{curr_day,1} = quiescence_t';
            bhv(curr_animal).iti_t{curr_day,1} = iti_t';
            bhv(curr_animal).wheel_bias(curr_day,1) = wheel_bias;

            bhv(curr_animal).prestim_move_pause{curr_day,1} = prestim_move_pause;
            bhv(curr_animal).all_move_pause{curr_day,1} = all_move_pause;
            bhv(curr_animal).wheel_iti_move_rate{curr_day,1} = wheel_iti_move_rate;
            bhv(curr_animal).nonresponse_move_rate(curr_day,1) = nonresponse_move_rate;

            % (stim-aligned movement)
            bhv(curr_animal).stim_surround_t = stim_surround_t_centers;
            bhv(curr_animal).stim_surround_wheel{curr_day,1} = stim_surround_move;
            
            % (stim vs null response)
            bhv(curr_animal).alt_stim_move_t{curr_day,1} = alt_stim_to_move;
            bhv(curr_animal).alt_stim_trialparams{curr_day,1} = alt_stimOn_trialparams;
            
            bhv(curr_animal).stim_move_max{curr_day,1} = stim_move_max;
            bhv(curr_animal).alt_stim_move_max{curr_day,1} = alt_stim_move_max;
            bhv(curr_animal).stim_response_idx(curr_day,1) = stim_response_idx;
            
            bhv(curr_animal).stim_move_max_daysplit{curr_day,1} = stim_move_max_daysplit;
            bhv(curr_animal).alt_stim_move_max_daysplit{curr_day,1} = alt_stim_move_max_daysplit;
            bhv(curr_animal).stim_response_idx_daysplit{curr_day,1} = stim_response_idx_daysplit;
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
        
    end
    
    if plot_flag
        
        % Plot summary
        day_num = cellfun(@(x) datenum(x),{experiments.day}');
        day_labels = cellfun(@(day,protocol) [day(6:end)], ...
            {experiments.day}',bhv(curr_animal).protocol,'uni',false);
        
        [unique_protocols,~,protocol_idx] = unique(bhv(curr_animal).protocol);
        protocol_col = hsv(length(unique_protocols));
        
        figure('Name',animal)
        
        % Trials and water
        subplot(2,3,1); hold on;
        yyaxis left
        % plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'linewidth',2);
        % ylabel('Trials/min');
        plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'linewidth',2);
        ylabel('Trials/min');
        yyaxis right
        plot(day_num,bhv(curr_animal).total_water,'linewidth',2);
        ylabel('Total water');
        xlabel('Session');
        set(gca,'XTick',day_num);
        set(gca,'XTickLabel',day_labels);
        set(gca,'XTickLabelRotation',90);
        
        protocol_plot = gscatter(day_num,zeros(size(day_num)),[bhv(curr_animal).protocol]);
        legend off;
        
        imaging_days = day_num([experiments.imaging]);
        for i = 1:length(imaging_days)
            line(repmat(imaging_days(i),1,2),ylim,'color','k');
        end
        
        ephys_days = day_num([experiments.ephys]);
        for i = 1:length(ephys_days)
            line(repmat(ephys_days(i),1,2),ylim,'color','r','linestyle','--');
        end
        
        % Wheel movement and bias
        subplot(2,3,2);
        yyaxis left
        % plot(day_num,bhv(curr_animal).wheel_velocity./bhv(curr_animal).session_duration,'linewidth',2);
        % ylabel('Wheel movement / min');
        plot(day_num,bhv(curr_animal).wheel_velocity,'linewidth',2);
        ylabel('Wheel movement');
        yyaxis right
        plot(day_num,bhv(curr_animal).wheel_bias,'linewidth',2);
        ylim([-1,1]);
        line(xlim,[0,0]);
        ylabel('Wheel bias');
        xlabel('Session');
        set(gca,'XTick',day_num);
        set(gca,'XTickLabel',day_labels);
        set(gca,'XTickLabelRotation',90);
        
        imaging_days = day_num([experiments.imaging]);
        for i = 1:length(imaging_days)
            line(repmat(imaging_days(i),1,2),ylim,'color','k');
        end
        
        ephys_days = day_num([experiments.ephys]);
        for i = 1:length(ephys_days)
            line(repmat(ephys_days(i),1,2),ylim,'color','r','linestyle','--');
        end
        
        % Stim-to-reward time
        subplot(2,3,3);
        plot(day_num,cellfun(@nanmedian,bhv(curr_animal).stim_feedback_t),'k','linewidth',2);
        ylabel('Stim to reward time (s)');
        xlabel('Session')
        set(gca,'XTick',day_num);
        set(gca,'XTickLabel',day_labels);
        set(gca,'XTickLabelRotation',90);
        
        % Move fraction heatmap
        stim_surround_wheel_cat = ...
            cell2mat(cellfun(@(x) nanmean(x,1),bhv(curr_animal).stim_surround_wheel,'uni',false));
        subplot(2,3,4);
        imagesc(bhv(curr_animal).stim_surround_t,[],stim_surround_wheel_cat);
        colormap(gca,brewermap([],'Greys'));
        caxis([0,1]);
        xlabel('Time from stim (s)');
        set(gca,'YTick',1:length(day_num));
        set(gca,'YTickLabel',day_labels);
        
        % Move fraction lineplot
        subplot(2,3,5); hold on;
        set(gca,'ColorOrder',copper(size(stim_surround_wheel_cat,1)));
        plot(bhv(curr_animal).stim_surround_t,stim_surround_wheel_cat');
        ylim([0,1]);
        ylabel('Fraction move');
        xlabel('Time from stim (s)');
        
        % Move fraction pre/post stim
        stim_surround_t = bhv(curr_animal).stim_surround_t;
        move_prestim_max = ...
            cell2mat(cellfun(@(x) max(nanmean(x(:,stim_surround_t < 0),1)), ...
            bhv(curr_animal).stim_surround_wheel,'uni',false));
        move_poststim_max = ...
            cell2mat(cellfun(@(x) max(nanmean(x(:,stim_surround_t > 0),1)), ...
            bhv(curr_animal).stim_surround_wheel,'uni',false));
        move_prepost_max_ratio = ...
            (move_poststim_max-move_prestim_max)./(move_poststim_max+move_prestim_max);
        
        subplot(2,3,6); hold on;
        plot(day_num,move_prepost_max_ratio,'k','linewidth',2);
        ylabel({'Post:pre max move ratio'});
        ylim([-1,1]);
        line(xlim,[0,0],'color','k','linestyle','--');
        set(gca,'XTick',day_num);
        set(gca,'XTickLabel',day_labels);
        set(gca,'XTickLabelRotation',90);
        
        drawnow;
    end
    
    clearvars('-except',preload_vars{:});
    
end

if save_bhv
    data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\cstr_data_test';
    save_fn = [data_path filesep 'bhv_' animal_group];
    save(save_fn,'bhv');
    disp(['Saved ' save_fn]);
end

%% Find "learned" days
% Days when response to stimulus was significantly different from null

% Load behavior from above
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\cstr_data_test';
% bhv_fn = [data_path filesep 'bhv_cstr'];
bhv_fn = [data_path filesep 'bhv_cstr_fixedq'];
load(bhv_fn);

% Define "learned" days by median reaction time
% (exclude rxn < 0.1: too fast for stim, doesn't change with learning)
% (keeping window defined, but not used for stats)
n_sample = 10000;
learned_days = cell(size(bhv));
for curr_animal = 1:length(bhv)

    curr_n_days = length(bhv(curr_animal).stim_move_t);
    learned_days{curr_animal} = nan(curr_n_days,1);

    for curr_day = 1:curr_n_days

        % (use trials only with alt - no alt if couldn't recapture stim
        % time with params, happens with some skipped wheel clicks)
        use_trials = ~cellfun(@isempty,bhv(curr_animal).alt_stim_move_t{curr_day});
        
        curr_rxn = bhv(curr_animal).stim_move_t{curr_day}(use_trials);
        curr_alt_rxn = cell2mat(cellfun(@(x) datasample(x,n_sample)', ...
            bhv(curr_animal).alt_stim_move_t{curr_day}(use_trials),'uni',false));     
        
        % Get measured and null medians (exclude rxn < 0.1)
%         rxn_stat = nanmedian(curr_rxn.*AP_nanout(curr_rxn < 0.1),1);
%         alt_rxn_stat = nanmedian(curr_alt_rxn.*AP_nanout(curr_alt_rxn < 0.1),1);

        rxn_stat = mad(curr_rxn.*AP_nanout(curr_rxn < 0),1);
        alt_rxn_stat = mad(curr_alt_rxn.*AP_nanout(curr_alt_rxn < 0),1);
        
        rxn_stat_rank = tiedrank([rxn_stat,alt_rxn_stat]);
        rxn_stat_p = rxn_stat_rank(1)./(n_sample+1);
        
        % (null rejected at 5%)
        rxn_window_frac_h = rxn_stat_p < 0.05;

        learned_days{curr_animal}(curr_day) = rxn_window_frac_h;
    end
end

% Plot to check
learned_days_padcat = AP_padcatcell(learned_days);

figure; 
imagesc(learned_days_padcat,'AlphaData',~isnan(learned_days_padcat));
set(gca,'Color',[0.5,0.5,0.5])
title('Sig. reaction time stat days');
set(gca,'XTick',1:length({bhv.animal}),'XTickLabel',{bhv.animal})
xlabel('Animal');
ylabel('Day');

% Put learned days into behavior structure and save
[bhv.learned_days] = learned_days{:};
save(bhv_fn,'bhv');
disp(['Saved learned days ' bhv_fn]);



%% Passive widefield

disp('Passive trial activity')

% Initialize
clear all
trial_data_all = struct;

% animals = {'AP092','AP094','AP096'}; % learned operant
% animals = {'AP093','AP095','AP097'}; % didn't learn operant?
% animals = {'AP089','AP090','AP091'}; % fixed quiescence 

% (combined)
animals = {'AP092','AP094','AP096','AP093','AP095','AP097','AP089','AP090','AP091'};

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    passive_experiments = AP_find_experiments_externalhd(animal,protocol);
    
    % Get days after muscimol starts (to exclude)
    data_path = 'E:\Cortexlab_documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    if any(muscimol_animal_idx)
        muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
        muscimol_experiments = datenum({passive_experiments.day})' >= datenum(muscimol_start_day);
    else
        muscimol_experiments = false(size({passive_experiments.day}));
    end

    % Get days with standard task (to exclude other tasks)
    task_protocol = 'AP_stimWheelRight';
    task_experiments = AP_find_experiments_externalhd(animal,task_protocol);
    standard_task_experiments =  datenum({passive_experiments.day})' <= ...
        datenum(task_experiments(end).day);
    
    % Set experiments to use (imaging, not muscimol, standard task)
    experiments = passive_experiments([passive_experiments.imaging] & ...
        ~muscimol_experiments & standard_task_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment_externalhd;
        
        % Pull out trial data
        operant_grab_trial_data_externalhd;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\cstr_data_test';
save_fn = fullfile(save_path,'trial_activity_passive_cstr');
save(save_fn,'-v7.3');
disp(['Saved: ' save_fn])


%% Load passive data

% Load data
trial_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\cstr_data_test';
data_fn = 'trial_activity_passive_cstr';
AP_load_trials_operant_externalhd;
n_naive = 3; % (number of naive passive-only days, just hard-coding)

% Load behavior, exclude animals not in dataset, get learned day
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\cstr_data_test';
bhv_fn = [data_path filesep 'bhv_cstr'];
load(bhv_fn);
bhv = bhv(ismember({bhv.animal},animals)); 
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days})';

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Get "learned day" for each trial
trial_learned_day = cell2mat(cellfun(@(x,ld) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell([1:length(x)]-ld)',x,'uni',false)), ...
    wheel_all,num2cell(learned_day + n_naive),'uni',false));

% Define learning stages
% (learning stages: day 1-3 passive-only, then pre/post learning)
trial_stage = 1*(trial_day <= 3) + ...
    2*(trial_day > 3 & trial_learned_day < 0) + ...
    3*(trial_learned_day >= 0);
n_stages = length(unique(trial_stage));

% Get trials with movement during stim to exclude
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);

% (print total fraction of quiescent trials)
q_exp = mat2cell(quiescent_trials,cellfun(@(x) size(x,1),vertcat(wheel_all{:})),1);
s_exp = mat2cell(trial_stim_allcat,cellfun(@(x) size(x,1),vertcat(wheel_all{:})),1);
q_exp_frac = cell2mat(cellfun(@(q,s) grpstats(q,s),q_exp,s_exp,'uni',false)');
fprintf('Quiescent trial frac: \n%.2f+-%.2f (left),%.2f+-%.2f (center),%.2f+-%.2f (right)\n', ...
    reshape([nanmean(q_exp_frac,2),nanstd(q_exp_frac,[],2)]',[],1))

% Turn values into IDs for grouping
stim_unique = unique(trial_stim_allcat);
[~,trial_stim_id] = ismember(trial_stim_allcat,stim_unique);

learned_day_unique = unique(trial_learned_day);
[~,trial_learned_day_id] = ismember(trial_learned_day,learned_day_unique);

% Get average fluorescence by animal/day/stim
stim_v_avg = cell(length(animals),1);
stim_roi_avg = cell(length(animals),1);
for curr_animal = 1:length(animals)
    for curr_day = 1:max(trial_day(trial_animal == curr_animal))
        for curr_stim_idx = 1:length(stim_unique)

            use_trials = quiescent_trials & ...
                trial_animal == curr_animal & ...
                trial_day == curr_day & ...
                trial_stim_allcat == stim_unique(curr_stim_idx);

            stim_v_avg{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean(fluor_allcat_deconv(use_trials,:,:),1),[3,2,1]);

            stim_roi_avg{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean(fluor_roi_deconv(use_trials,:,:),1),[3,2,1]);
            
        end
    end
end

% Plot average across animals
use_trials = quiescent_trials & trial_stim_allcat == 1;
[trial_groups,~,trial_group_idx] = unique([trial_animal(use_trials),trial_learned_day(use_trials)],'rows');
fluor_trialavg = ap.groupfun(@nanmean,fluor_allcat_deconv(use_trials,:,:),trial_group_idx,[],[]);
fluor_animalavg = ap.groupfun(@nanmean,fluor_trialavg,trial_groups(:,2),[],[]);

ld_x = unique(trial_groups(:,2));

fluor_animalavg_px = plab.wf.svd2px(U_master(:,:,1:n_vs),permute(fluor_animalavg,[3,2,1]));

ap.imscroll(fluor_animalavg_px(:,:,:,ismember(ld_x,-3:3)))
axis image;


% (manual - avg and plot)

use_trials = trial_day < 5 & trial_stim_allcat == 1 & quiescent_trials;
% use_trials = trial_day > 8 & trial_stim_allcat == 1 & quiescent_trials;

px = plab.wf.svd2px(U_master(:,:,1:n_vs),permute(nanmean(fluor_allcat_deconv(use_trials,:,:),1),[3,2,1]));
ap.imscroll(px,t); 
axis image;
colormap(AP_colormap('PWG',[],1.5));
clim(max(abs(clim)).*[-1,1]*0.8);l

use_trials = trial_day < 5 & trial_stim_allcat == 0 & quiescent_trials;
% use_trials = trial_day > 8 & trial_stim_allcat == 0 & quiescent_trials;

px = plab.wf.svd2px(U_master(:,:,1:n_vs),permute(nanmean(fluor_allcat_deconv(use_trials,:,:),1),[3,2,1]));
ap.imscroll(px,t); 
axis image;
colormap(AP_colormap('PWG',[],1.5));
clim(max(abs(clim)).*[-1,1]*0.8);l





