
%% Plot behavior of alt tasks

% animals = ...
%     {'AM001','AM002','AM004','AM005','AM006','AM007','AM008','AM009','AM011','AM012', ...
%     'AM013','AM014','AM015','AM016','AM017','AM018','AM019','AM021','AM022','AM023', ...
%     'AM024','AM025','AM026','AM027','AM029','AM030','AP004','AP005','AP006','AP007', ...
%     'AP008','AP009','AP010','AP011','AP012','AP013','AP014','AP015','AP016','AP017', ...
%     'AP018','AP020','AP021','AP022','AP023','AP025','DS001','DS007','DS009','DS010', ...
%     'DS011'};

animals = {'AP019', ... % no change
    'AP027','AP028','AP029', ... % opacity
    'HA000','HA001','HA002'}; % angle

% Set reaction statistic to use
use_stat = 'median';

% Grab reaction stats and learning day for each mouse
bhv = struct;

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};
    fprintf('%s (%d/%d)\n',animal,curr_animal_idx,length(animals))

    use_workflow = 'stim_wheel_right_stage\d';

    use_workflow = {'*no_change*','*opacity*','*angle'};

    recordings = plab.find_recordings(animal,[],use_workflow);
    relative_day = days(datetime({recordings.day}) - datetime({recordings(1).day}))+1;
    bhv(curr_animal_idx).relative_day = relative_day;

    rxn_stat_p = nan(length(recordings),1);
    rxn_stat = nan(length(recordings),1);
    rxn_null_stat = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get association stat
        % (skip if only a few trials)
        if n_trials < 10
            continue
        end
       
        [rxn_stat_p(curr_recording), ...
            rxn_stat(curr_recording),rxn_null_stat(curr_recording)] = ...
            AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});
        AP_print_progress_fraction(curr_recording,length(recordings));

    end

    % Define learned day from reaction stat p-value and reaction time
    learned_day = rxn_stat_p < 0.05;

    % Store behavior across animals
    bhv(curr_animal_idx).rxn_stat = rxn_stat;
    bhv(curr_animal_idx).rxn_null_stat = rxn_null_stat;
    bhv(curr_animal_idx).learned_day = learned_day;

end


% Identify animals that had other training before task
workflow = 'stim_wheel_right_stage1';
alt_task_animal = false(size(animals));
passive_workflows = {'lcr_passive','hml_passive','sparse_noise'};
for curr_animal = 1:length(animals)
    recordings = plab.find_recordings(animals{curr_animal});
    workflow_cat = vertcat(recordings.workflow);
    workflow_cat_sanitized = workflow_cat(~cellfun(@isempty,workflow_cat));

    alt_task_animal(curr_animal) = ...
        all(ismember(workflow_cat_sanitized(1: ...
        find(strcmp(workflow_cat_sanitized,workflow),1)-1), ...
        passive_workflows));
end

learned_animals = find(cellfun(@any,{bhv.learned_day}));
learned_day = cellfun(@(x) nanmax([NaN,find(x,1)]),{bhv.learned_day});
relative_learned_day = cell2mat(cellfun(@(x,y) x-y,{bhv.relative_day}, ...
    num2cell(learned_day),'uni',false))';

max_plot_day = 7; 

% Plot days (only with <2d break)
figure; hold on;
arrayfun(@(x) plot(bhv(x).relative_day(1:min([length(bhv(x).relative_day),find(diff(bhv(x).relative_day)>2)-1])), ...
    bhv(x).rxn_stat(1:min([length(bhv(x).relative_day),find(diff(bhv(x).relative_day)>2)-1]))),learned_animals);

[rxn_stat_mean,rxn_stat_day] = ap.groupfun(@nanmean,vertcat(bhv.rxn_stat),horzcat(bhv.relative_day));
plot_idx = rxn_stat_day <= max_plot_day;
[rxn_stat_sem,rxn_stat_day] = ap.groupfun(@AP_sem,vertcat(bhv.rxn_stat),horzcat(bhv.relative_day));
errorbar(rxn_stat_day(plot_idx),rxn_stat_mean(plot_idx), ...
    rxn_stat_sem(plot_idx),'k','linewidth',2)

% Plot only group measured/null average

[rxn_stat_ld_mean,ld_group] = ap.groupfun(@nanmean,vertcat(bhv.rxn_stat),relative_learned_day);
rxn_stat_ld_sem = ap.groupfun(@AP_sem,vertcat(bhv.rxn_stat),relative_learned_day);
plot_ld_idx = ld_group >= -max_plot_day/2 & ld_group <= max_plot_day/2;

rxn_null_stat_ld_mean = ap.groupfun(@nanmean,vertcat(bhv.rxn_null_stat),relative_learned_day);
rxn_null_stat_ld_sem = ap.groupfun(@AP_sem,vertcat(bhv.rxn_null_stat),relative_learned_day);

rxn_stat_null_mean = ap.groupfun(@nanmean,vertcat(bhv.rxn_null_stat),horzcat(bhv.relative_day));
rxn_stat_null_sem = ap.groupfun(@AP_sem,vertcat(bhv.rxn_null_stat),horzcat(bhv.relative_day));

figure; h = tiledlayout(2,2);

nexttile;
errorbar(rxn_stat_day(plot_idx),rxn_stat_mean(plot_idx), ...
    rxn_stat_sem(plot_idx),'k','linewidth',2)
xlim([0,max_plot_day+1]);

nexttile;
errorbar(rxn_stat_day(plot_idx),rxn_stat_null_mean(plot_idx), ...
    rxn_stat_null_sem(plot_idx),'r','linewidth',2)
xlim([0,max_plot_day+1]);

nexttile;
errorbar(ld_group(plot_ld_idx),rxn_stat_ld_mean(plot_ld_idx), ...
    rxn_stat_ld_sem(plot_ld_idx),'k','linewidth',2)

nexttile;
errorbar(ld_group(plot_ld_idx),rxn_null_stat_ld_mean(plot_ld_idx), ...
    rxn_null_stat_ld_sem(plot_ld_idx),'r','linewidth',2)

linkaxes(h.Children,'y');

figure;
h = tiledlayout(1,2);

nexttile;
rxn_idx = (vertcat(bhv.rxn_null_stat)-vertcat(bhv.rxn_stat))./ ...
    (vertcat(bhv.rxn_null_stat)+vertcat(bhv.rxn_stat));
rxn_idx_mean = ap.groupfun(@nanmean,rxn_idx,horzcat(bhv.relative_day));
rxn_idx_sem = ap.groupfun(@AP_sem,rxn_idx,horzcat(bhv.relative_day));
errorbar(rxn_stat_day(plot_idx),rxn_idx_mean(plot_idx), ...
    rxn_idx_sem(plot_idx),'k','linewidth',2)
ylabel('Reaction index')

nexttile;
[rxn_idx_mean,ld_group] = ap.groupfun(@nanmean,rxn_idx,relative_learned_day);
rxn_idx_sem = ap.groupfun(@AP_sem,rxn_idx,relative_learned_day);
errorbar(ld_group(plot_ld_idx),rxn_idx_mean(plot_ld_idx), ...
    rxn_idx_sem(plot_ld_idx),'k','linewidth',2)
ylabel('Reaction index')


figure;h = histogram(learned_day);
h.FaceAlpha = 1;
h.FaceColor = [0.5,0.5,0.5];




%% Plot behavior of av animals

animals = {'DS000','DS004','DS014','DS015','DS016'};

% Set reaction statistic to use
use_stat = 'median';

% Grab reaction stats and learning day for each mouse
bhv = struct;

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};
    fprintf('%s (%d/%d)\n',animal,curr_animal_idx,length(animals))

%     use_workflow = 'stim_wheel_right_stage\d';
    use_workflow = '*audio_volume*';

    recordings = plab.find_recordings(animal,[],use_workflow);
    relative_day = days(datetime({recordings.day}) - datetime({recordings(1).day}))+1;
    bhv(curr_animal_idx).relative_day = relative_day;

    rxn_stat_p = nan(length(recordings),1);
    rxn_stat = nan(length(recordings),1);
    rxn_null_stat = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get association stat
        % (skip if only a few trials)
        if n_trials < 10
            continue
        end
       
        [rxn_stat_p(curr_recording), ...
            rxn_stat(curr_recording),rxn_null_stat(curr_recording)] = ...
            AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});
        AP_print_progress_fraction(curr_recording,length(recordings));

    end

    % Define learned day from reaction stat p-value and reaction time
    learned_day = rxn_stat_p < 0.05;

    % Store behavior across animals
    bhv(curr_animal_idx).rxn_stat = rxn_stat;
    bhv(curr_animal_idx).rxn_null_stat = rxn_null_stat;
    bhv(curr_animal_idx).learned_day = learned_day;

end


learned_animals = find(cellfun(@any,{bhv.learned_day}));
learned_day = cellfun(@(x) nanmax([NaN,find(x,1)]),{bhv.learned_day});
relative_learned_day = cell2mat(cellfun(@(x,y) x-y,{bhv.relative_day}, ...
    num2cell(learned_day),'uni',false))';

max_plot_day = 7; 

% Plot days (only with <2d break)
figure; hold on;
arrayfun(@(x) plot(bhv(x).relative_day(1:min([length(bhv(x).relative_day),find(diff(bhv(x).relative_day)>2)-1])), ...
    bhv(x).rxn_stat(1:min([length(bhv(x).relative_day),find(diff(bhv(x).relative_day)>2)-1]))),learned_animals);

[rxn_stat_mean,rxn_stat_day] = ap.groupfun(@nanmean,vertcat(bhv.rxn_stat),horzcat(bhv.relative_day));
plot_idx = rxn_stat_day <= max_plot_day;
[rxn_stat_sem,rxn_stat_day] = ap.groupfun(@AP_sem,vertcat(bhv.rxn_stat),horzcat(bhv.relative_day));
errorbar(rxn_stat_day(plot_idx),rxn_stat_mean(plot_idx), ...
    rxn_stat_sem(plot_idx),'k','linewidth',2)

% Plot only group measured/null average

[rxn_stat_ld_mean,ld_group] = ap.groupfun(@nanmean,vertcat(bhv.rxn_stat),relative_learned_day);
rxn_stat_ld_sem = ap.groupfun(@AP_sem,vertcat(bhv.rxn_stat),relative_learned_day);
plot_ld_idx = ld_group >= -max_plot_day/2 & ld_group <= max_plot_day/2;

rxn_null_stat_ld_mean = ap.groupfun(@nanmean,vertcat(bhv.rxn_null_stat),relative_learned_day);
rxn_null_stat_ld_sem = ap.groupfun(@AP_sem,vertcat(bhv.rxn_null_stat),relative_learned_day);

rxn_stat_null_mean = ap.groupfun(@nanmean,vertcat(bhv.rxn_null_stat),horzcat(bhv.relative_day));
rxn_stat_null_sem = ap.groupfun(@AP_sem,vertcat(bhv.rxn_null_stat),horzcat(bhv.relative_day));

figure; h = tiledlayout(2,2);

nexttile;
errorbar(rxn_stat_day(plot_idx),rxn_stat_mean(plot_idx), ...
    rxn_stat_sem(plot_idx),'k','linewidth',2)
xlim([0,max_plot_day+1]);

nexttile;
errorbar(rxn_stat_day(plot_idx),rxn_stat_null_mean(plot_idx), ...
    rxn_stat_null_sem(plot_idx),'r','linewidth',2)
xlim([0,max_plot_day+1]);

nexttile;
errorbar(ld_group(plot_ld_idx),rxn_stat_ld_mean(plot_ld_idx), ...
    rxn_stat_ld_sem(plot_ld_idx),'k','linewidth',2)

nexttile;
errorbar(ld_group(plot_ld_idx),rxn_null_stat_ld_mean(plot_ld_idx), ...
    rxn_null_stat_ld_sem(plot_ld_idx),'r','linewidth',2)

linkaxes(h.Children,'y');

figure;
h = tiledlayout(1,2);

nexttile;
rxn_idx = (vertcat(bhv.rxn_null_stat)-vertcat(bhv.rxn_stat))./ ...
    (vertcat(bhv.rxn_null_stat)+vertcat(bhv.rxn_stat));
rxn_idx_mean = ap.groupfun(@nanmean,rxn_idx,horzcat(bhv.relative_day));
rxn_idx_sem = ap.groupfun(@AP_sem,rxn_idx,horzcat(bhv.relative_day));
errorbar(rxn_stat_day(plot_idx),rxn_idx_mean(plot_idx), ...
    rxn_idx_sem(plot_idx),'k','linewidth',2)
ylabel('Reaction index')

nexttile;
[rxn_idx_mean,ld_group] = ap.groupfun(@nanmean,rxn_idx,relative_learned_day);
rxn_idx_sem = ap.groupfun(@AP_sem,rxn_idx,relative_learned_day);
errorbar(ld_group(plot_ld_idx),rxn_idx_mean(plot_ld_idx), ...
    rxn_idx_sem(plot_ld_idx),'k','linewidth',2)
ylabel('Reaction index')


figure;h = histogram(learned_day);
h.FaceAlpha = 1;
h.FaceColor = [0.5,0.5,0.5];



