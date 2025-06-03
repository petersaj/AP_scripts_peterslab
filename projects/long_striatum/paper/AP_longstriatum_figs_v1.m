%% Figures for longitudinal striatum+cortex project
%
% Data packaged by AM save scripts
% Analysis taken fom `AP_longstriatum_exploratory_2`
%
% Uses `AP_longstriatum_load_data` to load/prep relevant data at the start
% of each figure (doesn't re-load if already overloaded - this lets each
% panel be run independently, or for the entire script to be run straight
% through)

tic

%% [Fig 1X] Reaction stat and learning histogram

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

% Plot reaction time/index, by training day and learning day
training_day = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));

rxn_cat = cell2mat(bhv.(['stimwheel_rxn_',use_stat]));
rxn_null_cat = cell2mat(bhv.(['stimwheel_rxn_null_',use_stat]));
rxn_idx_cat = (rxn_null_cat-rxn_cat)./(rxn_null_cat+rxn_cat);

figure; h = tiledlayout(2,2);

for day_align = 1:2

    switch day_align
        case 1
            use_plot_day_grp = training_day;
            plot_x = 1:7;
        case 2
            use_plot_day_grp = bhv.days_from_learning;
            plot_x = -3:2;
    end

    [rxn_idx_mean,rxn_group_x] = ap.groupfun(@mean,rxn_idx_cat,use_plot_day_grp);
    rxn_idx_sem = ap.groupfun(@AP_sem,rxn_idx_cat,use_plot_day_grp);

    rxn_mean = ap.groupfun(@mean,rxn_cat,use_plot_day_grp);
    rxn_sem = ap.groupfun(@AP_sem,rxn_cat,use_plot_day_grp);

    rxn_null_mean = ap.groupfun(@mean,rxn_null_cat,use_plot_day_grp);
    rxn_null_sem = ap.groupfun(@AP_sem,rxn_null_cat,use_plot_day_grp);

    plot_x_idx = ismember(rxn_group_x,plot_x);

    nexttile; hold on;
    % ap.errorfill(rxn_group_ld,rxn_null_mean,rxn_null_sem,'r',0.5,false);
    % errorbar(rxn_group_ld,rxn_mean,rxn_sem,'k','linewidth',2);
    plot(rxn_group_x(plot_x_idx),rxn_null_mean(plot_x_idx),'r','linewidth',2);
    plot(rxn_group_x(plot_x_idx),rxn_mean(plot_x_idx),'k','linewidth',2);
    set(gca,'YScale','log');
    xline(0);
    ylabel('Reaction index')

    nexttile
    errorbar(rxn_group_x(plot_x_idx),rxn_idx_mean(plot_x_idx),rxn_idx_sem(plot_x_idx),'k','linewidth',2);
    xline(0);
    ylabel('Reaction index')

end

linkaxes(h.Children(1:2:end),'y');
linkaxes(h.Children(2:2:end),'y');
ap.prettyfig;

% Plot histogram of learning days
n_learned_day = cellfun(@(x) max([0, ...
    find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'));
figure;
histogram(n_learned_day,-0.5:max(n_learned_day)+0.5,'FaceColor','k','EdgeColor','none');
ylabel('Number of mice');
xlabel('Days to learn');
ap.prettyfig;


%% [Fig 1X] K-means cluster means

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

figure;
h = tiledlayout(n_k,1,'tilespacing','none');
for curr_k = 1:n_k
    nexttile;
    imagesc(kmeans_cluster_mean(:,:,curr_k));
    axis image off
    clim(max(abs(clim)).*[-1,1]);
    colormap(ap.colormap('PWG'));
    ap.wf_draw('ccf',[0.5,0.5,0.5]);
end
ap.prettyfig;


%% [Fig 1X] Striatal MUA pre/post learning

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,0,Inf];
plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
learn_colormap = ap.colormap('BKR',3);
prepost_colormap = max(0,learn_colormap([1,end],:)-0.2);

% Plot average task activity (pre/post learning)
stim_x = [-0.2,0.3];
move_x = [-0.05,0.4];
outcome_x = [-0.1,0.5];

figure; h = tiledlayout(n_k,3,'TileSpacing','tight');
for curr_depth = 1:n_k

    curr_trials = striatum_mua_grp.kidx == curr_depth;

    % Stim
    nexttile; hold on; set(gca,'ColorOrder',prepost_colormap);
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean}, ...
        striatum_mua(curr_trials,:),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem}, ...
        striatum_mua(curr_trials,:),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(stim_x);
    xline(0);

    % Move
    nexttile; hold on; set(gca,'ColorOrder',prepost_colormap);
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean}, ...
        striatum_move_mua(curr_trials,:),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem}, ...
        striatum_move_mua(curr_trials,:),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(move_x);
    xline(0);

    % Outcome
    nexttile; hold on; set(gca,'ColorOrder',prepost_colormap);
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean}, ...
        striatum_outcome_mua(curr_trials,:),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem}, ...
        striatum_outcome_mua(curr_trials,:),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(outcome_x);
    xline(0);

end
% (link all y, and x of same-alignment)
linkaxes(h.Children,'y');
for ax = 1:3
    linkaxes(h.Children(ax:3:end),'x');
end
% (set same data aspect to have same x-size)
[h.Children.DataAspectRatio] = deal(min(vertcat(h.Children.DataAspectRatio),[],1));
ap.prettyfig;


%% [Fig 2X] Striatum task trial heatmap (reaction-sorted)

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

heatmap_smooth = [20,1]; % ([trials,time] to smooth for graphics)
figure; tiledlayout(n_k,max(plot_day_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(plot_day_grp)'
        nexttile;
        curr_trials = find(plot_day_grp == curr_day & striatum_mua_grp.kidx == curr_k);
       
        [~,sort_idx] = sort(striatum_mua_grp.rxn(curr_trials));      
        imagesc(psth_t,[],movmean(striatum_mua(curr_trials(sort_idx),:),heatmap_smooth));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;
        if curr_k == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end
ap.prettyfig;

%% [Fig 2X] Striatum task average PSTH w/ and w/o movement

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

% Plot average activity in trial (NaN-out movement times)
figure; h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','none');
move_leeway = 0.1; % time pre-movement to exclude
for curr_k = 1:n_k
    for curr_day = unique(plot_day_grp)'
        nexttile; hold on;
        curr_trials = find(plot_day_grp == curr_day & striatum_mua_grp.kidx == curr_k);      
        
        % Average all data
        curr_data = striatum_mua(curr_trials,:);
        curr_data_no_move = striatum_mua(curr_trials,:).* ...
            ap.nanout(psth_t > striatum_mua_grp.rxn(curr_trials)-move_leeway);

        curr_data_mean = mean(ap.groupfun(@nanmean,curr_data,striatum_mua_grp.animal(curr_trials)),1);
        curr_data_sem = AP_sem(ap.groupfun(@nanmean,curr_data,striatum_mua_grp.animal(curr_trials)),1);

        % Average data with movement removed
        % (only plot no-move timepoints with at least N trials and N animals)
        min_nomove_trials = 10;
        min_noanimals_trials = 5;

        curr_data_no_move_n = ap.groupfun(@(x) sum(~isnan(x)),curr_data_no_move,striatum_mua_grp.animal(curr_trials));
        curr_data_no_move_animal_mean = ap.groupfun(@nanmean,curr_data_no_move,striatum_mua_grp.animal(curr_trials));
        curr_data_no_move_animal_mean_mindata = curr_data_no_move_animal_mean.* ...
            ap.nanout(curr_data_no_move_n < min_nomove_trials).* ...
            ap.nanout(sum(curr_data_no_move_n >= min_nomove_trials) < min_noanimals_trials);
        curr_data_no_move_mean = nanmean(curr_data_no_move_animal_mean_mindata,1);
        curr_data_no_move_sem = AP_sem(curr_data_no_move_animal_mean_mindata,1);

        % Plot w/ movement (shaded) and w/o movement (line)
        ap.errorfill(psth_t,curr_data_mean,curr_data_sem,'k',0.3,false,2);
        plot(psth_t,curr_data_no_move_mean,'k','linewidth',2);
        axis off;
        if curr_k == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Fig 2X] Striatum task max stim activity split aross session

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

rxn_cutoff = 0.3; % only plot trials with slow reaction times

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

% Split by trial percentile within day
n_split = 3;
split_idx = cell2mat(cellfun(@(n_trials,n_mua) ...
    reshape(repmat(ap.quantile_bin(n_trials,n_split),1,n_mua),[],1), ...
    num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

% (mean across trials, max in time window, mean across animals)
use_trials = striatum_mua_grp.rxn > rxn_cutoff;
[activity_mean,activity_mean_grp] = ap.groupfun(@mean,striatum_mua(use_trials,:), ...
    [striatum_mua_grp.animal(use_trials),plot_day_grp(use_trials), ...
    split_idx(use_trials),striatum_mua_grp.kidx(use_trials)]);

max_t = psth_t > 0 & psth_t < 0.2;
activity_max = max(activity_mean(:,max_t),[],2);
[activity_max_avg,activity_max_grp] = ap.nestgroupfun({@mean,@mean},activity_max,activity_mean_grp(:,1),activity_mean_grp(:,2:end));
activity_max_sem = ap.nestgroupfun({@mean,@AP_sem},activity_max,activity_mean_grp(:,1),activity_mean_grp(:,2:end));

figure;
h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day = 1:length(plot_day_bins)-1
        nexttile; hold on;
        curr_data_idx = activity_max_grp(:,1) == curr_day & ...
            activity_max_grp(:,3) == curr_k;
        errorbar(activity_max_avg(curr_data_idx),activity_max_sem(curr_data_idx),'k','linewidth',2);
        axis padded
        if curr_k == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Fig 2X] Cortex ROIs

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

figure;tiledlayout(n_k,1,'tilespacing','none')
for curr_k = 1:n_k
    nexttile;
    imagesc(striatum_wf_roi(:,:,curr_k));
    axis image off; ap.wf_draw('ccf',[0.5,0.5,0.5]);
    colormap(ap.colormap('WG'));
end
ap.prettyfig;


%% [Fig 2X] Cortex task trial heatmap (reaction-sorted)

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% Plot heatmaps sorted by reaction times
heatmap_smooth = [20,1]; % ([trials,time] to smooth for graphics)
figure; tiledlayout(n_k,max(plot_day_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(plot_day_grp)'
        nexttile;
        curr_trials = find(plot_day_grp == curr_day);
       
        [~,sort_idx] = sort(wf_grp.rxn(curr_trials));
        imagesc(wf_t,[],movmean(wf_striatum_roi(curr_trials(sort_idx),:,curr_k),heatmap_smooth));
        colormap(ap.colormap('PWG'));
        clim(max(abs(clim)).*[-1,1]);
        axis off;
        if curr_k == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end
ap.prettyfig;


%% [Fig 2X] Cortex task average PSTH w/ and w/o movement

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% Plot average activity in trial (NaN-out movement times)
figure; h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','none');
move_leeway = 0.1;
for curr_k = 1:n_k
    for curr_day = unique(plot_day_grp)'
        nexttile; hold on;
        curr_trials = find(plot_day_grp == curr_day);      
        
        % Average all data
        curr_data = wf_striatum_roi(curr_trials,:,curr_k);

        curr_data_mean = mean(ap.groupfun(@nanmean,curr_data,wf_grp.animal(curr_trials)),1);
        curr_data_sem = AP_sem(ap.groupfun(@nanmean,curr_data,wf_grp.animal(curr_trials)),1);

        % Average data with movement removede
        % (only plot no-move timepoints with at least N trials and N animals)
        curr_data_no_move = wf_striatum_roi(curr_trials,:,curr_k).* ...
            AP_nanout(wf_t > wf_grp.rxn(curr_trials)-move_leeway);

        min_nomove_trials = 10;
        min_noanimals_trials = 5;

        curr_data_no_move_n = ap.groupfun(@(x) sum(~isnan(x)),curr_data_no_move,wf_grp.animal(curr_trials));
        curr_data_no_move_animal_mean = ap.groupfun(@nanmean,curr_data_no_move,wf_grp.animal(curr_trials));
        curr_data_no_move_animal_mean_mindata = curr_data_no_move_animal_mean.* ...
            AP_nanout(curr_data_no_move_n < min_nomove_trials).* ...
            AP_nanout(sum(curr_data_no_move_n >= min_nomove_trials) < min_noanimals_trials);
        curr_data_no_move_mean = nanmean(curr_data_no_move_animal_mean_mindata,1);
        curr_data_no_move_sem = AP_sem(curr_data_no_move_animal_mean_mindata,1);

        % Plot
        ap.errorfill(wf_t,curr_data_mean,curr_data_sem,[0,0.6,0],0.3,false,2);
        plot(wf_t,curr_data_no_move_mean,'color',[0,0.6,0],'linewidth',2);
        axis off;
        if curr_k == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Fig 2X] Cortex task max stim activity split aross session

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

rxn_cutoff = 0.3; % only plot trials with slow reaction times

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% Split by trial percentile within day
n_split = 3;
split_idx = cell2mat(cellfun(@(n_trials) ...
    ap.quantile_bin(n_trials,n_split), ...
    num2cell(cortex_trials_rec_n),'uni',false));

% (mean across trials, max in time window, mean across animals)
use_trials = wf_grp.rxn > rxn_cutoff;
[activity_mean,activity_mean_grp] = ap.groupfun(@mean,wf_striatum_roi(use_trials,:,:), ...
    [wf_grp.animal(use_trials),plot_day_grp(use_trials),split_idx(use_trials)]);

max_t = wf_t > 0 & wf_t < 0.2;
activity_max = permute(max(activity_mean(:,max_t,:),[],2),[1,3,2]);
[activity_max_avg,activity_max_grp] = ap.nestgroupfun({@mean,@mean},activity_max,activity_mean_grp(:,1),activity_mean_grp(:,2:end));
activity_max_sem = ap.nestgroupfun({@mean,@AP_sem},activity_max,activity_mean_grp(:,1),activity_mean_grp(:,2:end));

figure;
h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day = 1:length(plot_day_bins)-1
        nexttile; hold on;
        curr_data_idx = activity_max_grp(:,1) == curr_day;
        ap.errorfill(1:n_split,activity_max_avg(curr_data_idx,curr_k), ...
            activity_max_sem(curr_data_idx),[0,0.6,0],0.5,false,2);
        axis padded
        if curr_k == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Fig 3X] Striatum PSTH passive

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-1:1,Inf];
plot_days_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

[striatum_mua_avg,striatum_mua_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.kidx]);
striatum_mua_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.kidx]);

unique_stim = unique(striatum_mua_grp.stim);
figure;
h = tiledlayout(n_k,max(plot_days_grp)*length(unique_stim),'TileSpacing','tight');
for curr_k = unique(striatum_mua_avg_grp(:,3))'
    for curr_day = unique(striatum_mua_avg_grp(:,1))'
        for plot_stim = unique_stim'
            nexttile; axis off;
            plot_data = striatum_mua_avg_grp(:,1) == curr_day & ...
                striatum_mua_avg_grp(:,2) == plot_stim & ...
                striatum_mua_avg_grp(:,3) == curr_k;
            if ~any(plot_data)
                continue
            end
            ap.errorfill(psth_t,striatum_mua_avg(plot_data,:),striatum_mua_sem(plot_data,:), ...
                day_colormap(curr_day,:));
        end
    end
end
linkaxes(h.Children,'xy');
xlim([-0.2,0.8]);
title(h,'Striatum');
ap.prettyfig;


%% [Fig 3X] Striatum PSTH max passive

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2:2,Inf];
plot_days_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

[striatum_mua_dayavg,striatum_mua_dayavg_grp] = ...
    ap.groupfun(@mean,striatum_mua, ...
    [striatum_mua_grp.animal,plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.kidx]);

max_t = psth_t > 0 & psth_t < 0.2;
[striatum_mua_max,striatum_mua_max_grp] = ap.nestgroupfun({@mean,@mean}, ...
    max(striatum_mua_dayavg(:,max_t),[],2),striatum_mua_dayavg_grp(:,1), ...
    striatum_mua_dayavg_grp(:,2:end));
striatum_mua_max_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    max(striatum_mua_dayavg(:,max_t),[],2),striatum_mua_dayavg_grp(:,1), ...
    striatum_mua_dayavg_grp(:,2:end));

figure;
h = tiledlayout(n_k,1);
stim_colormap = ap.colormap('BKR',3);
binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');
for curr_k = unique(striatum_mua_grp.kidx)'
    nexttile; hold on; set(gca,'ColorOrder',stim_colormap);
    for curr_stim = unique(striatum_mua_grp.stim)'
        plot_grps = striatum_mua_max_grp(:,2) == curr_stim & ...
            striatum_mua_max_grp(:,3) == curr_k;

        errorbar(binned_days_x,striatum_mua_max(plot_grps), ...
            striatum_mua_max_sem(plot_grps),'linewidth',2);
    end
    axis padded
    xline(0);
end
linkaxes(h.Children,'xy');
title(h,'Striatum');
ap.prettyfig;


%% [Fig 3X] Striatum fraction stim-responsive units

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

% Set days to group
plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

% Grab "responsive" striatal units
unit_responsive_p_thresh = 0.05;
unit_responsive = cell2mat(cellfun(@(x) cell2mat(x)', ...
    horzcat(ephys.unit_resp_p_value{:})','uni',false)) < unit_responsive_p_thresh;
striatum_units_responsive = unit_responsive(cell2mat(striatum_units),:);

% Get fraction of responsive unit (restrict to MSNs)
use_units = striatum_sua_grp.msn;

[unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
    +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
    [plot_day_grp(use_units),striatum_sua_grp.kidx(use_units)]);

unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
    [plot_day_grp(use_units),striatum_sua_grp.kidx(use_units)]);

figure;
h = tiledlayout(n_k,1);
stim_colormap = ap.colormap('BKR',3);
binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');
for curr_k = 1:n_k
    nexttile; hold on;
    set(gca,'ColorOrder',stim_colormap);

    curr_data_idx = unit_responsive_mean_group(:,2) == curr_k;
    errorbar(binned_days_x(unit_responsive_mean_group(curr_data_idx,1)), ...
        unit_responsive_mean(curr_data_idx,:), ...
        unit_responsive_sem(curr_data_idx,:),'linewidth',2);
    ylabel('Frac. responsive units');
    xline(0);
end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Fig 3X] Cortex map passive max

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-1:1,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

stim_t = wf_t > 0 & wf_t < 0.2;
[wf_avg,wf_avg_grp] = ap.groupfun(@mean,cell2mat(wf.V_stim_align), ...
    [wf_grp.animal,plot_day_grp,cell2mat(wf.trial_stim_values)]);

wf_max_px = permute(max(plab.wf.svd2px(U_master,permute(wf_avg(:,stim_t,:),[3,2,1])),[],3),[1,2,4,3]);
[wf_max_px_avg,wf_max_px_avg_grp] = ap.nestgroupfun({@mean,@mean},reshape(wf_max_px,[],size(wf_max_px,3))', ...
    wf_avg_grp(:,1),wf_avg_grp(:,2:3));

figure;
h = tiledlayout(3,max(plot_day_grp),'TileIndexing','ColumnMajor','TileSpacing','none');
for curr_day = unique(wf_max_px_avg_grp(:,1))'
    for curr_stim = unique(wf_max_px_avg_grp(:,2))'
        nexttile;
        imagesc(reshape(wf_max_px_avg( ...
            ismember(wf_max_px_avg_grp,[curr_day,curr_stim],'rows'),:), ...
            size(U_master,[1,2])));
        colormap(ap.colormap('PWG',[],1.5));
        clim([-1,1]*1e-3*2);
        axis image off;
        ap.wf_draw('ccf',[0.5,0.5,0.5]);
    end
end
ap.prettyfig;


%% [Fig 3X] Cortex passive ROI average

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-1:1,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

[wf_striatum_roi_avg,wf_striatum_roi_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi, ...
    wf_grp.animal,[plot_day_grp,cell2mat(wf.trial_stim_values)]);

wf_striatum_roi_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi, ...
    wf_grp.animal,[plot_day_grp,cell2mat(wf.trial_stim_values)]);

figure;
unique_stim = unique(cell2mat(wf.trial_stim_values));
h = tiledlayout(n_k,max(plot_day_grp)*length(unique_stim),'TileSpacing','tight');
for curr_k = 1:size(striatum_wf_roi,3)
    for curr_ld = unique(plot_day_grp)'
        for curr_stim = unique_stim'
            nexttile; axis off;
            curr_data_idx = wf_striatum_roi_avg_grp(:,1) == curr_ld & ...
                wf_striatum_roi_avg_grp(:,2) == curr_stim;

            ap.errorfill(wf_t,wf_striatum_roi_avg(curr_data_idx,:,curr_k), ...
                wf_striatum_roi_sem(curr_data_idx,:,curr_k), ...
                day_colormap(curr_ld,:));
        end
    end
end
linkaxes(h.Children,'xy');
title(h,'Cortex');
ap.prettyfig;


%% [Fig 3X] Cortex passive ROI max

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

[wf_striatum_roi_avg,wf_striatum_roi_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi, ...
    wf_grp.animal,[plot_day_grp,cell2mat(wf.trial_stim_values)]);

wf_striatum_roi_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi, ...
    wf_grp.animal,[plot_day_grp,cell2mat(wf.trial_stim_values)]);

% Plot max in ROI by day
stim_t = wf_t >= 0 & wf_t <= 0.2;

[wf_striatum_roi_dayavg,wf_striatum_roi_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi, ...
    (1:size(wf_striatum_roi,1))', ...
    [wf_grp.animal,plot_day_grp,cell2mat(wf.trial_stim_values)]);

wf_striatum_roi_max = max(wf_striatum_roi_dayavg(:,stim_t,:),[],2);

[wf_striatum_roi_max_avg,wf_striatum_roi_max_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi_max, ...
    wf_striatum_roi_grp(:,1),wf_striatum_roi_grp(:,2:3));

wf_striatum_roi_max_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi_max, ...
    wf_striatum_roi_grp(:,1),wf_striatum_roi_grp(:,2:3));

figure;
h = tiledlayout(n_k,1,'TileSpacing','tight');
stim_colormap = ap.colormap('BKR',3);
binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');
for curr_k = 1:size(striatum_wf_roi,3)
    nexttile; hold on;
    set(gca,'ColorOrder',stim_colormap);
    for curr_stim = unique(unique(wf_striatum_roi_max_avg_grp(:,2)))'
        curr_data_idx = wf_striatum_roi_max_avg_grp(:,2) == curr_stim;

        errorbar(binned_days_x,wf_striatum_roi_max_avg(curr_data_idx,:,curr_k), ...
            wf_striatum_roi_max_sem(curr_data_idx,:,curr_k),'linewidth',2);
        ylabel('\DeltaF/F_0 max');
        axis padded
    end
end
linkaxes(h.Children,'xy');
title(h,'Cortex');
ap.prettyfig;

%% [Fig 4X] Striatal passive by cell type

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-inf,0,inf];
plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

stim_t = psth_t > 0.05 & psth_t < 0.15;

n_stim = size(striatum_sua,3);

% Plot heatmap
for curr_k = 1:n_k
    figure;
    colormap(ap.colormap('WK',[],2));
    h = tiledlayout(n_stim,max(plot_day_grp),'TileSpacing','compact');
    for curr_stim = 1:n_stim
        for curr_day_grp = 1:length(plot_day_bins)-1
            nexttile;

            curr_units = find(striatum_sua_grp.kidx == curr_k & ...
                plot_day_grp == curr_day_grp & ...
                striatum_sua_grp.tan);

            % (sort max across stim)
            [~,sort_idx] = sort(max(mean(striatum_sua(curr_units,stim_t,:),2),[],3),'descend');

            imagesc(psth_t,[],striatum_sua(curr_units(sort_idx),:,curr_stim));
            clim([-1,1])
            xlim([-0.2,0.8])
            title(plot_day_bins(curr_day_grp));

            % (unused: get rec/IDs of cells)
            curr_sorted_unit_coordinate = ...
                [ephys.animal(striatum_sua_grp.rec(curr_units(sort_idx))), ...
                ephys.rec_day(striatum_sua_grp.rec(curr_units(sort_idx))), ...
                num2cell(striatum_sua_grp.unit_id(curr_units(sort_idx)))];

        end
    end
    title(h,sprintf('Striatal cluster %d',curr_k));
    linkaxes(h.Children,'xy');
    ap.prettyfig;
end

% Plot average
figure;
h = tiledlayout(n_k,n_stim*max(plot_day_grp),'TileSpacing','compact');
stim_colormap = ap.colormap('BKR',n_stim);
for curr_k = 1:n_k
    for curr_day_grp = 1:length(plot_day_bins)-1

        curr_units = find(striatum_sua_grp.kidx == curr_k & ...
            plot_day_grp == curr_day_grp & ...
            striatum_sua_grp.tan);

        curr_sua_mean = mean(ap.groupfun(@mean, ...
            striatum_sua(curr_units,:,:),striatum_sua_grp.animal(curr_units)),1);
        curr_sua_sem = AP_sem(ap.groupfun(@mean, ...
            striatum_sua(curr_units,:,:),striatum_sua_grp.animal(curr_units)),1);

        for curr_stim = 1:3
            nexttile;
            ap.errorfill(psth_t,curr_sua_mean(:,:,curr_stim),curr_sua_sem(:,:,curr_stim),stim_colormap(curr_stim,:));
        end

    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;
xlim(h.Children(1),[-0.2,0.8]);

% Stim response by celltype
celltype_order = {'msn','fsi','tan'};
celltype_id = sum([striatum_sua_grp.(celltype_order{1}), ...
    striatum_sua_grp.(celltype_order{2}), ...
    striatum_sua_grp.(celltype_order{3})].*[1,2,3],2);

stim_colormap = ap.colormap('BWR',n_stim);

figure;
h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day_grp = 1:max(plot_day_grp)

        curr_units = find(striatum_sua_grp.kidx == curr_k & ...
            plot_day_grp == curr_day_grp & celltype_id ~= 0);

        curr_unit_mean = squeeze(mean(striatum_sua(curr_units,stim_t,:),2));

        curr_act_mean = ap.nestgroupfun({@nanmean,@nanmean},curr_unit_mean, ...
            striatum_sua_grp.animal(curr_units),celltype_id(curr_units));

        curr_act_sem = ap.nestgroupfun({@nanmean,@AP_sem},curr_unit_mean, ...
            striatum_sua_grp.animal(curr_units),celltype_id(curr_units));

        nexttile; hold on; set(gca,'ColorOrder',stim_colormap);
        celltype_x = reordercats(categorical(celltype_order),celltype_order);
        b = bar(celltype_x,curr_act_mean);
        errorbar(vertcat(b.XEndPoints)',curr_act_mean,curr_act_sem, ...
            'marker','none','linestyle','none','color','k','linewidth',1);

        % ~stats~ %
        rc_diff = ap.nestgroupfun({@nanmean,@nanmean}, ...
            diff(curr_unit_mean(:,[2,3]),[],2), ...
            striatum_sua_grp.animal(curr_units),celltype_id(curr_units));

        n_shuff = 1000;
        rc_diff_shuff = nan(3,n_shuff);
        for curr_shuff = 1:n_shuff
            data_shuff = ap.shake(curr_unit_mean(:,[2,3]),2);
            rc_diff_shuff(:,curr_shuff) = ...
                ap.nestgroupfun({@nanmean,@nanmean}, ...
                diff(data_shuff,[],2), ...
                striatum_sua_grp.animal(curr_units),celltype_id(curr_units));
        end
        stat_rank = tiedrank([rc_diff,rc_diff_shuff]');
        stat_p = 1-stat_rank(1,:)/(n_shuff+1);
        title(sprintf('RvC p = %.2f/%.2f/%.2f',stat_p));

    end
end
linkaxes(h.Children,'xy');
legend({'L','C','R'});
ap.prettyfig;


%% [Fig 4X] Striatal task heatmap by cell type

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

% plot_day_bins = [-Inf,-2:2,Inf];
plot_day_bins = [-Inf,0,Inf];
plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

plot_kidx = [1];

stim_t = psth_t > 0.05 & psth_t < 0.15;

% Plot grouped days
figure;
colormap(ap.colormap('BWR',[],2));
h = tiledlayout(2,max(plot_day_grp),'TileIndexing','column','TileSpacing','compact');
for curr_day_grp = 1:length(plot_day_bins)-1
    nexttile;

    curr_units = find(ismember(striatum_sua_grp.kidx,plot_kidx) & ...
        plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.tan);

    [~,sort_idx] = sort(max(striatum_sua(curr_units,stim_t),[],2),'descend');
    imagesc(psth_t,[],striatum_sua(curr_units(sort_idx),:));
    clim([-1,1])
    title(plot_day_bins(curr_day_grp));

    nexttile;
    curr_sua_mean = mean(ap.groupfun(@mean, ...
        striatum_sua(curr_units,:,:),striatum_sua_grp.animal(curr_units)),1);
    curr_sua_sem = AP_sem(ap.groupfun(@mean, ...
        striatum_sua(curr_units,:,:),striatum_sua_grp.animal(curr_units)),1);
    ap.errorfill(psth_t,curr_sua_mean,curr_sua_sem,'k');

end
linkaxes(h.Children(1:2:end),'xy');
linkaxes(h.Children(2:2:end),'xy'); 
ap.prettyfig;


%% (end timer)

toc







