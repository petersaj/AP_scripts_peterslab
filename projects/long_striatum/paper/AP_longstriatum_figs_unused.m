%% [Fig 1X] Cortex map depths pre/post learning

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

n_depths = 6;

ctx_map_depthmean = cellfun(@(x) ap.groupfun(@mean,x,[],[],discretize(1:size(x,3),n_depths)), ...
    all_ctx_maps_to_str.cortex_kernel_px,'uni',false);

ctx_map_depthmean_prelearn = nanmean(cat(4,ctx_map_depthmean{bhv.days_from_learning<0 & cellfun(@(x) size(x,3)==n_depths,ctx_map_depthmean)}),4);
ctx_map_depthmean_postlearn = nanmean(cat(4,ctx_map_depthmean{bhv.days_from_learning>=0 & cellfun(@(x) size(x,3)==n_depths,ctx_map_depthmean)}),4);
ap.imscroll(cat(4,ctx_map_depthmean_prelearn,ctx_map_depthmean_postlearn));
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG',[],1.5));


figure;tiledlayout(n_depths,2,'TileSpacing','none');
colormap(ap.colormap('PWG',[],1.5));
for curr_depth = 1:n_depths
    nexttile;
    imagesc(ctx_map_depthmean_prelearn(:,:,curr_depth));
    axis image off;
    ap.wf_draw('ccf',[0.5,0.5,0.5]);
    clim([-1,1]*0.05);

    nexttile;
    imagesc(ctx_map_depthmean_postlearn(:,:,curr_depth));
    axis image off;
    ap.wf_draw('ccf',[0.5,0.5,0.5]);
    clim([-1,1]*0.05);
end
ap.prettyfig;


%% [Fig 1X] Striatal MUA pre/post learning (using final domain_idx)

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

figure; h = tiledlayout(n_domains,3,'TileSpacing','tight');
for curr_depth = 1:n_domains

    curr_trials = striatum_mua_grp.domain_idx == curr_depth;

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


%% [Fig 1X] Striatal MUA pre/post learning (using starting domain_idx)

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

figure; h = tiledlayout(n_domains,3,'TileSpacing','tight');
for curr_depth = 1:n_domains

    curr_trials = striatum_mua_grp.domain_idx == curr_depth;

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


%% [Fig 2X] Striatum task average PSTH w/ and w/o movement

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,Inf];
plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

% Plot average activity in trial (NaN-out movement times)
figure; h = tiledlayout(n_domains,max(plot_day_grp),'TileSpacing','none');
move_leeway = 0.1; % time pre-movement to exclude
act_stim_max = cell(n_domains,length(unique(plot_day_grp)));
for curr_domain = 1:n_domains
    for curr_day = unique(plot_day_grp)'
        
        curr_trials = find(plot_day_grp == curr_day & striatum_mua_grp.domain_idx == curr_domain);      
        
        % Average all data
        curr_data = striatum_mua(curr_trials,:,1);
        curr_data_no_move = curr_data.* ...
            ap.nanout(psth_t > striatum_mua_grp.rxn(curr_trials)-move_leeway);

        curr_data_mean = mean(ap.groupfun(@nanmean,curr_data,striatum_mua_grp.animal(curr_trials)),1);
        curr_data_sem = AP_sem(ap.groupfun(@nanmean,curr_data,striatum_mua_grp.animal(curr_trials)),1);

        % Average data with movement removed 
        curr_data_no_move_animal_mean = ap.groupfun(@nanmean,curr_data_no_move,striatum_mua_grp.animal(curr_trials));
        curr_data_no_move_mean = nanmean(curr_data_no_move_animal_mean,1);
        curr_data_no_move_sem = AP_sem(curr_data_no_move_animal_mean,1);

        % Plot w/ movement (shaded) and w/o movement (line)
        nexttile; hold on;
        ap.errorfill(psth_t,curr_data_mean,curr_data_sem,'k',0.3,false,2);
        plot(psth_t,curr_data_no_move_mean,'k','linewidth',2);
        axis off;
        if curr_domain == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end

        % Store max response to stim for stats
        stim_t = [0,0.2];
        act_stim_max{curr_domain,curr_day} = ...
            max(curr_data_no_move_animal_mean(:, ...
            isbetween(psth_t,stim_t(1),stim_t(2))),[],2);

    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;

% ~stats~
% compare_days = [-inf,-2];
% [~,compare_day_grps] = ismember(compare_days,plot_day_bins);
compare_day_grps = [1,2];

stat_data = act_stim_max(:,compare_day_grps);
stat_meas = diff(cellfun(@mean,stat_data),[],2);

stat_domain_idx = cellfun(@(data,domain) repmat(domain,length(data),1), ...
    stat_data,repmat(num2cell(1:n_domains)',1,2),'uni',false);

n_shuff = 10000;
stat_null = nan(n_domains,n_shuff);
for curr_shuff = 1:n_shuff
    stat_data_shuff = reshape(mat2cell(ap.shake(vertcat(stat_data{:}),1, ...
        vertcat(stat_domain_idx{:})), ...
        cellfun(@length,stat_data(:))),size(stat_data));

    stat_null(:,curr_shuff) = diff(cellfun(@mean,stat_data_shuff),[],2);
end

stat_rank = tiedrank([stat_meas,stat_null]')';
stat_p = 1-stat_rank(:,1)/(n_shuff+1);
stat_sig = discretize(stat_p < 0.05,[0,1,Inf],{'','*'});
fprintf('Day grps %d,%d:\n',compare_day_grps);
for curr_domain = 1:n_domains
    fprintf('D%d p = %.2g%s\n',curr_domain,stat_p(curr_domain),stat_sig{curr_domain});
end


%% [Fig 2X] Cortex task average PSTH w/ and w/o movement

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% Plot average activity in trial (NaN-out movement times)
figure; h = tiledlayout(n_domains,max(plot_day_grp),'TileSpacing','none');
move_leeway = 0.1;
for curr_domain = 1:n_domains
    for curr_day = unique(plot_day_grp)'
        
        curr_trials = find(plot_day_grp == curr_day);      
        
        % Average all data
        curr_data = wf_striatum_roi(curr_trials,:,curr_domain,1);

        curr_data_mean = mean(ap.groupfun(@nanmean,curr_data,wf_grp.animal(curr_trials)),1);
        curr_data_sem = AP_sem(ap.groupfun(@nanmean,curr_data,wf_grp.animal(curr_trials)),1);

        % Average data with movement removede
        % (only plot no-move timepoints with at least N trials and N animals)
        curr_data_no_move = curr_data.* ...
            ap.nanout(wf_t > wf_grp.rxn(curr_trials)-move_leeway);
      
        curr_data_no_move_animal_mean = ap.groupfun(@nanmean,curr_data_no_move,wf_grp.animal(curr_trials));
        curr_data_no_move_mean = nanmean(curr_data_no_move_animal_mean,1);
        curr_data_no_move_sem = AP_sem(curr_data_no_move_animal_mean,1);

        % Plot
        nexttile; hold on;
        ap.errorfill(wf_t,curr_data_mean,curr_data_sem,[0,0.6,0],0.3,false,2);
        plot(wf_t,curr_data_no_move_mean,'color',[0,0.6,0],'linewidth',2);
        axis off;
        if curr_domain == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end

        % Store max response to stim for stats
        stim_t = [0,0.2];
        act_stim_max{curr_domain,curr_day} = ...
            max(curr_data_no_move_animal_mean(:, ...
            isbetween(wf_t,stim_t(1),stim_t(2))),[],2);

    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;


% ~stats~
% compare_days = [-inf,-2];
% [~,compare_day_grps] = ismember(compare_days,plot_day_bins);
compare_day_grps = [1,2];

stat_data = act_stim_max(:,compare_day_grps);
stat_meas = diff(cellfun(@mean,stat_data),[],2);

stat_domain_idx = cellfun(@(data,domain) repmat(domain,length(data),1), ...
    stat_data,repmat(num2cell(1:n_domains)',1,2),'uni',false);

n_shuff = 10000;
stat_null = nan(n_domains,n_shuff);
for curr_shuff = 1:n_shuff
    stat_data_shuff = reshape(mat2cell(ap.shake(vertcat(stat_data{:}),1, ...
        vertcat(stat_domain_idx{:})), ...
        cellfun(@length,stat_data(:))),size(stat_data));

    stat_null(:,curr_shuff) = diff(cellfun(@mean,stat_data_shuff),[],2);
end

stat_rank = tiedrank([stat_meas,stat_null]')';
stat_p = 1-stat_rank(:,1)/(n_shuff+1);
stat_sig = discretize(stat_p < 0.05,[0,1,Inf],{'','*'});
fprintf('Day grps %d,%d:\n',compare_day_grps);
for curr_domain = 1:n_domains
    fprintf('ROI%d p = %.2g%s\n',curr_domain,stat_p(curr_domain),stat_sig{curr_domain});
end



%% [Fig 2X] Cortex + striatal task max, split within session (OLD)
% (old version: splits all days and averages n = sessions if multiple days
% in a group)

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

rxn_cutoff = 0.3; % only plot trials with slow reaction times

plot_day_bins = [-Inf,-1,0,Inf];
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);


% Split by trial percentile within day
n_split = 1;

striatum_split_idx = cell2mat(cellfun(@(n_trials,n_mua) ...
    reshape(repmat(ap.quantile_bin(n_trials,n_split),1,n_mua),[],1), ...
    num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

cortex_split_idx = cell2mat(cellfun(@(n_trials) ...
    ap.quantile_bin(n_trials,n_split), ...
    num2cell(cortex_trials_rec_n),'uni',false));

% Get activity: trial mean > time max > animal mean 
stim_t = [0,0.2];

% (striatum)
striatum_use_trials = striatum_mua_grp.rxn > rxn_cutoff;
[striatum_activity_mean,striatum_activity_mean_grp] = ap.groupfun(@mean,striatum_mua(striatum_use_trials,:,1), ...
    [striatum_mua_grp.animal(striatum_use_trials),striatum_plot_day_grp(striatum_use_trials), ...
    striatum_split_idx(striatum_use_trials),striatum_mua_grp.domain_idx(striatum_use_trials)]);

striatum_activity_max = max(striatum_activity_mean(:,isbetween(psth_t,stim_t(1),stim_t(2))),[],2);
[striatum_activity_max_avg,striatum_activity_max_grp] = ap.nestgroupfun({@mean,@mean},striatum_activity_max,striatum_activity_mean_grp(:,1),striatum_activity_mean_grp(:,2:end));
striatum_activity_max_sem = ap.nestgroupfun({@mean,@AP_sem},striatum_activity_max,striatum_activity_mean_grp(:,1),striatum_activity_mean_grp(:,2:end));

% (cortex)
cortex_use_trials = wf_grp.rxn > rxn_cutoff;
[cortex_activity_mean,cortex_activity_mean_grp] = ap.groupfun(@mean,wf_striatum_roi(cortex_use_trials,:,:,1), ...
    [wf_grp.animal(cortex_use_trials),cortex_plot_day_grp(cortex_use_trials),cortex_split_idx(cortex_use_trials)]);

cortex_activity_max = permute(max(cortex_activity_mean(:,isbetween(wf_t,stim_t(1),stim_t(2)),:),[],2),[1,3,2]);
[cortex_activity_max_avg,cortex_activity_max_grp] = ap.nestgroupfun({@mean,@mean},cortex_activity_max,cortex_activity_mean_grp(:,1),cortex_activity_mean_grp(:,2:end));
cortex_activity_max_sem = ap.nestgroupfun({@mean,@AP_sem},cortex_activity_max,cortex_activity_mean_grp(:,1),cortex_activity_mean_grp(:,2:end));

% Plot cortex and striatum overlaid
figure;
h = tiledlayout(n_domains,length(plot_day_bins)-1,'TileSpacing','compact');
for curr_domain = 1:n_domains
    for curr_day = 1:length(plot_day_bins)-1
        nexttile; hold on;

        yyaxis left;
        curr_cortex_data_idx = cortex_activity_max_grp(:,1) == curr_day;
        ap.errorfill(1:n_split,cortex_activity_max_avg(curr_cortex_data_idx,curr_domain), ...
            cortex_activity_max_sem(curr_cortex_data_idx),[0,0.6,0],0.5,false,2);

        yyaxis right;
        curr_striatum_data_idx = striatum_activity_max_grp(:,1) == curr_day & ...
            striatum_activity_max_grp(:,3) == curr_domain;
        errorbar(striatum_activity_max_avg(curr_striatum_data_idx), ...
            striatum_activity_max_sem(curr_striatum_data_idx),'k','linewidth',2);

        if curr_domain == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end

% Set axis limits
for curr_ax = 1:length(h.Children)
    yyaxis(h.Children(curr_ax),'left');
end
linkaxes(h.Children);
ylim(h(1).Children,[0,7e-3]);
[h.Children.YColor] = deal([0,0.5,0]);

for curr_ax = 1:length(h.Children)
    yyaxis(h.Children(curr_ax),'right');
end
linkaxes(h.Children);
ylim(h(1).Children,[-0.2,6]);
[h.Children.YColor] = deal([0,0,0]);

xlim(h(1).Children,[0.8,n_split+0.2])

ap.prettyfig;

%% [Fig 2X] Cortex task trial heatmap (reaction-sorted)

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% Plot heatmaps sorted by reaction times
heatmap_smooth = [20,1]; % ([trials,time] to smooth for graphics)
figure; tiledlayout(n_domains,max(plot_day_grp),'TileSpacing','none');
for curr_domain = 1:n_domains
    for curr_day = unique(plot_day_grp)'
        nexttile;
        curr_trials = find(plot_day_grp == curr_day);
       
        [sorted_rxn,sort_idx] = sort(wf_grp.rxn(curr_trials));
        imagesc(wf_t,[],movmean(wf_striatum_roi(curr_trials(sort_idx),:,curr_domain,1),heatmap_smooth));
        colormap(ap.colormap('PWG'));
        clim(max(abs(clim)).*[-1,1]);
        axis off;

        hold on
        xline(0,'color','r');
        plot(sorted_rxn,1:length(curr_trials),'b');
        xlim(prctile(psth_t,[0,100]));

        if curr_domain == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end
ap.prettyfig;

%% [Fig 2X] Striatum task trial heatmap (reaction-sorted)

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,Inf];
plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

heatmap_smooth = [20,1]; % ([trials,time] to smooth for graphics)
figure; tiledlayout(n_domains,max(plot_day_grp),'TileSpacing','none');
for curr_domain = 1:n_domains
    for curr_day = unique(plot_day_grp)'
        nexttile;
        curr_trials = find(plot_day_grp == curr_day & striatum_mua_grp.domain_idx == curr_domain);
       
        [sorted_rxn,sort_idx] = sort(striatum_mua_grp.rxn(curr_trials));      
        imagesc(psth_t,[],movmean(striatum_mua(curr_trials(sort_idx),:,1),heatmap_smooth));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;

        hold on
        xline(0,'color','r');
        plot(sorted_rxn,1:length(curr_trials),'b');
        xlim(prctile(psth_t,[0,100]));

        if curr_domain == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end
ap.prettyfig;




%% [Fig 2X] Cortex ROIs

%%% Load data for figure
load_dataset = 'noact';
AP_longstriatum_load_data;
%%%

figure;tiledlayout(n_domains,1,'tilespacing','none')
for curr_domain = 1:n_domains
    nexttile;
    imagesc(striatum_wf_roi(:,:,curr_domain));
    axis image off; ap.wf_draw('ccf',[0.5,0.5,0.5]);
    colormap(ap.colormap('WG'));
end
ap.prettyfig;

%% [Fig 3X] Cortex passive ROI max

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,2,Inf];
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
h = tiledlayout(n_domains,1,'TileSpacing','tight');
stim_colormap = ap.colormap('BKR',3);
binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');
for curr_domain = 1:size(striatum_wf_roi,3)
    nexttile; hold on;
    set(gca,'ColorOrder',stim_colormap);
    for curr_stim = unique(unique(wf_striatum_roi_max_avg_grp(:,2)))'
        curr_data_idx = wf_striatum_roi_max_avg_grp(:,2) == curr_stim;

        errorbar(binned_days_x,wf_striatum_roi_max_avg(curr_data_idx,:,curr_domain), ...
            wf_striatum_roi_max_sem(curr_data_idx,:,curr_domain),'linewidth',2);
        ylabel('\DeltaF/F_0 max');
        axis padded
    end
end
linkaxes(h.Children,'xy');
title(h,'Cortex');
ap.prettyfig;

% ~stats~
% compare_days = [-inf,-2];
% [~,compare_day_grps] = ismember(compare_days,plot_day_bins);
compare_day_grps = [1,2];

curr_data_idx = ismember(wf_striatum_roi_grp(:,2),compare_day_grps);

[stat_meas,stat_grp] = ap.nestgroupfun({@mean,@diff},wf_striatum_roi_max(curr_data_idx,:,:), ...
    wf_striatum_roi_grp(curr_data_idx,2),wf_striatum_roi_grp(curr_data_idx,3));

[~,~,shuff_grp] = unique(wf_striatum_roi_grp(curr_data_idx,[1,3]),'rows');
n_shuff = 1000;
stat_null = nan(size(stat_meas,1),n_shuff,n_domains);
for curr_shuff = 1:n_shuff
    curr_data_shuff = ap.shake(wf_striatum_roi_max(curr_data_idx,:,:),1,shuff_grp);
    stat_null(:,curr_shuff,:) = ap.nestgroupfun({@mean,@diff},curr_data_shuff, ...
        wf_striatum_roi_grp(curr_data_idx,2),wf_striatum_roi_grp(curr_data_idx,3));
end

stat_rank = permute(tiedrank(permute([stat_meas,stat_null],[2,1,3])),[2,1,3]);
stat_p = 1-stat_rank(:,1,:)/(n_shuff+1);

fprintf('Stats: day grps %d vs %d\n',compare_days);
stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
for curr_domain = 1:n_domains
    for curr_stim = unique(stat_grp(:,1))'
        curr_stat_idx = ismember(stat_grp,curr_stim,'rows');
        fprintf('ROI%d, Stim %3.f, p = %.2g%s\n', ...
            curr_domain,stat_grp(curr_stat_idx),stat_p(curr_stat_idx,1,curr_domain),stat_sig(curr_stat_idx,1,curr_domain));
    end
end

%% [Fig 3X] Striatum PSTH max passive

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,2,Inf];
plot_days_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

[striatum_mua_dayavg,striatum_mua_dayavg_grp] = ...
    ap.groupfun(@mean,striatum_mua, ...
    [striatum_mua_grp.animal,plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);

max_t = psth_t >= 0 & psth_t <= 0.2;
striatum_mua_dayavg_tmax = max(striatum_mua_dayavg(:,max_t),[],2);

[striatum_mua_max,striatum_mua_max_grp] = ap.nestgroupfun({@mean,@mean}, ...
    striatum_mua_dayavg_tmax,striatum_mua_dayavg_grp(:,1), ...
    striatum_mua_dayavg_grp(:,2:end));
striatum_mua_max_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    striatum_mua_dayavg_tmax,striatum_mua_dayavg_grp(:,1), ...
    striatum_mua_dayavg_grp(:,2:end));

figure;
h = tiledlayout(n_domains,1);
stim_colormap = ap.colormap('BKR',3);
binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');
for curr_domain = unique(striatum_mua_grp.domain_idx)'
    nexttile; hold on; set(gca,'ColorOrder',stim_colormap);
    for curr_stim = unique(striatum_mua_grp.stim)'
        plot_grps = striatum_mua_max_grp(:,2) == curr_stim & ...
            striatum_mua_max_grp(:,3) == curr_domain;

        errorbar(binned_days_x,striatum_mua_max(plot_grps), ...
            striatum_mua_max_sem(plot_grps),'linewidth',2);
    end
    axis padded
    xline(0);
end
linkaxes(h.Children,'xy');
title(h,'Striatum');
ap.prettyfig;


% ~stats~
% compare_days = [-inf,-2];
% [~,compare_day_grps] = ismember(compare_days,plot_day_bins);
compare_day_grps = [1,2];

curr_data_idx = ismember(striatum_mua_dayavg_grp(:,2),compare_day_grps);

[stat_meas,stat_grp] = ap.nestgroupfun({@mean,@diff},striatum_mua_dayavg_tmax(curr_data_idx), ...
    striatum_mua_dayavg_grp(curr_data_idx,2),striatum_mua_dayavg_grp(curr_data_idx,3:end));

[~,~,shuff_grp] = unique(striatum_mua_dayavg_grp(curr_data_idx,[1,3,4]),'rows');
n_shuff = 1000;
stat_null = nan(length(stat_meas),n_shuff);
for curr_shuff = 1:n_shuff
    curr_data_shuff = ap.shake(striatum_mua_dayavg_tmax(curr_data_idx),1,shuff_grp);
    stat_null(:,curr_shuff) = ap.nestgroupfun({@mean,@diff},curr_data_shuff, ...
        striatum_mua_dayavg_grp(curr_data_idx,2),striatum_mua_dayavg_grp(curr_data_idx,3:end));
end

stat_rank = tiedrank([stat_meas,stat_null]')';
stat_p = 1-stat_rank(:,1)/(n_shuff+1);

fprintf('Stats: day grps %d vs %d\n',compare_days);
stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
for curr_domain = 1:n_domains
    for curr_stim = unique(stat_grp(:,1))'
        curr_stat_idx = ismember(stat_grp,[curr_stim,curr_domain],'rows');
        fprintf('D%d, Stim %3.f, p = %.2g%s\n', ...
            curr_domain,stat_grp(curr_stat_idx,1),stat_p(curr_stat_idx),stat_sig(curr_stat_idx));
    end
end