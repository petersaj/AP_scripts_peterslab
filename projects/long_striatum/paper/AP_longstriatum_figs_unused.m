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

%% [Fig 3X] Cortex passive ROI average

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,2,Inf];
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
h = tiledlayout(n_domains,max(plot_day_grp)*length(unique_stim),'TileSpacing','tight');
for curr_domain = 1:size(striatum_wf_roi,3)
    for curr_ld = unique(plot_day_grp)'
        for curr_stim = unique_stim'
            nexttile; axis off;
            curr_data_idx = wf_striatum_roi_avg_grp(:,1) == curr_ld & ...
                wf_striatum_roi_avg_grp(:,2) == curr_stim;

            ap.errorfill(wf_t,wf_striatum_roi_avg(curr_data_idx,:,curr_domain), ...
                wf_striatum_roi_sem(curr_data_idx,:,curr_domain), ...
                day_colormap(curr_ld,:));
        end
    end
end
linkaxes(h.Children,'xy');
title(h,'Cortex');
ap.prettyfig;

%% [Fig 3X] Striatum PSTH passive

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,2,Inf];
plot_days_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

[striatum_mua_avg,striatum_mua_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);
striatum_mua_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);

unique_stim = unique(striatum_mua_grp.stim);
figure;
h = tiledlayout(n_domains,max(plot_days_grp)*length(unique_stim),'TileSpacing','tight');
for curr_domain = unique(striatum_mua_avg_grp(:,3))'
    for curr_day = unique(striatum_mua_avg_grp(:,1))'
        for plot_stim = unique_stim'
            nexttile; axis off;
            plot_data = striatum_mua_avg_grp(:,1) == curr_day & ...
                striatum_mua_avg_grp(:,2) == plot_stim & ...
                striatum_mua_avg_grp(:,3) == curr_domain;
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


%% [Fig 3X] Striatum fraction stim-responsive units

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

% Set days to group
plot_day_bins = [-Inf,-2,0,2,Inf];
plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

% Grab "responsive" striatal units
unit_responsive_p_thresh = 0.95;
unit_responsive = cell2mat(horzcat(ephys.unit_resp_p_value{:})') > unit_responsive_p_thresh;
striatum_units_responsive = unit_responsive(cell2mat(striatum_units),:);

% Get fraction of responsive unit (restrict to MSNs)
use_units = striatum_sua_grp.msn;

[unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
    +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
    [plot_day_grp(use_units),striatum_sua_grp.domain_idx(use_units)]);

unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
    [plot_day_grp(use_units),striatum_sua_grp.domain_idx(use_units)]);

figure;
h = tiledlayout(n_domains,1);
stim_colormap = ap.colormap('BKR',3);
binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');
for curr_domain = 1:n_domains
    nexttile; hold on;
    set(gca,'ColorOrder',stim_colormap);

    curr_data_idx = unit_responsive_mean_group(:,2) == curr_domain;
    errorbar(binned_days_x(unit_responsive_mean_group(curr_data_idx,1)), ...
        unit_responsive_mean(curr_data_idx,:), ...
        unit_responsive_sem(curr_data_idx,:),'linewidth',2);
    ylabel('Frac. responsive units');
    xline(0);
end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Fig 4X] Passive striatum unit heatmap

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

% Set striatal domains to plot (combines if multiple)
plot_domains = 1:2;

% Set days to group
plot_day_bins = [-inf,-2,0,Inf];
plot_day_grp = discretize(max(-inf,striatum_sua_grp.ld),plot_day_bins);

% Get mean activity in window after stim onset
stim_t = [0,0.2];
stim_use_t = isbetween(psth_t,stim_t(1),stim_t(2));
striatum_sua_tavg = permute(mean(striatum_sua(:,stim_use_t,:),2),[1,3,2]);

% Plot heatmap
celltypes = ["msn","fsi","tan"];
n_stim = size(striatum_sua,3);

figure;
h = tiledlayout(1,length(celltypes),'TileSpacing','tight');
for curr_celltype = celltypes
    h_sub = tiledlayout(h,max(plot_day_grp),n_stim,'TileSpacing','tight');
    [~,curr_celltype_idx] = ismember(curr_celltype,celltypes);
    h_sub.Layout.Tile = curr_celltype_idx;

    colormap(ap.colormap('BWR',[],2));
    title(h_sub,curr_celltype);
    for curr_day_grp = 1:length(plot_day_bins)-1

        % Get units and sorting
        curr_units = find( ...
            ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            plot_day_grp == curr_day_grp & ...
            striatum_sua_grp.(curr_celltype));

        % (sort max stim, then max within stim)
        sort_idx_cell = cell(3,1);
        [~,max_stim] = max(striatum_sua_tavg,[],2);
        for curr_stim_sort = 1:3
            curr_stim_units = find(max_stim(curr_units)==curr_stim_sort);
            [~,curr_sort_idx] = sort(max(striatum_sua_tavg(curr_units(curr_stim_units),curr_stim_sort),[],2),'descend');
            sort_idx_cell{curr_stim_sort} = curr_stim_units(curr_sort_idx);
        end
        sort_idx = cell2mat(sort_idx_cell);
        plot_units = curr_units(sort_idx);

        % (get rec/IDs of cells for single-unit PSTH plotting)
        curr_sorted_unit_coordinate = ...
            [ephys.animal(striatum_sua_grp.rec(plot_units)), ...
            ephys.rec_day(striatum_sua_grp.rec(plot_units)), ...
            num2cell(striatum_sua_grp.unit_id(plot_units))];

        % (y-smooth heatmaps with large number of units)
        max_n_cells = max(ap.groupfun(@sum, ...
            +(ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            striatum_sua_grp.(curr_celltype)),plot_day_grp));
        smooth_n = 5*max_n_cells/500;

        for curr_stim = 1:n_stim
            nexttile(h_sub); 
            imagesc(psth_t,[],movmean(striatum_sua(plot_units,:,curr_stim),[smooth_n,0],1));
            clim([-1,1])
            yline(cumsum(cellfun(@length,sort_idx_cell(1:end-1))),'k');
            axis off;
        end

    end
end
linkaxes(vertcat(h.Children.Children),'x');
xlim(vertcat(h.Children.Children),[0,0.5]);
ap.prettyfig;


% Plot max by stim response
[~,unit_max_stim] = max(striatum_sua_tavg,[],2);

figure;
h = tiledlayout(1,length(celltypes),'TileSpacing','tight');
for curr_celltype = celltypes

    h_sub = tiledlayout(h,n_stim*max(plot_day_grp),1,'TileSpacing','tight');
    [~,curr_celltype_idx] = ismember(curr_celltype,celltypes);
    h_sub.Layout.Tile = curr_celltype_idx;
    title(h_sub,curr_celltype);

    stim_colormap = ap.colormap('BKR',n_stim);
    for curr_day_grp = 1:length(plot_day_bins)-1
        for curr_max_stim = 1:3

            curr_units = ...
                ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
                plot_day_grp == curr_day_grp & striatum_sua_grp.(curr_celltype) & ...
                unit_max_stim == curr_max_stim;

            curr_act_mean = ap.nestgroupfun({@nanmean,@nanmean},striatum_sua_tavg(curr_units,:), ...
                striatum_sua_grp.animal(curr_units),plot_day_grp(curr_units));

            curr_act_sem = ap.nestgroupfun({@nanmean,@AP_sem},striatum_sua_tavg(curr_units,:), ...
                striatum_sua_grp.animal(curr_units),plot_day_grp(curr_units));

            nexttile(h_sub); hold on;
            b = bar(curr_act_mean,'FaceColor','flat','CData',stim_colormap);
            errorbar(curr_act_mean,curr_act_sem, ...
                'marker','none','linestyle','none','color','k','linewidth',1);

        end
    end
    linkaxes(h_sub.Children,'y');

end
linkaxes(vertcat(h.Children.Children),'x');
ap.prettyfig;


% Average stim response grouped by celltype
stim_colormap = ap.colormap('BWR',n_stim);

stim_color = {'KB';'KW';'KR'};
stim_day_colormap = reshape(mat2cell(permute(cell2mat(permute(cellfun(@(stim_color) ...
    ap.colormap(stim_color,length(plot_day_bins)-1),stim_color,'uni',false), ...
    [2,3,1])),[3,2,1]),n_stim,3,ones(length(plot_day_bins)-1,1)),[],1);

figure;
h = tiledlayout(1,length(celltypes),'TileSpacing','compact');
for curr_celltype = celltypes

    curr_units =  ...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        striatum_sua_grp.(curr_celltype);

    curr_act_mean = ap.nestgroupfun({@nanmean,@nanmean},striatum_sua_tavg(curr_units,:), ...
        striatum_sua_grp.animal(curr_units),plot_day_grp(curr_units));

    curr_act_sem = ap.nestgroupfun({@nanmean,@AP_sem},striatum_sua_tavg(curr_units,:), ...
        striatum_sua_grp.animal(curr_units),plot_day_grp(curr_units));

    nexttile; hold on;
    b = bar(curr_act_mean','FaceColor','flat');
    [b.CData] = stim_day_colormap{:};
    errorbar(vertcat(b.XEndPoints)',curr_act_mean',curr_act_sem', ...
        'marker','none','linestyle','none','color','k','linewidth',1);
    title(curr_celltype);

end
linkaxes(h.Children,'x');
ap.prettyfig;

% ~~~ STATS ~~~
fprintf('---- STATS ----\n')

% Shuffle day group for each cell independently (within animal)
for curr_day_grp = 1:length(plot_day_bins)-2

    compare_day_grps = curr_day_grp + [0,1];

    celltype_idx = sum(cell2mat(cellfun(@(x) ...
        striatum_sua_grp.(x),num2cell(celltypes),'uni',false)).* ...
        (1:length(celltypes)),2);

    stat_use_units = celltype_idx ~= 0 & ismember(plot_day_grp,compare_day_grps);

    [stat_data,stat_data_grp] = ap.nestgroupfun({@nanmean,@nanmean},striatum_sua_tavg(stat_use_units,:), ...
        striatum_sua_grp.animal(stat_use_units),[plot_day_grp(stat_use_units),celltype_idx(stat_use_units)]);

    stat_meas = permute(stat_data(stat_data_grp(:,1) == compare_day_grps(2),:) - ...
        stat_data(stat_data_grp(:,1) == compare_day_grps(1),:),[3,2,1]);

    n_shuff = 10000;
    [~,~,shuff_grp] = unique([striatum_sua_grp.animal(stat_use_units),celltype_idx(stat_use_units)],'rows');
    stat_null = nan(n_shuff,n_stim,length(celltypes));
    for curr_shuff = 1:n_shuff

        curr_data_shuff = ap.shake(striatum_sua_tavg(stat_use_units,:),1,shuff_grp);

        [stat_data,stat_data_grp] = ap.nestgroupfun({@nanmean,@nanmean},curr_data_shuff, ...
            striatum_sua_grp.animal(stat_use_units),[plot_day_grp(stat_use_units),celltype_idx(stat_use_units)]);

        stat_null(curr_shuff,:,:) = permute(stat_data(stat_data_grp(:,1) == compare_day_grps(2),:) - ...
            stat_data(stat_data_grp(:,1) == compare_day_grps(1),:),[3,2,1]);

    end

    stat_rank = tiedrank([stat_meas;stat_null]);
    stat_p = 1-stat_rank(1,:,:)/(n_shuff+1);

    unique_stim = unique(striatum_mua_grp.stim)';
    fprintf('Cell types: day grps %d vs %d\n',compare_day_grps);
    stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
    for curr_celltype_idx = 1:length(celltypes)
        for curr_stim_idx = 1:length(unique_stim)
            fprintf('%s, Stim %3.f, p = %.2g%s\n', ...
                celltypes(curr_celltype_idx),unique_stim(curr_stim_idx), ...
                stat_p(1,curr_stim_idx,curr_celltype_idx), ...
                stat_sig(1,curr_stim_idx,curr_celltype_idx));
        end
    end
end

% % ALT STAT unused: shuffle day group for cell type average (within animal)
% compare_day_grps = [2,3];
% 
% celltype_idx = sum(cell2mat(cellfun(@(x) ...
%     striatum_sua_grp.(x),num2cell(celltypes),'uni',false)).* ...
%     (1:length(celltypes)),2);
% 
% stat_use_units = celltype_idx ~= 0 & ismember(plot_day_grp,compare_day_grps);
% 
% [stat_data,stat_data_grp] = ap.groupfun(@nanmean,striatum_sua_tavg(stat_use_units,:), ...
%     [striatum_sua_grp.animal(stat_use_units),plot_day_grp(stat_use_units),celltype_idx(stat_use_units)]);
% 
% stat_meas = permute(diff(reshape(ap.groupfun(@mean,stat_data, ...
%     stat_data_grp(:,2:3))',n_stim,length(celltypes),2),[],3),[3,1,2]);
% 
% n_shuff = 10000;
% [~,~,shuff_grp] = unique(stat_data_grp(:,[1,3]),'rows');
% stat_null = nan(n_shuff,n_stim,length(celltypes));
% for curr_shuff = 1:n_shuff
%     curr_data_shuff = ap.shake(stat_data,1,shuff_grp);
%     stat_null(curr_shuff,:,:) = ...
%         permute(diff(reshape(ap.groupfun(@mean,curr_data_shuff, ...
%         stat_data_grp(:,2:3))',n_stim,length(celltypes),2),[],3),[3,1,2]);
% end
% 
% stat_rank = tiedrank([stat_meas;stat_null]);
% stat_p = 1-stat_rank(1,:,:)/(n_shuff+1);
% 
% unique_stim = unique(striatum_mua_grp.stim)';
% fprintf('Cell types: day grps %d vs %d\n',compare_day_grps);
% stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
% for curr_celltype_idx = 1:length(celltypes)
%     for curr_stim_idx = 1:length(unique_stim)
%         fprintf('%s, Stim %3.f, p = %.2g%s\n', ...
%             celltypes(curr_celltype_idx),unique_stim(curr_stim_idx), ...
%             stat_p(1,curr_stim_idx,curr_celltype_idx), ...
%             stat_sig(1,curr_stim_idx,curr_celltype_idx));              
%     end
% end


%% [Fig 4X] Passive striatum units R vs C response (old version)

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_domains = 1:3;

% Set stim to compare, and colors
compare_stim = [2,3]; % C, R stim
stim_1_color = [0,0,0.8];
stim_2_color = [0.8,0,0];

% Get max response in window
use_t = psth_t > 0.05 & psth_t < 0.15;
unit_psth_max = permute(max(striatum_sua(:,use_t,:),[],2),[1,3,2]);

plot_day_bins = [-Inf,-2,0,Inf];
unit_plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

% Get responsive units
unit_responsive_p_thresh = 0.99;
unit_responsive = cell2mat(horzcat(ephys.unit_resp_p_value{:})') > unit_responsive_p_thresh;
striatum_units_responsive = unit_responsive(cell2mat(striatum_units),:);

% Loop through celltypes, plot scatter and stim responsive overlap
striatum_celltypes = ["msn","fsi","tan"];

figure;
h = tiledlayout(length(striatum_celltypes),1,'TileSpacing','tight');
example_units = nan(length(striatum_celltypes),2);
for curr_celltype = striatum_celltypes

    h_sub = tiledlayout(h,1,length(plot_day_bins),'TileSpacing','tight');
    [~,curr_celltype_idx] = ismember(curr_celltype,striatum_celltypes);
    h_sub.Layout.Tile = curr_celltype_idx;
    title(h_sub,curr_celltype);

    % Scatter
    h_scatter = gobjects(length(plot_day_bins)-1,1);
    for curr_day_grp = 1:length(plot_day_bins)-1
        curr_units = ...
            ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            unit_plot_day_grp == curr_day_grp & ...
            any(striatum_units_responsive(:,compare_stim),2) & ...
            striatum_sua_grp.(curr_celltype);

        curr_dot_size = 10;
        curr_dot_color = ...
            stim_1_color.*striatum_units_responsive(curr_units,compare_stim(1)) + ...
            stim_2_color.*striatum_units_responsive(curr_units,compare_stim(2));

        h_scatter(curr_day_grp) = nexttile(h_sub); hold on;
        scatter(unit_psth_max(curr_units,2),unit_psth_max(curr_units,3), ...
            curr_dot_size,curr_dot_color,'filled');
    end

    linkaxes(h_scatter,'xy');
    ax_range = [min([xlim,ylim]),max([xlim,ylim])];
    xlim(ax_range);ylim(ax_range);
    axis(h_scatter,'square')
    arrayfun(@(x) line(h_scatter(x),ax_range,ax_range,'color','k'),1:length(h_scatter));

    % Grab example units to plot
    plot_unit_prctile = 90;
    plot_day_grp = 2;
    for curr_stim = compare_stim

        curr_units = ...
            ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            unit_plot_day_grp == plot_day_grp & ...
            striatum_units_responsive(:,curr_stim) & ...
            ~striatum_units_responsive(:,setxor(curr_stim,compare_stim)) & ...
            striatum_sua_grp.(curr_celltype);

        curr_unit_idx = find(curr_units);
        [~,sort_idx] = sort(unit_psth_max(curr_units,curr_stim));
        
        curr_example_unit = curr_unit_idx(sort_idx(ceil(sum(curr_units).*plot_unit_prctile/100)));

        % (circle example units on scatter)
        scatter(h_scatter(plot_day_grp), ...
            unit_psth_max(curr_example_unit,2),unit_psth_max(curr_example_unit,3), ...
            45,[0.5,0.5,0.5],'linewidth',2);

        [~,curr_stim_idx] = ismember(curr_stim,compare_stim);
        example_units(curr_celltype_idx,curr_stim_idx) = curr_example_unit;
        
    end

    % Frac R/C/R+C as stacked barplot
    use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        striatum_sua_grp.(curr_celltype);

    striatum_units_responsive_1_12_2 = permute(all(striatum_units_responsive(:,compare_stim) == ...
        cat(3,[1,0],[1,1],[0,1]),2),[1,3,2]);

    [unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
        +striatum_units_responsive_1_12_2(use_units,:), ...
        striatum_sua_grp.animal(use_units), ...
        unit_plot_day_grp(use_units));

    unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
        +striatum_units_responsive_1_12_2(use_units,:), ...
        striatum_sua_grp.animal(use_units), ...
        unit_plot_day_grp(use_units));

    nexttile(h_sub); hold on;
    set(gca,'ColorOrder',[stim_1_color;stim_1_color+stim_2_color;stim_2_color]);
    bar(unit_responsive_mean,'stacked');
    ylabel('Frac. responsive units');
    legend(["C","C+R","R"]);

    % % ~~~ STATS ~~~
    % fprintf('---- STATS ----\n')
    % n_shuff = 10000;
    % 
    % % Compare R+C-responsive overlap to shuffling R/C responsiveness
    % use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    %     striatum_sua_grp.(curr_celltype);
    % 
    % [~,~,overlap_shuff_grp] = unique([striatum_sua_grp.animal, ...
    %     unit_plot_day_grp, ...
    %     ismember(striatum_sua_grp.domain_idx,plot_domains), ...
    %     any(striatum_units_responsive(:,compare_stim),2) & ...
    %     striatum_sua_grp.(curr_celltype)],'rows');
    % 
    % responsive_overlap_meas = ap.nestgroupfun({@mean,@mean}, ...
    %     +(all(striatum_units_responsive(use_units,compare_stim),2)), ...
    %     striatum_sua_grp.animal(use_units),unit_plot_day_grp(use_units));
    % 
    % responsive_overlap_shuff = nan(length(plot_day_bins)-1,n_shuff);
    % for curr_shuff = 1:n_shuff
    %     striatum_units_responsive_shuff = ap.shake(striatum_units_responsive(use_units,:),1,overlap_shuff_grp(use_units));
    %     responsive_overlap_shuff(:,curr_shuff) = ap.nestgroupfun({@mean,@mean}, ...
    %         +(all(striatum_units_responsive_shuff(:,compare_stim),2)),striatum_sua_grp.animal(use_units), ...
    %         unit_plot_day_grp(use_units));
    % end
    % 
    % stat_rank = tiedrank([responsive_overlap_meas,responsive_overlap_shuff]')';
    % stat_p = stat_rank(:,1)/(n_shuff+1);
    % 
    % fprintf('Overlap stim %d + %d < shuffle:\n',compare_stim);
    % stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
    % for curr_day = 1:length(plot_day_bins)-1
    %     fprintf('%s day grp %d, p = %.2g%s\n',curr_celltype,curr_day, ...
    %         stat_p(curr_day),stat_sig(curr_day));
    % end
    % 
    % 
    % % Compare R-responsive fraction across days
    % for curr_compare_day = 1:length(plot_day_bins)-2
    % 
    %     compare_day_grps = curr_compare_day+[0,1];
    % 
    %     use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    %         ismember(unit_plot_day_grp,compare_day_grps) & ...
    %         striatum_sua_grp.(curr_celltype);
    % 
    %     [~,~,r_shuff_grp] = unique([striatum_sua_grp.animal, ...
    %         ismember(striatum_sua_grp.domain_idx,plot_domains), ...
    %         striatum_sua_grp.(curr_celltype)],'rows');
    % 
    %     r_frac_meas = diff(ap.nestgroupfun({@mean,@mean}, ...
    %         +(all(striatum_units_responsive(use_units,3),2)), ...
    %         striatum_sua_grp.animal(use_units),unit_plot_day_grp(use_units)));
    % 
    %     r_frac_shuff = nan(n_shuff,1);
    %     for curr_shuff = 1:n_shuff
    %         unit_plot_day_grp_shuff = ap.shake(unit_plot_day_grp(use_units,:),1,r_shuff_grp(use_units));
    %         r_frac_shuff(curr_shuff) = diff(ap.nestgroupfun({@mean,@mean}, ...
    %             +(all(striatum_units_responsive(use_units,3),2)), ...
    %             striatum_sua_grp.animal(use_units),unit_plot_day_grp_shuff));
    %     end
    % 
    %     stat_rank = tiedrank([r_frac_meas;r_frac_shuff]);
    %     stat_p = 1-stat_rank(1)/(n_shuff+1);
    % 
    %     stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
    %     fprintf('%s R-frac day %d vs %d: p = %.2g%s\n',curr_celltype,compare_day_grps,stat_p,stat_sig);
    % 
    % end
end
ap.prettyfig;


% Plot example PSTHs
figure;
h_units = tiledlayout(1,numel(example_units));
for curr_unit = reshape(example_units',1,[])
    animal = ephys.animal{striatum_sua_grp.rec(curr_unit)};
    rec_day = ephys.rec_day{striatum_sua_grp.rec(curr_unit)};
    unit_id = striatum_sua_grp.unit_id(curr_unit);
    AP_longstriatum_psth_fig(animal,rec_day,unit_id,true,h_units);    
end
linkaxes(vertcat(h_units.Children.Children),'x');
ap.prettyfig;


%% [Supp. Fig 1x] Cortex image passive max

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,2,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

stim_t = wf_t > 0 & wf_t < 0.2;
[wf_avg,wf_avg_grp] = ap.groupfun(@mean,cell2mat(wf.V_event_align), ...
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
        clim([-1,1]*3e-3);
        axis image off;
        ap.wf_draw('ccf',[0.5,0.5,0.5]);
    end
end
ap.prettyfig;