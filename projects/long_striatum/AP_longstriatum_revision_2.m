%% Analyses for reviewers (cleaned up from _1 scratch)

%% R1 p1: Correct performance rate by day

%%% Load non-activity data
load_dataset = 'noact';
Marica_2025.figures.load_data;
%%%

% Plot reaction time and association index, split within day
n_daysplit = 3;

outcome_mean_daysplit = cell2mat(cellfun(@(x,idx) ap.groupfun(@mean,+x, ...
    ap.quantile_bin(length(x),n_daysplit)),bhv.trial_outcome,'uni',false)')';

[outcome_daysplit_mean,outcome_group_x] = ...
    ap.groupfun(@mean,outcome_mean_daysplit,bhv.days_from_learning);
outcome_daysplit_sem = ap.groupfun(@AP_sem,outcome_mean_daysplit, ...
    bhv.days_from_learning);

plot_days = -3:2;
plot_day_idx = ismember(outcome_group_x,plot_days);

figure;
outcome_group_x_daysplit = outcome_group_x+(0:n_daysplit)./n_daysplit;
errorbar(reshape(outcome_group_x_daysplit(plot_day_idx,:)',[],1), ...
    reshape(padarray(outcome_daysplit_mean(plot_day_idx,:),[0,1],nan,'post')',[],1), ...
    reshape(padarray(outcome_daysplit_sem(plot_day_idx,:),[0,1],nan,'post')',[],1),'k','linewidth',2);
xline(0,'r');
ylabel('Fraction correct');
xlabel('Day from learning');
ap.prettyfig;

%% R1 p3 / R3 M1: mPFC-striatum timing
% CLEAN UP

load_dataset_retain = true;

% ~~ Set up data structure params
load_dataset = 'passive';
Marica_2025.figures.load_data;

% Set up parameters for activity grids [animal x day x domain x stim]
data_grid_params = struct;

data_grid_params.stim_t = [0,0.2];
data_grid_params.cortex_stim_t = isbetween(wf_t,data_grid_params.stim_t(1),data_grid_params.stim_t(2));
data_grid_params.striatum_stim_t = isbetween(psth_t,data_grid_params.stim_t(1),data_grid_params.stim_t(2));
data_grid_params.rxn_cutoff = 0.3;

data_grid_params.ld_bins = [-Inf,-2:0,1,Inf];
data_grid_params.ld_unique = 1:(length(data_grid_params.ld_bins)-1);
data_grid_params.grid_size = [length(unique(bhv.animal)),length(data_grid_params.ld_unique),n_domains];

% Set up data structure
data_grids = struct;

clearvars -except load_dataset_retain data_grid_params data_grids

% ~~ Get task activity
load_dataset = 'task';
Marica_2025.figures.load_data;

% (striatum)
striatum_use_trials = striatum_mua_grp.rxn > data_grid_params.rxn_cutoff;
striatum_ld_idx = discretize(striatum_mua_grp.ld,data_grid_params.ld_bins);

[striatum_rec,striatum_rec_grp] = ap.groupfun(@nanmean,striatum_mua(:,:,1), ...
    [striatum_mua_grp.animal,striatum_ld_idx,striatum_mua_grp.domain_idx].* ...
    ap.nanout(~(striatum_use_trials & ~isnan(striatum_ld_idx))));
striatum_rec_tmax = max(striatum_rec(:,data_grid_params.striatum_stim_t),[],2);

data_grids.striatum_task = accumarray(striatum_rec_grp,striatum_rec_tmax,data_grid_params.grid_size,[],NaN);

% (widefield)
wf_use_trials = wf_grp.rxn > data_grid_params.rxn_cutoff;
wf_ld_idx = discretize(wf_grp.ld,data_grid_params.ld_bins);

[wf_roi_rec,wf_roi_rec_grp] = ap.groupfun(@nanmean,wf_striatum_roi(:,:,:,1), ...
    [wf_grp.animal,wf_ld_idx].* ...
    ap.nanout(~(wf_use_trials & ~isnan(wf_ld_idx))));
wf_roi_rec_tmax = permute(max(wf_roi_rec(:,data_grid_params.cortex_stim_t,:),[],2),[1,3,2]);

data_grids.wf_roi_task = cell2mat(permute(arrayfun(@(domain) accumarray(wf_roi_rec_grp, ...
    wf_roi_rec_tmax(:,domain),data_grid_params.grid_size(1:2),[],NaN('single')),1:n_domains,'uni',false),[1,3,2]));

clearvars -except load_dataset_retain data_grid_params data_grids

% Get passive activity
load_dataset = 'passive';
Marica_2025.figures.load_data;

% (striatum)
striatum_ld_idx = discretize(striatum_mua_grp.ld,data_grid_params.ld_bins);
[~,striatum_stim_idx] = ismember(striatum_mua_grp.stim,unique(striatum_mua_grp.stim));

[striatum_rec,striatum_rec_grp] = ap.groupfun(@nanmean,striatum_mua(:,:,1), ...
    [striatum_mua_grp.animal,striatum_ld_idx,striatum_mua_grp.domain_idx,striatum_stim_idx].* ...
    ap.nanout(isnan(striatum_ld_idx)));
striatum_rec_tmax = max(striatum_rec(:,data_grid_params.striatum_stim_t),[],2);

data_grids.striatum_passive = accumarray(striatum_rec_grp,striatum_rec_tmax,[data_grid_params.grid_size,max(striatum_stim_idx)],[],NaN);

% (widefield)
wf_ld_idx = discretize(wf_grp.ld,data_grid_params.ld_bins);
[~,wf_stim_idx] = ismember(wf_grp.stim,unique(wf_grp.stim));

[wf_roi_rec,wf_roi_rec_grp] = ap.groupfun(@nanmean,wf_striatum_roi(:,:,:,1), ...
    [wf_grp.animal,wf_ld_idx,wf_stim_idx].* ...
    ap.nanout(isnan(wf_ld_idx)));
wf_roi_rec_tmax = permute(max(wf_roi_rec(:,data_grid_params.cortex_stim_t,:),[],2),[1,3,2]);

data_grids.wf_roi_passive = permute(cell2mat(permute(arrayfun(@(domain) accumarray(wf_roi_rec_grp, ...
    wf_roi_rec_tmax(:,domain),[data_grid_params.grid_size(1:2),length(unique(wf_grp.stim))], ...
    [],NaN('single')),1:n_domains,'uni',false),[1,3,4,2])),[1,2,4,3]);

clearvars -except load_dataset_retain data_grid_params data_grids

% ~~ Get widefield stim kernels
load_dataset = 'noact';
Marica_2025.figures.load_data;

U_master = plab.wf.load_master_U;
load(fullfile(data_path,'wf_kernels'));

n_vs = size(wf_kernels.task_kernels{1},1);

wf_grid_idx = [grp2idx(bhv.animal),discretize(bhv.days_from_learning,data_grid_params.ld_bins)];
wf_grid_idx_use = ~any(isnan(wf_grid_idx),2);

% (task)
wf_kernel_roi_task = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi),[3,2,1]), ...
    wf_kernels.task_kernels,'uni',false));

wf_kernel_roi_task_tmax = permute(max(wf_kernel_roi_task,[],2),[1,3,2]);

data_grids.wf_kernel_roi_task = ...
    cell2mat(permute(arrayfun(@(domain) accumarray(wf_grid_idx(wf_grid_idx_use,:), ...
    wf_kernel_roi_task_tmax(wf_grid_idx_use,domain),data_grid_params.grid_size(1:2),@nanmean,NaN('single')), ...
    1:n_domains,'uni',false),[1,3,2]));

% (passive)
wf_kernel_roi_passive = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi),[4,2,1,3]), ...
    wf_kernels.passive_kernels,'uni',false));
wf_kernel_roi_passive_tmax = permute(max(wf_kernel_roi_passive,[],2),[1,3,4,2]);

data_grids.wf_kernel_roi_passive = cell2mat(permute(arrayfun(@(stim) ...
    cell2mat(permute(arrayfun(@(domain) accumarray(wf_grid_idx(wf_grid_idx_use,:), ...
    wf_kernel_roi_passive_tmax(wf_grid_idx_use,domain,stim),data_grid_params.grid_size(1:2),@nanmean,NaN('single')), ...
    1:n_domains,'uni',false),[1,3,2])),1:size(wf_kernel_roi_passive_tmax,3),'uni',false), [1,3,4,2]));

clearvars -except load_dataset_retain data_grid_params data_grids


% ~~ Plot

% Normalize to task LD0
str_normval = data_grids.striatum_task(:,find(data_grid_params.ld_bins>=0,1),:);
wf_normval = data_grids.wf_roi_task(:,find(data_grid_params.ld_bins>=0,1),:);
wf_kernel_normval = data_grids.wf_kernel_roi_task(:,find(data_grid_params.ld_bins>=0,1),:);

% Set days to plot (n>3)
% plot_ld_idx = sum(~isnan(data_grids.striatum_task(:,:,1))) > 3;
plot_ld_idx = 1:size(data_grids.striatum_task,2);

% Plot task vs passive for mPFC and striatum
max_ld = max(abs(data_grid_params.ld_unique(plot_ld_idx)));
ld_colors = ap.colormap('BKR',max_ld+(1-mod(max_ld,2)));
% plot_ld_colors = ld_colors(ismember(-max_ld:max_ld, ...
%     data_grid_params.ld_unique(plot_ld_idx)),:);
plot_ld_colors = ld_colors(1:max_ld,:);

plot_str = 1;
plot_wf_roi = 2;
plot_stim = 3;

figure; tiledlayout(1,2);
nexttile; hold on;
scatter(reshape(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,plot_wf_roi,plot_stim),[],1), ...
    reshape(data_grids.wf_kernel_roi_task(:,plot_ld_idx,plot_wf_roi),[],1), ...
    20,repelem(plot_ld_colors,size(data_grids.striatum_task,1),1),'filled');
plot(nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,plot_wf_roi,plot_stim),1), ...
    nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,plot_wf_roi),1),'color',[0.5,0.5,0.5],'linewidth',2);
scatter(nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,plot_wf_roi,plot_stim),1), ...
    nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,plot_wf_roi),1),80,plot_ld_colors,'linewidth',4);
xlim(prctile([xlim,ylim],[0,100]));ylim(xlim)
axis square
line(ylim,ylim,'linestyle','--','color',[0.5,0.5,0.5]);
xlabel('Passive');ylabel('Task');
title(sprintf('Cortex kernel %d',plot_wf_roi));

nexttile; hold on;
scatter(reshape(data_grids.striatum_passive(:,plot_ld_idx,plot_str,plot_stim),[],1), ...
    reshape(data_grids.striatum_task(:,plot_ld_idx,plot_str),[],1), ...
    20,repelem(plot_ld_colors,size(data_grids.striatum_task,1),1),'filled');
plot(nanmean(data_grids.striatum_passive(:,plot_ld_idx,plot_str,plot_stim),1), ...
    nanmean(data_grids.striatum_task(:,plot_ld_idx,plot_str),1),'color',[0.5,0.5,0.5],'linewidth',2);
scatter(nanmean(data_grids.striatum_passive(:,plot_ld_idx,plot_str,plot_stim),1), ...
    nanmean(data_grids.striatum_task(:,plot_ld_idx,plot_str),1),80,plot_ld_colors,'linewidth',4);
xlim(prctile([xlim,ylim],[0,100]));ylim(xlim)
axis square
line(ylim,ylim,'linestyle','--','color',[0.5,0.5,0.5]);
xlabel('Passive');ylabel('Task');
title(sprintf('Striatum %d',plot_str));
ap.prettyfig;

% Plot striatum vs mPFC for task and passive
outline_ld_cols = [0.5,0.5,0.5;0,0,0];
outline_ld = @(h) scatter(h.XData(find(data_grid_params.ld_bins==0)+[-1,0]), ...
    h.YData(find(data_grid_params.ld_bins==0)+[-1,0]), ...
    100,outline_ld_cols,'linewidth',3);

% (wf kernel, LD-0 norm)
plot_task = @(roi,str,col) ...
    plot(nanmean(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
    nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,roi)./wf_kernel_normval(:,:,roi),1), ...
    '.-','color',col,'linewidth',2,'MarkerSize',30);
plot_passive = @(roi,str,plot_stim,col) ...
    plot(nanmean(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
    nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_kernel_normval(:,:,roi),1), ...
    '.-','color',col,'linewidth',2,'MarkerSize',30);

% % (wf kernel, LD-0 norm,errorbar)
% plot_task = @(roi,str,col) ...
%     errorbar(nanmean(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
%     nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,roi)./wf_kernel_normval(:,:,roi),1), ...
%     ap.sem(data_grids.wf_kernel_roi_task(:,plot_ld_idx,roi)./wf_kernel_normval(:,:,roi),1), ...
%     ap.sem(data_grids.wf_kernel_roi_task(:,plot_ld_idx,roi)./wf_kernel_normval(:,:,roi),1), ...
%     ap.sem(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
%     ap.sem(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30,'capsize',0);
% plot_passive = @(roi,str,plot_stim,col) ...
%     errorbar(nanmean(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
%     nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_kernel_normval(:,:,roi),1), ...
%     ap.sem(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_kernel_normval(:,:,roi),1), ...
%     ap.sem(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_kernel_normval(:,:,roi),1), ...
%     ap.sem(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
%     ap.sem(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30,'capsize',0);

% % (wf dff, non-norm)
% plot_task = @(roi,str,col) ...
%     plot(nanmean(data_grids.striatum_task(:,plot_ld_idx,str),1), ...
%     nanmean(data_grids.wf_roi_task(:,plot_ld_idx,roi),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30);
% plot_passive = @(roi,str,plot_stim,col) ...
%     plot(nanmean(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim),1), ...
%     nanmean(data_grids.wf_roi_passive(:,plot_ld_idx,roi,plot_stim),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30);

% % (wf dff, norm)
% plot_task = @(roi,str,col) ...
%     plot(nanmean(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
%     nanmean(data_grids.wf_roi_task(:,plot_ld_idx,roi)./wf_normval(:,:,roi),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30);
% plot_passive = @(roi,str,plot_stim,col) ...
%     plot(nanmean(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
%     nanmean(data_grids.wf_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_normval(:,:,roi),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30);

plot_str = 1;
plot_wf_roi = 2;

figure; hold on;
h = plot_task(plot_wf_roi,plot_str,[0.5,0,0]); outline_ld(h);
stim_col = [0.8,0.8,0.3;0.3,0.3,0.8;0.8,0.3,0.3];
for curr_stim = 1:3
    h = plot_passive(plot_wf_roi,plot_str,curr_stim,stim_col(curr_stim,:)); outline_ld(h);
end
xlabel(sprintf('Striatum %d (LD0-norm)',plot_str));
ylabel(sprintf('Cortex kernel %d (LD0-norm)',plot_wf_roi));
ap.prettyfig;


figure; hold on;
plot_str = [1,2];
plot_wf_roi = 2;
h1 = plot_task(plot_wf_roi,plot_str(1),[0.5,0,0]);
h2 = plot_task(plot_wf_roi,plot_str(2),'k');
outline_ld(h1); outline_ld(h2);
axis equal square
line(xlim,xlim,'color',[0.5,0.5,0.5]);
legend(arrayfun(@(x) sprintf('Str %d',x),plot_str,'uni',false));
ylabel(sprintf('Cortex kernel %d (LD0-norm)',plot_wf_roi));
xlabel('Striatum (LD0-norm)')
ap.prettyfig;


% ~~~ STATS ~~~
n_shuff = 10000;
sig_flag = @(p) discretize(p < 0.05,[0,1,Inf],["","*"]);

plot_str = 1;
plot_wf_roi = 2;
plot_stim = 3;

pre_days =  data_grid_params.ld_unique == 1;
post_days = ismember(data_grid_params.ld_unique,3);

stat_label = {...
    sprintf('striatum %d task',plot_str), ...
    sprintf('striatum %d passive',plot_str), ...
    sprintf('cortex kernel %d task',plot_wf_roi), ...
    sprintf('cortex_kernel %d passive',plot_wf_roi)};
pre_data = [...
    nanmean(data_grids.striatum_task(:,pre_days,plot_str)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.striatum_passive(:,pre_days,plot_str,plot_stim)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.wf_kernel_roi_task(:,pre_days,plot_wf_roi)./wf_kernel_normval(:,:,plot_wf_roi),2), ...
    nanmean(data_grids.wf_kernel_roi_passive(:,pre_days,plot_wf_roi,plot_stim)./wf_kernel_normval(:,:,plot_wf_roi),2)];
post_data = [...
    nanmean(data_grids.striatum_task(:,post_days,plot_str)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.striatum_passive(:,post_days,plot_str,plot_stim)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.wf_kernel_roi_task(:,post_days,plot_wf_roi)./wf_kernel_normval(:,:,plot_wf_roi),2), ...
    nanmean(data_grids.wf_kernel_roi_passive(:,post_days,plot_wf_roi,plot_stim)./wf_kernel_normval(:,:,plot_wf_roi),2)];

data_meas = nanmean(diff(cat(3,pre_data,post_data),[],3),1);
data_shuff = nan(n_shuff,size(data_meas,2));
for curr_shuff = 1:n_shuff
    data_shuff(curr_shuff,:) = nanmean(diff(ap.shake(cat(3,pre_data,post_data),3),[],3),1);
end
stat_rank = tiedrank([data_meas;data_shuff]);
stat_p = 1-stat_rank(1,:)/(n_shuff+1);
for curr_stat = 1:length(stat_p)
    fprintf('Shuffle: %s p = %.2g%s\n', ...
        stat_label{curr_stat},stat_p(curr_stat),sig_flag(stat_p(curr_stat)));
end

for curr_stat = 1:size(pre_data,2)
    stat_signrank = ranksum(pre_data(:,curr_stat),post_data(:,curr_stat));
    fprintf('Ranksum: %s p = %.2g%s\n', ...
        stat_label{curr_stat},stat_signrank,sig_flag(stat_signrank));
end


%% R2 p4: CCF-aligned probe histology

[av,tv,st] = ap_histology.load_ccf;

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

% Set atlas bins to plot through striatum
striatum_ccf_id = find(contains(lower(st.safe_name),'caudoputamen'));
striatum_ccf_ap = prctile(find(max(av == striatum_ccf_id,[],[2,3])),[0,100]);

n_atlas_bins = 10;
atlas_bins = round(linspace(striatum_ccf_ap(1),striatum_ccf_ap(2),n_atlas_bins+1));

histology_atlas_bin_max = cell(length(animals),n_atlas_bins);

for curr_animal_idx = 1:length(animals)

    histology_path = plab.locations.filename('server',animals{curr_animal_idx},[],[],'histology');
    load(fullfile(histology_path,'AP_histology_processing.mat'));

    % Load images
    image_path = histology_path;
    image_dir = dir(fullfile(image_path,'*.tif'));
    image_filenames = cellfun(@(path,name) fullfile(path,name), ...
        {image_dir.folder},{image_dir.name},'uni',false);
    [~,sort_idx] = ap_histology.natsortfiles(image_filenames);

    images = cell(size(image_dir));
    for curr_im = 1:length(sort_idx)
        images{curr_im} = tiffreadVolume( ...
            image_filenames{sort_idx(curr_im)});
    end

    % Grab atlas images
    n_slices = length(images);
    slice_atlas = struct('tv',cell(n_slices,1), 'av',cell(n_slices,1));
    slice_atlas_ccf = struct('ap',cell(n_slices,1),'ml',cell(n_slices,1),'dv',cell(n_slices,1));
    for curr_slice = 1:length(images)
        [slice_atlas(curr_slice),slice_atlas_ccf(curr_slice)] = ...
            ap_histology.grab_atlas_slice(av,tv, ...
            AP_histology_processing.histology_ccf.slice_vector, ...
            AP_histology_processing.histology_ccf.slice_points(curr_slice,:), 1);
    end

    % Build volume of histology images
    histology_volume = zeros(size(tv),'single');
    probe_channel = 2;
    for curr_im_idx = 1:length(images)

        % Rigid transform
        im_rigid_transformed = ap_histology.rigid_transform( ...
            images{curr_im_idx}(:,:,probe_channel),curr_im_idx,AP_histology_processing);

        % Affine/nonlin transform
        if isfield(AP_histology_processing.histology_ccf,'control_points') && ...
                (size(AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx},1) == ...
                size(AP_histology_processing.histology_ccf.control_points.atlas{curr_im_idx},1)) && ...
                size(AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx},1) >= 3
            % Manual alignment (if >3 matched points)
            histology2atlas_tform = fitgeotform2d( ...
                AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx}, ...
                AP_histology_processing.histology_ccf.control_points.atlas{curr_im_idx},'pwl');
        elseif isfield(AP_histology_processing.histology_ccf,'atlas2histology_tform')
            % Automatic alignment
            histology2atlas_tform = invert(AP_histology_processing.histology_ccf.atlas2histology_tform{curr_im_idx});
        end

        atlas_slice_aligned = imwarp(im_rigid_transformed, ...
            histology2atlas_tform,'nearest','OutputView', ...
            imref2d(size(slice_atlas(curr_im_idx).av)));

        % % Check match
        % figure; imshowpair(slice_atlas(curr_im_idx).av,atlas_slice_aligned);

        % Add points to volume in CCF space
        curr_ccf_idx = sub2ind(size(tv), ...
            round(slice_atlas_ccf(curr_im_idx).ap(:)), ...
            round(slice_atlas_ccf(curr_im_idx).dv(:)), ...
            round(slice_atlas_ccf(curr_im_idx).ml(:)));

        histology_volume(curr_ccf_idx) = histology_volume(curr_ccf_idx) + ...
            single(atlas_slice_aligned(:));

    end

    % Get max of histology volume in atlas bins
    for curr_atlas_bin = 1:n_atlas_bins
        plot_atlas_ap = round(mean(atlas_bins(curr_atlas_bin+[0,1])));
        histology_atlas_bin_max{curr_animal_idx,curr_atlas_bin} = ...
            permute(max(histology_volume(atlas_bins(curr_atlas_bin): ...
            atlas_bins(curr_atlas_bin+1),:,:),[],1),[2,3,1]);
    end
    
    ap.print_progress_fraction(curr_animal_idx,length(animals));
end

% Plot probe channel overlay by colored animals
animal_colors = ap.colormap('tube',14);
overlay_dilation = 1;
histology_clim = repelem({[200,500]},length(animals),1);

figure; tiledlayout('TileSpacing','none');
for curr_atlas_bin = 1:n_atlas_bins

    % Plot CCF borders from the middle of the bin
    plot_atlas_ap = round(mean(atlas_bins(curr_atlas_bin+[0,1])));
    curr_ccf_borders = imdilate(boundarymask(permute(av(plot_atlas_ap,:,:),[2,3,1])),ones(overlay_dilation));

    % Max, color, and flip contrast
    curr_histology_volume_max_gray = cellfun(@(x,c) ...
        mat2gray(x,c),histology_atlas_bin_max(:,curr_atlas_bin),histology_clim,'uni',false);
        
    curr_histology_volume_max_colored = ...
        cat(4,curr_histology_volume_max_gray{:}).*permute(animal_colors,[3,4,2,1]);
    
    curr_histology_volume_max_colored_white = ...
        curr_histology_volume_max_colored+(1-cat(4,curr_histology_volume_max_gray{:}));

    curr_histology_combined = min(curr_histology_volume_max_colored_white,[],4);

    % Plot CCF over combined colored image
    curr_overlay = imoverlay(curr_histology_combined,curr_ccf_borders,'k');
    nexttile;imagesc(curr_overlay);axis image off;
    drawnow;
end

%% R3 M3: Fraction of quiescent passive trials

%%% Load non-activity data
load_dataset = 'noact';
Marica_2025.figures.load_data;
%%%

% Get which trials in passive were quiescent
animals = unique(bhv.animal,'stable');

passive_quiescent_trials = struct( ...
    'stim_x',cell(length(animals),1), ...
    'quiescent',cell(length(animals),1));

for animal_idx=1:length(animals)
   
    animal = animals{animal_idx};

    % Find passive recording days that also have task
    workflow_passive = {'lcr_passive'};
    recordings_passive = plab.find_recordings(animal, [], workflow_passive);
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);
    training_days = ismember({recordings_passive.day}, {recordings_task.day});
    train_rec_passive = recordings_passive(training_days);
    bhv_days = {train_rec_passive.day};
    ephys_days =  bhv_days([train_rec_passive.ephys]);    

    for use_rec = 1:length(ephys_days)

        rec_day = train_rec_passive(use_rec).day;
        rec_time = train_rec_passive(use_rec).recording{end};
        load_parts.behavior = true;
        ap.load_recording
   
        % Trial stim values
        trial_stim_values = vertcat(trial_events.values.TrialStimX);
        trial_stim_values = trial_stim_values(1:length(stimOn_times));

        % Quiescent trials
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');

        passive_quiescent_trials(animal_idx).stim_x{use_rec,1} = trial_stim_values;
        passive_quiescent_trials(animal_idx).quiescent{use_rec,1} = quiescent_trials;

    end
    ap.print_progress_fraction(animal_idx,length(animals));
end

stim_unique = unique(passive_quiescent_trials(1).stim_x{1});
stim_color = ap.colormap('BKR',length(stim_unique));

% Plot quiescent trials pre/post learning
passive_quiescent_stim = cell2mat(cellfun(@(stim,quiescent) ...
    ap.groupfun(@mean,+quiescent,stim)', ...
    cat(1,passive_quiescent_trials.stim_x), ...
    cat(1,passive_quiescent_trials.quiescent),'uni',false));

[passive_quiescent_stim_ld,passive_quiescent_stim_ld_grp] = ...
    ap.nestgroupfun({@mean,@mean},passive_quiescent_stim, ...
    grp2idx(bhv.animal),bhv.days_from_learning >= 0);
passive_quiescent_stim_ld_sem = ...
    ap.nestgroupfun({@mean,@ap.sem},passive_quiescent_stim, ...
    grp2idx(bhv.animal),bhv.days_from_learning >= 0);

figure;
x_labels = ["Pre-learn","Post-learn"];
errorbar(reordercats(categorical(x_labels),x_labels), ...
    passive_quiescent_stim_ld,passive_quiescent_stim_ld_sem, ...
    'linewidth',2,'capsize',0)
axis padded;
ylabel('Frac. quiescent trials');
set(gca,'ColorOrder',stim_color)
ap.prettyfig;

% Plot quiescent trials binned by trial number
n_trialsplit = length(stim_unique)*5; % n presentations of each stim
passive_trialsplit_idx = cell2mat(cellfun(@(stim) ...
    (floor((0:length(stim)-1)/n_trialsplit)+1)', ...
    cat(1,passive_quiescent_trials.stim_x),'uni',false));

passive_ld_idx = cell2mat(cellfun(@(rec,stim) repelem(rec,length(stim))', ...
    num2cell(bhv.days_from_learning), ...
    cat(1,passive_quiescent_trials.stim_x),'uni',false));

[passive_quiescent_stim_trialsplit,passive_quiescent_stim_trialsplit_grp] = ...
    ap.groupfun(@mean,+cell2mat(cat(1,passive_quiescent_trials.quiescent)), ...
    [passive_ld_idx,passive_trialsplit_idx,cell2mat(cat(1,passive_quiescent_trials.stim_x))]);

passive_quiescent_stim_trialsplit_sem = ...
    ap.groupfun(@ap.sem,+cell2mat(cat(1,passive_quiescent_trials.quiescent)), ...
    [passive_ld_idx,passive_trialsplit_idx,cell2mat(cat(1,passive_quiescent_trials.stim_x))]);

passive_quiescent_stim_trialsplit_grid = ...
    accumarray([grp2idx(passive_quiescent_stim_trialsplit_grp(:,1)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,2)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,3))], ...
    passive_quiescent_stim_trialsplit,[],[],nan);

passive_quiescent_stim_trialsplit_sem_grid = ...
    accumarray([grp2idx(passive_quiescent_stim_trialsplit_grp(:,1)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,2)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,3))], ...
    passive_quiescent_stim_trialsplit_sem,[],[],nan);

passive_quiescent_stim_daysplit_x = reshape((unique(passive_quiescent_stim_trialsplit_grp(:,1)) + ...
    linspace(0,1,max(passive_quiescent_stim_trialsplit_grp(:,2)+1)))',[],1);

figure;
plot_x_idx = isbetween(passive_quiescent_stim_daysplit_x,-3,3);

passive_quiescent_stim_trialsplit_grid_reshape = ...
    reshape(permute(padarray(passive_quiescent_stim_trialsplit_grid,[0,1],nan,'post'),[2,1,3]),[],3);
passive_quiescent_stim_trialsplit_sem_grid_reshape = ...
    reshape(permute(padarray(passive_quiescent_stim_trialsplit_sem_grid,[0,1],nan,'post'),[2,1,3]),[],3);

ap.errorfill(passive_quiescent_stim_daysplit_x(plot_x_idx), ...
    passive_quiescent_stim_trialsplit_grid_reshape(plot_x_idx,:), ...
    passive_quiescent_stim_trialsplit_sem_grid_reshape(plot_x_idx,:),stim_color);
xline(0,'k');
xlabel('Days from learning');
ylabel('Frac. quiescent trials');
axis padded;
ap.prettyfig;

% ~~~ STATS ~~~
sig_flag = @(p) discretize(p < 0.05,[0,1,Inf],["","*"]);

[passive_quiescent_stim_ld_animal,passive_quiescent_stim_ld_animal_grp] = ...
    ap.groupfun(@mean,passive_quiescent_stim, ...
    [grp2idx(bhv.animal),bhv.days_from_learning >= 0]);

for curr_stim = 1:length(stim_unique)
    stat_p = ranksum( ...
        passive_quiescent_stim_ld_animal(passive_quiescent_stim_ld_animal_grp(:,2) == 0,curr_stim), ...
        passive_quiescent_stim_ld_animal(passive_quiescent_stim_ld_animal_grp(:,2) == 1,curr_stim));

    print_stat('Ranksum pre/post stim %d: p = %.2g%s\n', ...
        stim_unique(curr_stim),stat_p,sig_flag(stat_p));
end

%% R3m1: Rate of non-stim movements

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

use_stat = 'mean';

% Grab learning day for each mouse
bhv = struct;

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};

    use_workflow = 'stim_wheel_right_stage\d';
    recordings = plab.find_recordings(animal,[],use_workflow);

    n_bins = 3;

    move_rate = nan(length(recordings),n_bins);
    move_rate_stim = nan(length(recordings),n_bins);
    move_rate_nonstim = nan(length(recordings),n_bins);
    p_cue_given_move = nan(length(recordings),n_bins);
    rxn_stat_p = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get P(cue|wheel_move)
        move_onsets = timelite.timestamps(diff(wheel_move) == 1);
        move_onset_prev_stim = interp1(stimOn_times,stimOn_times,move_onsets,'previous','extrap');
        move_prevstim_t = move_onsets - move_onset_prev_stim;

        move_poststim = move_prevstim_t < 0.4;

        % Bin by time within session
        session_bin_edges = linspace(timelite.timestamps(1),timelite.timestamps(end),n_bins+1);
        move_session_bins = discretize(move_onsets,session_bin_edges);

        move_rate(curr_recording,:) = ap.groupfun(@length,+move_poststim,move_session_bins)./diff(session_bin_edges)';
        move_rate_stim(curr_recording,:) = ap.groupfun(@sum,+move_poststim,move_session_bins)./diff(session_bin_edges)';
        move_rate_nonstim(curr_recording,:) = ap.groupfun(@sum,+~move_poststim,move_session_bins)./diff(session_bin_edges)';
        p_cue_given_move(curr_recording,:) = ap.groupfun(@mean,+move_poststim,move_session_bins);

        % Get association stat
        [rxn_stat_p(curr_recording), ...
            rxn_stat(curr_recording),rxn_null_stat(curr_recording)] = ...
            AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    bhv(curr_animal_idx).move_rate = move_rate;
    bhv(curr_animal_idx).move_rate_stim = move_rate_stim;
    bhv(curr_animal_idx).move_rate_nonstim = move_rate_nonstim;
    bhv(curr_animal_idx).p_cue_given_move = p_cue_given_move;
    bhv(curr_animal_idx).rxn_stat_p = rxn_stat_p;

    % Clear vars except pre-load for next loop
    clearvars('-except',preload_vars{:});
    AP_print_progress_fraction(curr_recording,length(recordings));

end

ld = cellfun(@(x) ((1:size(x,1)) - find(x<0.05,1))',{bhv.rxn_stat_p}','uni',false);
use_animals = ~cellfun(@isempty,ld);

move_rate_cat = cell2mat({bhv(use_animals).move_rate}');
move_rate_stim_cat = cell2mat({bhv(use_animals).move_rate_stim}');
move_rate_nonstim_cat = cell2mat({bhv(use_animals).move_rate_nonstim}');
p_c2m_cat = cell2mat({bhv(use_animals).p_cue_given_move}');

ld_cat = cell2mat(ld(use_animals));
ld_split = ld_cat + linspace(0,(n_bins-1)/n_bins,n_bins);

[move_rate_avg,move_rate_grp] = ap.groupfun(@mean,move_rate_cat(:),ld_split(:));
move_rate_sem = ap.groupfun(@AP_sem,move_rate_cat(:),ld_split(:));

[move_rate_stim_avg,move_rate_stim_grp] = ap.groupfun(@mean,move_rate_stim_cat(:),ld_split(:));
move_rate_stim_sem = ap.groupfun(@AP_sem,move_rate_stim_cat(:),ld_split(:));

[move_rate_nonstim_avg,move_rate_nonstim_grp] = ap.groupfun(@mean,move_rate_nonstim_cat(:),ld_split(:));
move_rate_nonstim_sem = ap.groupfun(@AP_sem,move_rate_nonstim_cat(:),ld_split(:));

[p_c2m_avg,p_c2m_avg_grp] = ap.groupfun(@mean,p_c2m_cat(:),ld_split(:));
p_c2m_sem = ap.groupfun(@AP_sem,p_c2m_cat(:),ld_split(:));

figure('Name','Fig S1 p stim move'); tiledlayout(2,1);

ax1 = nexttile; hold on;
h1 = errorbar( ...
    padarray(reshape(move_rate_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_sem,n_bins,[]),[1,0],NaN,'post'),'color',[0.5,0.5,0.5],'linewidth',2);
h2 = errorbar( ...
    padarray(reshape(move_rate_stim_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_stim_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_stim_sem,n_bins,[]),[1,0],NaN,'post'),'color',[0,0.5,0],'linewidth',2);
h3 = errorbar( ...
    padarray(reshape(move_rate_nonstim_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_nonstim_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_nonstim_sem,n_bins,[]),[1,0],NaN,'post'),'color',[0.5,0,0],'linewidth',2);
xlabel('Learned day')
ylabel(' Moves/s');
xline(0,'r');
legend([h1(1),h2(1),h3(1)],{'All','Non-stim','Stim'});

ax2 = nexttile;
errorbar( ...
    padarray(reshape(p_c2m_avg_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(p_c2m_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(p_c2m_sem,n_bins,[]),[1,0],NaN,'post'),'k','linewidth',2);
xlabel('Learned day')
ylabel('P(stim|move)');
xline(0,'r');

linkaxes([ax1,ax2],'x');
ap.prettyfig;
xlim([-3.1,2.9])


%% R3m3: Non-stim move activity

%%% Load non-activity data
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

% Load nonstim move activity
U_master = plab.wf.load_master_U;
load(fullfile(data_path,'nonstim_move'));

use_nostim_move_recordings = ~cellfun(@isempty,nonstim_move.V_move_nostim_align);
nostim_move_wheel_t = nonstim_move.wheel_align_time{find(use_nostim_move_recordings,1)};

%%% NON-STIM MOVE DATA PREPROCESSING
% Get widefield ROIs for no stim moves
wf_nostim_move_striatum_roi = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master,x,[],[],striatum_wf_roi), ...
    [3,2,1]),nonstim_move.V_move_nostim_align(use_nostim_move_recordings), ...
    'uni',false));
baseline_t = wf_t < 0;
wf_nostim_move_striatum_roi = wf_nostim_move_striatum_roi - nanmean(wf_nostim_move_striatum_roi(:,baseline_t,:,:),2);

% Get nonstim move ephys
% (sum into domain multiunit)
striatum_nostim_move_mua_sum = cellfun(@(mua,domain_idx) ...
    permute(ap.groupfun(@sum,mua,domain_idx,[]),[3,2,1]), ...
    nonstim_move.binned_msn_spikes_move_nostim_align,domain_idx_rec,'uni',false);
% (smooth and normalize - sub baseline = trial, div baseline = day) 
baseline_t = psth_t < 0;
striatum_nostim_sub_baseline = cellfun(@(mua) ...
        mean(mua(:,baseline_t,:,1),[2]), ...
        striatum_nostim_move_mua_sum,'uni',false,'ErrorHandler',@(varargin) NaN);
striatum_nostim_move_mua = cellfun(mua_norm_smooth_reshape_fcn, ...
    striatum_nostim_move_mua_sum,cellfun(@(x) nanmean(x,1), ...
    striatum_nostim_sub_baseline,'uni',false),mua_div_baseline,'uni',false, ...
    'ErrorHandler',@(varargin) []); %#ok<FUNFUN>
%%%

% Plot move-aligned (stim, non-stim) binned by day
% (non-stim not saved by trial, so day bins by weighted average)
plot_day_bins = [-Inf,0,Inf];

plot_day_colors = [0,0,0;0.7,0,0];

day_grp = discretize(max(bhv.days_from_learning,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

figure; h_striatum = tiledlayout(n_domains,3); title(h_striatum,'Striatum');
figure; h_cortex = tiledlayout(n_domains,3); title(h_cortex,'Cortex');
figure; h_rxn = axes; hold on; set(gca,'ColorOrder',plot_day_colors);

figure; h_wheel = tiledlayout(1,2);
nexttile(h_wheel); hold on; set(gca,'ColorOrder',plot_day_colors);
nexttile(h_wheel); hold on; set(gca,'ColorOrder',plot_day_colors);

for curr_domain = 1:n_domains

    h_striatum_sub = gobjects(3,1);
    h_cortex_sub = gobjects(3,1);
    for sub = 1:3
        h_striatum_sub(sub) = nexttile(h_striatum); hold on;
        set(gca,'ColorOrder',plot_day_colors);

        h_cortex_sub(sub) = nexttile(h_cortex); hold on;
        set(gca,'ColorOrder',plot_day_colors);
    end

    for curr_day_grp = 1:length(plot_day_bins)-1

        % Get day-binned (weighted average) activity
        curr_nomove_idx = use_nostim_move_recordings & day_grp==curr_day_grp;
        n_moves = cellfun(@(x) size(x,1),nonstim_move.move_nostim_wheel(curr_nomove_idx));

        animals = unique(bhv.animal);
        [~,curr_nostim_move_animal_idx] = ismember(nonstim_move.animal(curr_nomove_idx),animals);

        % (striatum)
        curr_nostim_striatum = cellfun(@(x,domain) x(domain==curr_domain,:), ...
            striatum_nostim_move_mua(curr_nomove_idx), ...
            striatum_domain_idx(curr_nomove_idx),'uni',false);
        
        curr_striatum_nostim_weighted = cellfun(@(act,n) act.*n,curr_nostim_striatum,num2cell(n_moves),'uni',false);
        
        curr_striatum_nostim_wavg = ap.groupfun(@sum,cell2mat(curr_striatum_nostim_weighted), ...
            curr_nostim_move_animal_idx(~cellfun(@isempty,curr_striatum_nostim_weighted)))./ ...
            ap.groupfun(@sum,n_moves(~cellfun(@isempty,curr_striatum_nostim_weighted)), ...
             curr_nostim_move_animal_idx(~cellfun(@isempty,curr_striatum_nostim_weighted)));

        % (cortex)
        curr_nostim_cortex = wf_nostim_move_striatum_roi( ...
            day_grp(use_nostim_move_recordings)==curr_day_grp,:,curr_domain);

        curr_nostim_cortex_wavg = ...
            ap.groupfun(@sum,curr_nostim_cortex.*n_moves,curr_nostim_move_animal_idx)./ ...
            ap.groupfun(@sum,n_moves,curr_nostim_move_animal_idx);
        
        % Get stim move activity 
        % (only use animals from above to pair data)
        
        % (striatum)
        curr_nomove_animals_striatum = unique(curr_nostim_move_animal_idx(~cellfun(@isempty,curr_striatum_nostim_weighted)));
        curr_striatum_trials_idx = ismember(striatum_mua_grp.animal,curr_nomove_animals_striatum) & ...
            striatum_mua_grp.domain_idx == curr_domain & ...
            striatum_plot_day_grp == curr_day_grp;

        curr_striatum_stimmove = ap.groupfun(@mean, ...
            striatum_mua(curr_striatum_trials_idx,:,2), ...
            striatum_mua_grp.animal(curr_striatum_trials_idx));

        % (cortex) 
        curr_nomove_animals_cortex = unique(curr_nostim_move_animal_idx);
        curr_cortex_trials_idx = ismember(wf_grp.animal,curr_nomove_animals_cortex) & ...
            cortex_plot_day_grp == curr_day_grp;
       
        curr_cortex_stimmove = ap.groupfun(@mean, ...
            wf_striatum_roi(curr_cortex_trials_idx,:,curr_domain,2), ...
            wf_grp.animal(curr_cortex_trials_idx));

        % Plot
        axes(h_striatum_sub(1));ap.errorfill(psth_t,nanmean(curr_striatum_stimmove,1),ap.sem(curr_striatum_stimmove,1));
        axes(h_striatum_sub(2));ap.errorfill(psth_t,nanmean(curr_striatum_nostim_wavg,1),ap.sem(curr_striatum_nostim_wavg,1));
        axes(h_striatum_sub(3));ap.errorfill(psth_t,nanmean(curr_striatum_stimmove-curr_striatum_nostim_wavg,1), ...
            ap.sem(curr_striatum_stimmove-curr_striatum_nostim_wavg,1));

        axes(h_cortex_sub(1));ap.errorfill(wf_t,nanmean(curr_cortex_stimmove,1),ap.sem(curr_cortex_stimmove,1));
        axes(h_cortex_sub(2));ap.errorfill(wf_t,nanmean(curr_nostim_cortex_wavg,1),ap.sem(curr_nostim_cortex_wavg,1));
        axes(h_cortex_sub(3));ap.errorfill(wf_t,nanmean(curr_cortex_stimmove-curr_nostim_cortex_wavg,1), ...
            ap.sem(curr_cortex_stimmove-curr_nostim_cortex_wavg,1));

        if curr_domain == 1
            % Plot histogram of stim relative to move onset (on first domain)
            rxn_bins = [-Inf,-0.5:0.05:1,Inf];
            rxn_bin_x = [rxn_bins(2),rxn_bins(2:end-2)+diff(rxn_bins(2:end-1))/2,rxn_bins(end-1)];
            rxn_histogram = cell2mat(arrayfun(@(x) histcounts(-striatum_mua_grp.rxn(curr_striatum_trials_idx & ...
                striatum_mua_grp.animal == x),rxn_bins,'Normalization','probability'), ...
                (1:length(unique(striatum_mua_grp.animal)))','uni',false));
            axes(h_rxn);ap.errorfill(rxn_bin_x,nanmean(rxn_histogram,1),ap.sem(rxn_histogram,1));

            % Plot wheel
            wheel_stim_animal = arrayfun(@(x) nanmean(cell2mat(nonstim_move.move_stim_wheel( ...
                (grp2idx(bhv.animal) == x) & (day_grp == curr_day_grp))),1), ...
                unique(grp2idx(bhv.animal)),'uni',false);  
            nexttile(h_wheel,1);
            ap.errorfill(nonstim_move.wheel_align_time{end}, ...
                nanmean(vertcat(wheel_stim_animal{:}),1),ap.sem(vertcat(wheel_stim_animal{:}),1));      

            wheel_nostim_animal = arrayfun(@(x) nanmean(cell2mat(nonstim_move.move_nostim_wheel( ...
                (grp2idx(bhv.animal) == x) & (day_grp == curr_day_grp))),1), ...
                unique(grp2idx(bhv.animal)),'uni',false);        
            nexttile(h_wheel,2);
            ap.errorfill(nonstim_move.wheel_align_time{end}, ...
                nanmean(vertcat(wheel_nostim_animal{:}),1),ap.sem(vertcat(wheel_nostim_animal{:}),1));        
        end

    end
end
linkaxes(h_striatum.Children,'xy');
linkaxes(h_cortex.Children,'xy');
linkaxes(h_wheel.Children,'xy');

ap.prettyfig([],h_striatum.Parent);
ap.prettyfig([],h_cortex.Parent);
ap.prettyfig([],h_wheel.Parent);

