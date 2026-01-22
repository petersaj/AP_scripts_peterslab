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
ylabel('Association index');
xlabel('Day from learning');
ap.prettyfig;


%% R2 p4: histological alignment

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

histology_clim = arrayfun(@(x) ...
    double(max(cat(3,histology_atlas_bin_max{x,:}),[],'all')./[8,5]), ...
    (1:length(animals))','uni',false);

figure; tiledlayout('TileSpacing','none');
for curr_atlas_bin = 1:n_atlas_bins

    plot_atlas_ap = round(mean(atlas_bins(curr_atlas_bin+[0,1])));
    curr_ccf_borders = imdilate(boundarymask(permute(av(plot_atlas_ap,:,:),[2,3,1])),ones(overlay_dilation));

    curr_histology_volume_max_gray = cellfun(@(x,c) ...
        mat2gray(x,c),histology_atlas_bin_max(:,curr_atlas_bin),histology_clim,'uni',false);
        
    curr_histology_volume_max_colored = ...
        cat(4,curr_histology_volume_max_gray{:}).*permute(animal_colors,[3,4,2,1]);
    
    % (make this colored on white - not sure how to combined though)
    % curr_histology_volume_max_colored_white = ...
    %     curr_histology_volume_max_colored+(1-cat(4,curr_histology_volume_max_gray{:}));

    curr_histology_combined = max(curr_histology_volume_max_colored,[],4);

    curr_overlay = imoverlay(curr_histology_combined,curr_ccf_borders,'w');

    nexttile;imagesc(curr_overlay);axis image off;
    drawnow;
end


%% R3 maj1: Plot striatum vs mPFC in task/passive

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

data_grid_params.ld_unique = unique((bhv.days_from_learning(~isnan(bhv.days_from_learning))));
data_grid_params.grid_size = [length(unique(bhv.animal)),length(data_grid_params.ld_unique),n_domains];

% Set up data structure
data_grids = struct;

clearvars -except load_dataset_retain data_grid_params data_grids

% ~~ Get task activity
load_dataset = 'task';
Marica_2025.figures.load_data;

% (striatum)
striatum_use_trials = striatum_mua_grp.rxn > data_grid_params.rxn_cutoff;
[~,striatum_ld_idx] = ismember(striatum_mua_grp.ld,data_grid_params.ld_unique);

[striatum_rec,striatum_rec_grp] = ap.groupfun(@nanmean,striatum_mua(:,:,1), ...
    [striatum_mua_grp.animal,striatum_ld_idx,striatum_mua_grp.domain_idx].* ...
    ap.nanout(~(striatum_use_trials & striatum_ld_idx)));
striatum_rec_tmax = max(striatum_rec(:,data_grid_params.striatum_stim_t),[],2);

data_grids.striatum_task = accumarray(striatum_rec_grp,striatum_rec_tmax,data_grid_params.grid_size,[],NaN);

% (widefield)
wf_use_trials = wf_grp.rxn > data_grid_params.rxn_cutoff;
[~,wf_ld_idx] = ismember(wf_grp.ld,data_grid_params.ld_unique);

[wf_roi_rec,wf_roi_rec_grp] = ap.groupfun(@nanmean,wf_striatum_roi(:,:,:,1), ...
    [wf_grp.animal,wf_ld_idx].* ...
    ap.nanout(~(wf_use_trials & wf_ld_idx)));
wf_roi_rec_tmax = permute(max(wf_roi_rec(:,data_grid_params.cortex_stim_t,:),[],2),[1,3,2]);

data_grids.wf_roi_task = cell2mat(permute(arrayfun(@(domain) accumarray(wf_roi_rec_grp, ...
    wf_roi_rec_tmax(:,domain),data_grid_params.grid_size(1:2),[],NaN('single')),1:n_domains,'uni',false),[1,3,2]));

clearvars -except load_dataset_retain data_grid_params data_grids

% Get passive activity
load_dataset = 'passive';
Marica_2025.figures.load_data;

% (striatum)
[~,striatum_ld_idx] = ismember(striatum_mua_grp.ld,data_grid_params.ld_unique);
[~,striatum_stim_idx] = ismember(striatum_mua_grp.stim,unique(striatum_mua_grp.stim));

[striatum_rec,striatum_rec_grp] = ap.groupfun(@nanmean,striatum_mua(:,:,1), ...
    [striatum_mua_grp.animal,striatum_ld_idx,striatum_mua_grp.domain_idx,striatum_stim_idx].* ...
    ap.nanout(~(striatum_ld_idx)));
striatum_rec_tmax = max(striatum_rec(:,data_grid_params.striatum_stim_t),[],2);

data_grids.striatum_passive = accumarray(striatum_rec_grp,striatum_rec_tmax,[data_grid_params.grid_size,max(striatum_stim_idx)],[],NaN);

% (widefield)
[~,wf_ld_idx] = ismember(wf_grp.ld,data_grid_params.ld_unique);
[~,wf_stim_idx] = ismember(wf_grp.stim,unique(wf_grp.stim));

[wf_roi_rec,wf_roi_rec_grp] = ap.groupfun(@nanmean,wf_striatum_roi(:,:,:,1), ...
    [wf_grp.animal,wf_ld_idx,wf_stim_idx].* ...
    ap.nanout(~(wf_ld_idx)));
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

wf_grid_idx = [grp2idx(bhv.animal),grp2idx(bhv.days_from_learning)];
wf_grid_idx_use = ~any(isnan(wf_grid_idx),2);

% (task)
wf_kernel_roi_task = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi),[3,2,1]), ...
    wf_kernels.task_kernels,'uni',false));

wf_kernel_roi_task_tmax = permute(max(wf_kernel_roi_task,[],2),[1,3,2]);

data_grids.wf_kernel_roi_task = ...
    cell2mat(permute(arrayfun(@(domain) accumarray(wf_grid_idx(wf_grid_idx_use,:), ...
    wf_kernel_roi_task_tmax(wf_grid_idx_use,domain),data_grid_params.grid_size(1:2),[],NaN('single')), ...
    1:n_domains,'uni',false),[1,3,2]));

% (passive)
wf_kernel_roi_passive = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi),[4,2,1,3]), ...
    wf_kernels.passive_kernels,'uni',false));
wf_kernel_roi_passive_tmax = permute(max(wf_kernel_roi_passive,[],2),[1,3,4,2]);

data_grids.wf_kernel_roi_passive = cell2mat(permute(arrayfun(@(stim) ...
    cell2mat(permute(arrayfun(@(domain) accumarray(wf_grid_idx(wf_grid_idx_use,:), ...
    wf_kernel_roi_passive_tmax(wf_grid_idx_use,domain,stim),data_grid_params.grid_size(1:2),[],NaN('single')), ...
    1:n_domains,'uni',false),[1,3,2])),1:size(wf_kernel_roi_passive_tmax,3),'uni',false), [1,3,4,2]));

clearvars -except load_dataset_retain data_grid_params data_grids


% ~~ Plot

% Normalize to task LD0
str_normval = data_grids.striatum_task(:,data_grid_params.ld_unique==0,:);
wf_normval = data_grids.wf_roi_task(:,data_grid_params.ld_unique==0,:);
wf_kernel_normval = data_grids.wf_kernel_roi_task(:,data_grid_params.ld_unique==0,:);

% Set days to plot (n>3)
plot_ld_idx = sum(~isnan(data_grids.striatum_task(:,:,1))) > 3;

% Plot task vs passive for mPFC and striatum
max_ld = max(abs(data_grid_params.ld_unique(plot_ld_idx)));
ld_colors = ap.colormap('BKR',max_ld*2+1);
plot_ld_colors = ld_colors(ismember(-max_ld:max_ld, ...
    data_grid_params.ld_unique(plot_ld_idx)),:);

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
outline_ld = @(h) scatter(h.XData(ismember(data_grid_params.ld_unique(plot_ld_idx),[-1,0])), ...
    h.YData(ismember(data_grid_params.ld_unique(plot_ld_idx),[-1,0])), ...
    100,outline_ld_cols,'linewidth',3);

plot_task = @(roi,str,col) ...
    plot(nanmean(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
    nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,roi)./wf_kernel_normval(:,:,roi),1), ...
    '.-','color',col,'linewidth',2,'MarkerSize',30);

plot_passive = @(roi,str,plot_stim,col) ...
    plot(nanmean(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
    nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_kernel_normval(:,:,roi),1), ...
    '.-','color',col,'linewidth',2,'MarkerSize',30);

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

% ~~~ STATS ~~~
n_shuff = 10000;
sig_flag = @(p) discretize(p < 0.05,[0,1,Inf],["","*"]);

stat_day = -1; % day to compare vs mean(<day)
plot_str = 1;
plot_wf_roi = 2;
plot_stim = 3;

pre_days =  data_grid_params.ld_unique < -1;
post_days = ismember(data_grid_params.ld_unique,-1);

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
    fprintf('Shuffle day %d vs mean(<): %s p = %.2g%s\n', ...
        stat_day,stat_label{curr_stat},stat_p(curr_stat),sig_flag(stat_p(curr_stat)));
end



