
%% Load ephys

data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');

load(fullfile(data_path,'swr_bhv'));
load(fullfile(data_path,'ctx_maps_to_str'));
load(fullfile(data_path,'ephys'));
U_master = plab.wf.load_master_U;

%% Load widefield

data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');

load(fullfile(data_path,'swr_bhv'));
load(fullfile(data_path,'ctx_wf'));
U_master = plab.wf.load_master_U;

% need to load ephys too because that's where stim info is, change this?
load(fullfile(data_path,'ephys'));


%% K-means on maps 

% K-means cluster maps
maps_cat = cat(3,all_ctx_maps_to_str.cortex_kernel_px{:});
expl_var_cat = vertcat(all_ctx_maps_to_str.explained_var{:});

kmeans_starting = mean(cell2mat(permute(cellfun(@(x) x(:,:,round(linspace(1,size(x,3),4))), ...
    all_ctx_maps_to_str.cortex_kernel_px(~cellfun(@isempty, ...
    all_ctx_maps_to_str.cortex_kernel_px)),'uni',false),[2,3,4,1])),4);

n_k = 4;
[kidx,kmeans_map] = kmeans(...
    reshape(maps_cat,prod(size(U_master,[1,2])),[])',n_k, ...
    'Distance','Correlation','start',reshape(kmeans_starting,[],n_k)');

kmeans_map = reshape(kmeans_map',size(U_master,1),size(U_master,2),[]);

kidx_rec = mat2cell(kidx,cellfun(@(x) size(x,3).*(size(x,1)>0), ...
    all_ctx_maps_to_str.cortex_kernel_px));

figure;
h = tiledlayout(n_k,1,'tilespacing','none');
for curr_k = 1:n_k
    nexttile;
    imagesc(kmeans_map(:,:,curr_k));
    axis image off
    clim(max(abs(clim)).*[-1,1]);
    colormap(ap.colormap('PWG'));
    ap.wf_draw('ccf','k');
end

figure; hold on; axis image; set(gca,'ydir','reverse')
max_map = imgaussfilt(maps_cat,10).*round(linspace(1,0,size(maps_cat,2)));
[~,kernel_max_idx_full] = max(max_map,[],[1,2],'linear');
[kernel_max_y,kernel_max_x,~] = ind2sub(size(maps_cat),kernel_max_idx_full);
for plot_kidx = 1:n_k
    plot(squeeze(kernel_max_x(kidx==plot_kidx)),squeeze(kernel_max_y(kidx==plot_kidx)),'x','linewidth',2.5);
end
ap.wf_draw('ccf','k');
legend(string(1:n_k))



%% Map group by max within ROI

% maps_cat = cat(3,all_ctx_maps_to_str.cortex_kernel_px{:}).^2;
% maps_cat = maps_cat - imgaussfilt(maps_cat,15);

maps_cat = cat(3,all_ctx_maps_to_str.cortex_kernel_px{:});


% Plot sum of maps
figure;
imagesc(mean(max(0,maps_cat),3));
colormap(ap.colormap('WG'));
ap.wf_draw('ccf','k');
axis image; 


% Overlay max weight location
% (exclude right half of image)
max_map = imgaussfilt(maps_cat,5).*round(linspace(1,0,size(maps_cat,2)));
[~,kernel_max_idx_full] = max(max_map,[],[1,2],'linear');
[kernel_max_y,kernel_max_x,~] = ind2sub(size(maps_cat),kernel_max_idx_full);
hold on;
plot(squeeze(kernel_max_x),squeeze(kernel_max_y),'x','LineWidth',2,'color',[0.5,0.5,0.5]);


% Draw ROIs (manually do for N ROIs)
roi = images.roi.Polygon;roi(:) = [];

roi(end+1) = drawpolygon;

roi_masks = arrayfun(@(x) roi(x).createMask,1:length(roi),'uni',false);

% Label masks with maxima within ROIs
roi_maps = cell2mat(cellfun(@(x) squeeze(ismember(sub2ind(size(maps_cat,[1,2]), ...
    kernel_max_y,kernel_max_x),find(x))),roi_masks,'uni',false));

[~,roi_maps_idx] = max(roi_maps,[],2);
kidx = nan(size(roi_maps_idx));
kidx(any(roi_maps,2)) = roi_maps_idx(any(roi_maps,2));

n_k = size(roi_maps,2);
kidx_rec = mat2cell(kidx,cellfun(@(x) size(x,3).*(size(x,1)>0), ...
    all_ctx_maps_to_str.cortex_kernel_px));

% Plot average maps and maxima
maps_cat_mean = ap.groupfun(@mean,maps_cat,[],[],kidx);

maps_com = sum(maps_cat_mean.*permute(1:n_k,[1,3,2]),3)./sum(maps_cat_mean,3);
figure;
im_h = imagesc(maps_com);
ap.wf_draw('ccf','k');
set(im_h,'AlphaData',mat2gray(max(maps_cat_mean,[],3), ...
    [0,double(prctile(reshape(max(maps_cat_mean,[],3),[],1),95))]));
colormap(jet); clim([1,n_k]);
axis image off;

figure; colormap(ap.colormap('KWG'));
h = tiledlayout(n_k,1,'tilespacing','none');
for curr_k = 1:n_k
    nexttile;
    imagesc(maps_cat_mean(:,:,curr_k));
    axis image off
    clim(max(abs(clim)).*[-1,1]); ap.wf_draw('ccf','r');
end

figure; hold on; axis image; set(gca,'ydir','reverse')
plot(squeeze(kernel_max_x(isnan(kidx))),squeeze(kernel_max_y(isnan(kidx))),'x','color',[0.5,0.5,0.5]);
for plot_kidx = 1:n_k
    plot(squeeze(kernel_max_x(kidx==plot_kidx)),squeeze(kernel_max_y(kidx==plot_kidx)),'x');
end
ap.wf_draw('ccf','k');
legend(["Excluded",string(1:n_k)])



%% Unit analysis

% Make index of recording per unit
[~,rec_idx] = unique([grp2idx(bhv.animal),bhv.days_from_learning],'rows');
unit_rec_idx = cell2mat(cellfun(@(rec,units) repmat(rec,length(units),1), ...
    num2cell(rec_idx),ephys.single_unit_idx,'uni',false));

% Concatenate unit data (responsive, single, psth, mean)
unit_resp_p_thresh = 0.05;
unit_resp_cat = cell2mat(cellfun(@(x) cell2mat(x)', ...
    horzcat(ephys.unit_resp_p_value{:})','uni',false)) < unit_resp_p_thresh;

unit_single_cat = cell2mat(cellfun(@logical,ephys.single_unit_idx,'uni',false));

unit_psth_cat = arrayfun(@(stim) cell2mat(cellfun(@(x) x{stim}, ...
    ephys.unit_smooth_event_psths,'ErrorHandler',@(varargin) [],'uni',false)), ...
    1:size(unit_resp_cat,2),'uni',false);

unit_mean_post_stim_cat = cell2mat(vertcat(ephys.mean_post_stim{:}));

% Make index of recording per unit
unit_animal = cell2mat(cellfun(@(grp,units) repmat(grp,length(units),1), ...
    num2cell(grp2idx(bhv.animal)),ephys.single_unit_idx,'uni',false));

% Make learning day per unit
unit_ld = cell2mat(cellfun(@(grp,units) repmat(grp,length(units),1), ...
    num2cell(bhv.days_from_learning),ephys.single_unit_idx,'uni',false));

% Make kmeans cluster per unit
unit_kidx_subset = cell2mat(cellfun(@(depth,kidx) kidx(depth(~isnan(depth))), ...
    ephys.unit_depth_group,kidx_rec,'uni',false));
unit_kidx = nan(size(unit_rec_idx));
unit_kidx(~isnan(cell2mat(ephys.unit_depth_group))) = unit_kidx_subset;


% Group data and plot
group_labels = [unit_rec_idx];
split_labels = [unit_ld,unit_kidx];

% (frac responsive cells)
curr_stim = 3;
use_units = true(size(unit_single_cat));
% use_units = unit_single_cat;
% use_units = logical(vertcat(ephys.str_fsi_idx{:}));

[unit_group_mean,groups] = ap.nestgroupfun({@mean,@mean}, ...
    +unit_resp_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

unit_group_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    +unit_resp_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

figure; 
h = tiledlayout(n_k,1);
plot_days = -3:2;
for curr_k = 1:n_k
    nexttile;
    plot_data = ismember(groups(:,1),plot_days) & groups(:,2) == curr_k;
    errorbar(groups(plot_data,1),unit_group_mean(plot_data,1), ...
        unit_group_sem(plot_data,1),'k','linewidth',2);
    ylabel('Frac responsive cells')
end
linkaxes(h.Children,'xy');

% (response amplitude)
curr_stim = 3;
% use_units = true(size(unit_kidx));
% use_units = unit_single_cat;
use_units = unit_resp_cat(:,3);

[unit_group_mean,groups] = ap.nestgroupfun({@mean,@mean}, ...
    unit_mean_post_stim_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

unit_group_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    unit_mean_post_stim_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

figure; hold on;
h = arrayfun(@(x) ap.errorfill(groups(groups(:,2) == x,1),unit_group_mean(groups(:,2) == x,1),unit_group_sem(groups(:,2) == x,1)),unique(groups(:,2)));
legend(h,string(num2cell(unique(groups(:,2)))));
ylabel('Amplitude responsive cells')


% (psth: average of single units)
curr_stim = 3;
% use_units = true(size(unit_kidx));
% use_units = unit_single_cat;
% use_units = unit_resp_cat(:,3);
use_units = logical(vertcat(ephys.str_msn_idx{:}));

[unit_group_mean,groups] = ap.nestgroupfun({@mean,@mean}, ...
    unit_psth_cat{curr_stim}(use_units,:),group_labels(use_units,:),split_labels(use_units,:));

unit_group_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    unit_psth_cat{curr_stim}(use_units,:),group_labels(use_units,:),split_labels(use_units,:));

plot_days = -3:2;

day_colormap_sym = ap.colormap('BKR',max(abs(plot_days))*2+1);
day_colormap = day_colormap_sym(ismember( ...
    -max(abs(plot_days)):max(abs(plot_days)),plot_days),:);

figure; 
h = tiledlayout(n_k,length(plot_days),'TileSpacing','none');
for curr_k = unique(groups(:,2))'
    for curr_day = plot_days
        nexttile; axis off;
        ap.errorfill([], ...
            unit_group_mean(groups(:,1) == curr_day & groups(:,2) == curr_k,:)', ...
            unit_group_sem(groups(:,1) == curr_day & groups(:,2) == curr_k,:)','k');
    end
end
linkaxes(h.Children,'xy');

figure; 
h = tiledlayout(n_k,1);
for curr_k = unique(groups(:,2))'
    nexttile; set(gca,'ColorOrder',day_colormap);
    arrayfun(@(x) ...
        ap.errorfill([], ...
        unit_group_mean(groups(:,1) == x & groups(:,2) == curr_k,:)', ...
        unit_group_sem(groups(:,1) == x & groups(:,2) == curr_k,:)'), ...
        plot_days);
end
linkaxes(h.Children,'xy');


%% Plot cell type heatmap

plot_celltype = vertcat(ephys.str_tan_idx{:});

plot_kidx = [1,2];
plot_days = -3:2;
plot_stim_idx = 3;

psth_t = -0.5:0.001:1;
stim_t = psth_t > 0 & psth_t < 0.2;
basieline_t = psth_t > -0.2 & psth_t < 0;

figure;
colormap(ap.colormap('BWR',[],2));
h = tiledlayout(2,length(plot_days),'TileIndexing','column');
title(h,'Learned day');
for curr_ld = plot_days
    nexttile;

    curr_data = unit_psth_cat{plot_stim_idx}(ismember(unit_kidx,plot_kidx) & unit_ld == curr_ld & ...
        unit_single_cat & plot_celltype,:);
    [~,sort_idx] = sort(max(curr_data(:,stim_t),[],2),'descend');
    imagesc(psth_t,[],curr_data(sort_idx,:));
    clim([-10,10])

    nexttile;
    curr_data = unit_psth_cat{plot_stim_idx}(ismember(unit_kidx,plot_kidx) & unit_ld == curr_ld & ...
        unit_single_cat & plot_celltype,:);
    ap.errorfill(psth_t,mean(curr_data,1),AP_sem(curr_data,1),'k');
end

linkaxes(h.Children(1:2:end),'xy');
linkaxes(h.Children(2:2:end),'xy'); 


unit_ld_prepost = unit_ld >= 0;
figure;
colormap(ap.colormap('BWR',[],2));
h = tiledlayout(2,2,'TileIndexing','column');
title(h,'Pre/post');
for curr_ld = unique(unit_ld_prepost)'
    nexttile;

    curr_data = unit_psth_cat{plot_stim_idx}(ismember(unit_kidx,plot_kidx) & unit_ld_prepost == curr_ld & ...
        unit_single_cat & plot_celltype,:);
    [~,sort_idx] = sort(max(curr_data(:,stim_t),[],2),'descend');
    imagesc(psth_t,[],curr_data(sort_idx,:));
    clim([-10,10])

    nexttile;
    curr_data = unit_psth_cat{plot_stim_idx}(ismember(unit_kidx,plot_kidx) & unit_ld_prepost == curr_ld & ...
        unit_single_cat & plot_celltype,:);
    ap.errorfill(psth_t,mean(curr_data,1),AP_sem(curr_data,1),'k');
end
linkaxes(h.Children(1:2:end),'xy');
linkaxes(h.Children(2:2:end),'xy'); 


%% MUA

animal_groupfun = @mean;

[psth_sum,psth_sum_kidx] = cellfun(@(mua,kidx) ap.groupfun(@sum,mua,[],[],kidx), ...
    ephys.binned_spikes_stim_align,kidx_rec,'uni',false);

% (currently psth time is hard coded: save this somewhere)
psth_t = -0.5:0.001:1;
baseline_t = psth_t < 0;
softnorm = 10;
psth_norm = cell(size(psth_sum));
psth_norm(~cellfun(@isempty,psth_sum)) = cellfun(@(mua) ...
    (mua - mean(mua(:,baseline_t,:),[1,2])) ./ ...
    (mean(mua(:,baseline_t,:),[1,2]) + softnorm), ...
    psth_sum(~cellfun(@isempty,psth_sum)),'uni',false);

psth_norm_cat_smooth = smoothdata(cell2mat(cellfun(@(x) reshape(permute(x,[1,3,2]),[], ...
    length(psth_t)),psth_norm,'uni',false)),2,'gaussian',[50,0]);

[stim_grp_sq,kidx_grp_sq] = cellfun(@(stim,kidx) ndgrid(stim,kidx), ...
    ephys.trial_stim_values,psth_sum_kidx,'uni',false);

animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(grp2idx(ephys.animal)),stim_grp_sq,'uni',false));

ld_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(bhv.days_from_learning),stim_grp_sq,'uni',false));

stim_grp = cell2mat(cellfun(@(x) reshape(x,[],1),stim_grp_sq,'uni',false));
kidx_grp = cell2mat(cellfun(@(x) reshape(x,[],1),kidx_grp_sq,'uni',false));

[psth_avg,psth_grp] = ap.nestgroupfun({@mean,animal_groupfun},psth_norm_cat_smooth,animal_grp, ...
    [ld_grp,stim_grp,kidx_grp]);
psth_sem = ap.nestgroupfun({@mean,@AP_sem},psth_norm_cat_smooth,animal_grp, ...
    [ld_grp,stim_grp,kidx_grp]);

max_t = psth_t > 0 & psth_t < 0.2;
psth_max = ap.nestgroupfun({@mean,animal_groupfun}, ...
    max(psth_norm_cat_smooth(:,max_t),[],2),animal_grp, ...
    [ld_grp,stim_grp,kidx_grp]);
psth_max_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    max(psth_norm_cat_smooth(:,max_t),[],2),animal_grp, ...
    [ld_grp,stim_grp,kidx_grp]);

plot_days = -3:2;
day_colormap_sym = ap.colormap('BKR',max(abs(plot_days))*2+1);
day_colormap = day_colormap_sym(ismember( ...
    -max(abs(plot_days)):max(abs(plot_days)),plot_days),:);

% Plot overlaid PSTH by learning day for stim
figure;
h = tiledlayout(n_k,length(unique(stim_grp)),'TileIndexing','column');
for curr_stim = unique(stim_grp)'
    for curr_k = unique(kidx_grp)'
        nexttile; hold on; set(gca,'ColorOrder',day_colormap);
        arrayfun(@(x) ...
            plot(psth_avg(psth_grp(:,1) == x & psth_grp(:,2) == curr_stim & psth_grp(:,3) == curr_k,:)'), ...
            plot_days);
    end
end
linkaxes(h.Children,'xy');

% Plot PSTH by learning day for stim
plot_stim = 90;
figure('Name','MUA by LD');
h = tiledlayout(n_k,length(plot_days),'TileSpacing','none');
for curr_k = unique(groups(:,2))'
    for curr_day = plot_days
        nexttile; axis off;
        plot_data = psth_grp(:,1) == curr_day & ...
            psth_grp(:,2) == plot_stim & ...
            psth_grp(:,3) == curr_k;
        if ~any(plot_data)
            continue
        end
        ap.errorfill(psth_t,psth_avg(plot_data,:),psth_sem(plot_data,:), ...
            day_colormap(curr_day==plot_days,:));
    end
end
linkaxes(h.Children,'xy');

% Plot PSTH max by learning day
figure;
h = tiledlayout(n_k,length(unique(stim_grp)),'TileIndexing','column');
for curr_stim = unique(stim_grp)'
    for curr_k = unique(kidx_grp)'
        nexttile; hold on; set(gca,'ColorOrder',day_colormap);
        plot_grps = ismember(psth_grp(:,1),plot_days) & ...
            psth_grp(:,2) == curr_stim & ...
            psth_grp(:,3) == curr_k;

        errorbar(plot_days,psth_max(plot_grps),psth_max_sem(plot_grps),'k','linewidth',2);
    end
end
linkaxes(h.Children,'xy');

% Plot pre/post-learn PSTH
ld_prepost = +(ld_grp >= 0) + 1;
[psth_learn_avg,psth_learn_grp] = ap.nestgroupfun({@mean,animal_groupfun},psth_norm_cat_smooth,animal_grp, ...
    [ld_prepost,stim_grp,kidx_grp]);
psth_learn_sem = ap.nestgroupfun({@mean,@AP_sem},psth_norm_cat_smooth,animal_grp, ...
    [ld_prepost,stim_grp,kidx_grp]);
prepost_cmap = [0,0,0.7;0.7,0,0.];

plot_stim = 90;
figure('Name','Post-learn MUA'); 
h = tiledlayout(n_k,2,'TileSpacing','none');
for curr_k = unique(groups(:,2))'
    for curr_ld = 1:2
        nexttile; axis off;
        plot_data = psth_learn_grp(:,1) == curr_ld & ...
            psth_learn_grp(:,2) == plot_stim & ...
            psth_learn_grp(:,3) == curr_k;
        if ~any(plot_data)
            continue
        end
        ap.errorfill(psth_t,psth_learn_avg(plot_data,:), ...
            psth_learn_sem(plot_data,:),prepost_cmap(curr_ld,:));
    end
end
linkaxes(h.Children,'xy');


%% Widefield


wf.V_stim_align








