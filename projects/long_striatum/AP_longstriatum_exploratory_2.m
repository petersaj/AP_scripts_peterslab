
%% Load ephys

data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');

load(fullfile(data_path,'swr_bhv'));
load(fullfile(data_path,'ctx_maps_to_str'));

load_workflow = 'task';
switch load_workflow
    case 'passive'
        load(fullfile(data_path,'ephys'));
    case 'task'
        load(fullfile(data_path,'task_ephys'));
        ephys = task_ephys;
        clear task_ephys;
end

U_master = plab.wf.load_master_U;

%% Load widefield

data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');

load(fullfile(data_path,'swr_bhv'));
load(fullfile(data_path,'ctx_wf'));
U_master = plab.wf.load_master_U;

%% Behavior

% Plot reaction time/index, by training day and learning day
training_day = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));

rxn_cat = cell2mat(bhv.stimwheel_rxn_mean);
rxn_null_cat = cell2mat(bhv.stimwheel_rxn_null_mean);
rxn_idx_cat = (rxn_null_cat-rxn_cat)./(rxn_null_cat+rxn_cat);

figure; h = tiledlayout(2,2);

for day_align = 1:2

    switch day_align
        case 1
            use_day_grp = training_day;
            plot_x = 1:7;
        case 2
            use_day_grp = bhv.days_from_learning;
            plot_x = -3:2;
    end

    [rxn_idx_mean,rxn_group_x] = ap.groupfun(@mean,rxn_idx_cat,use_day_grp);
    rxn_idx_sem = ap.groupfun(@AP_sem,rxn_idx_cat,use_day_grp);

    rxn_mean = ap.groupfun(@mean,rxn_cat,use_day_grp);
    rxn_sem = ap.groupfun(@AP_sem,rxn_cat,use_day_grp);

    rxn_null_mean = ap.groupfun(@mean,rxn_null_cat,use_day_grp);
    rxn_null_sem = ap.groupfun(@AP_sem,rxn_null_cat,use_day_grp);

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

% Plot histogram of learning days
n_learned_day = cell2mat(cellfun(@(x) find(bhv.learned_days(strcmp(bhv.animal,x)),1), ...
    unique(bhv.animal,'stable'),'uni',false));
figure;
histogram(n_learned_day,0.5:max(n_learned_day)+0.5)
ylabel('Number of mice');
xlabel('Days to learn');


%% K-means on maps 

% K-means cluster maps
maps_cat = cat(3,all_ctx_maps_to_str.cortex_kernel_px{:});
expl_var_cat = vertcat(all_ctx_maps_to_str.explained_var{:});

n_k = 4;
kmeans_starting = nanmean(cell2mat(permute(cellfun(@(x) ...
    x(:,:,round(linspace(1,size(x,3),n_k))), ...
    all_ctx_maps_to_str.cortex_kernel_px(~cellfun(@isempty, ...
    all_ctx_maps_to_str.cortex_kernel_px)),'uni',false),[2,3,4,1])),4);

[kidx,kmeans_centroid_flat] = kmeans(...
    reshape(maps_cat,prod(size(U_master,[1,2])),[])',n_k, ...
    'Distance','Correlation','start',reshape(kmeans_starting,[],n_k)');
kmeans_centroid = reshape(kmeans_centroid_flat',size(U_master,1),size(U_master,2),[]);

% %%%%%% (enforce increasing depth order)
% kidx_rec = mat2cell(kidx,cellfun(@(x) size(x,3).*(size(x,1)>0), ...
%     all_ctx_maps_to_str.cortex_kernel_px));
% kidx = cell2mat(cellfun(@cummax,kidx_rec,'uni',false));
% %%%%%%

kmeans_cluster_mean = ap.groupfun(@nanmean,maps_cat,[],[],kidx);

kidx_rec = mat2cell(kidx,cellfun(@(x) size(x,3).*(size(x,1)>0), ...
    all_ctx_maps_to_str.cortex_kernel_px));

figure;
h = tiledlayout(n_k,2,'tilespacing','none');
for curr_k = 1:n_k
    nexttile;
    imagesc(kmeans_centroid(:,:,curr_k));
    axis image off
    clim(max(abs(clim)).*[-1,1]);
    colormap(ap.colormap('PWG'));
    ap.wf_draw('ccf','k');

    nexttile;
    imagesc(kmeans_cluster_mean(:,:,curr_k));
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

plot_day_bins = [-Inf,-2:1,Inf];
unit_plot_days_grp = discretize(max(unit_ld,-inf),plot_day_bins);

% Make kmeans cluster per unit
unit_kidx_subset = cell2mat(cellfun(@(depth,kidx) kidx(depth(~isnan(depth))), ...
    ephys.unit_depth_group,kidx_rec,'uni',false));
unit_kidx = nan(size(unit_rec_idx));
unit_kidx(~isnan(cell2mat(ephys.unit_depth_group))) = unit_kidx_subset;


% Group data and plot
group_labels = [unit_rec_idx];
split_labels = [unit_plot_days_grp,unit_kidx];

% (frac responsive cells)
curr_stim = 3;
% use_units = true(size(unit_single_cat));
% use_units = unit_single_cat;
use_units = logical(vertcat(ephys.str_fsi_idx{:}));

[unit_group_mean,groups] = ap.nestgroupfun({@mean,@mean}, ...
    +unit_resp_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

unit_group_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    +unit_resp_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

[n_units_ld,n_units_ld_group] = ap.nestgroupfun({@sum,@mean},+use_units,group_labels,split_labels);

figure; 
h = tiledlayout(n_k,2);
for curr_k = 1:n_k
    nexttile;
    plot_data = groups(:,2) == curr_k;
    errorbar(groups(plot_data,1),unit_group_mean(plot_data,1), ...
        unit_group_sem(plot_data,1),'k','linewidth',2);
    ylabel('Frac responsive cells');

    nexttile;
    plot_data = n_units_ld_group(:,2) == curr_k;
    plot(n_units_ld_group(plot_data,1),n_units_ld(plot_data),'k');
    ylabel('N units');
end
linkaxes(h.Children(1:2:end),'xy');
linkaxes(h.Children(2:2:end),'xy');


% (psth: average of single units)
curr_stim = 2;
use_units = true(size(unit_kidx));
% use_units = unit_single_cat;
% use_units = unit_resp_cat(:,3);
% use_units = logical(vertcat(ephys.str_msn_idx{:}));

[unit_group_mean,groups] = ap.nestgroupfun({@mean,@mean}, ...
    unit_psth_cat{curr_stim}(use_units,:),group_labels(use_units,:),split_labels(use_units,:));

unit_group_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    unit_psth_cat{curr_stim}(use_units,:),group_labels(use_units,:),split_labels(use_units,:));

day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

figure; 
h = tiledlayout(n_k,max(plot_days_grp),'TileSpacing','none');
for curr_k = unique(groups(:,2))'
    for curr_day = unique(plot_days_grp)'
        nexttile; axis off;
        ap.errorfill([], ...
            unit_group_mean(groups(:,1) == curr_day & groups(:,2) == curr_k,:)', ...
            unit_group_sem(groups(:,1) == curr_day & groups(:,2) == curr_k,:)','k');
    end
end
linkaxes(h.Children,'xy');




%% Plot cell type heatmap

plot_celltype = logical(vertcat(ephys.str_tan_idx{:}));
% plot_celltype = vertcat(ephys.str_fsi_idx{:}) & vertcat(ephys.single_unit_idx{:});

plot_kidx = [1:2];
plot_stim_idx = 2;

psth_t = -0.5:0.001:1;
stim_t = psth_t > 0 & psth_t < 0.2;


% Plot grouped days

% unit_ld_prepost = unit_ld >= 0;
plot_day_bins = [-Inf,-2:1,Inf];
unit_plot_day_grp = discretize(unit_ld,plot_day_bins);

figure;
colormap(ap.colormap('BWR',[],2));
h = tiledlayout(2,max(unit_plot_day_grp),'TileIndexing','column','TileSpacing','compact');
for curr_ld = unique(unit_plot_day_grp(~isnan(unit_plot_day_grp)))'
    nexttile;

    curr_data = unit_psth_cat{plot_stim_idx}(ismember(unit_kidx,plot_kidx) & unit_plot_day_grp == curr_ld & ...
        plot_celltype,:);
    [~,sort_idx] = sort(max(curr_data(:,stim_t),[],2),'descend');
    imagesc(psth_t,[],curr_data(sort_idx,:));
    clim([-10,10])

    nexttile;
    curr_data = unit_psth_cat{plot_stim_idx}(ismember(unit_kidx,plot_kidx) & unit_plot_day_grp == curr_ld & ...
        plot_celltype,:);
    ap.errorfill(psth_t,mean(curr_data,1),AP_sem(curr_data,1),'k');
end
linkaxes(h.Children(1:2:end),'xy');
linkaxes(h.Children(2:2:end),'xy'); 



%% MUA (passive)

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

plot_day_bins = [-Inf,-2:1,Inf];
plot_days_grp = discretize(max(ld_grp,-inf),plot_day_bins);

stim_grp = cell2mat(cellfun(@(x) reshape(x,[],1),stim_grp_sq,'uni',false));
kidx_grp = cell2mat(cellfun(@(x) reshape(x,[],1),kidx_grp_sq,'uni',false));

[psth_avg,psth_grp] = ap.nestgroupfun({@mean,animal_groupfun},psth_norm_cat_smooth,animal_grp, ...
    [plot_days_grp,stim_grp,kidx_grp]);
psth_sem = ap.nestgroupfun({@mean,@AP_sem},psth_norm_cat_smooth,animal_grp, ...
    [plot_days_grp,stim_grp,kidx_grp]);

max_t = psth_t > 0 & psth_t < 0.2;
psth_max = ap.nestgroupfun({@mean,animal_groupfun}, ...
    max(psth_norm_cat_smooth(:,max_t),[],2),animal_grp, ...
    [plot_days_grp,stim_grp,kidx_grp]);
psth_max_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    max(psth_norm_cat_smooth(:,max_t),[],2),animal_grp, ...
    [plot_days_grp,stim_grp,kidx_grp]);


day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

% Plot overlaid PSTH by learning day for stim
figure;
h = tiledlayout(n_k,length(unique(stim_grp)),'TileIndexing','column');
for curr_stim = unique(stim_grp)'
    for curr_k = unique(kidx_grp)'
        nexttile; hold on; set(gca,'ColorOrder',day_colormap);
        plot(psth_avg(psth_grp(:,2) == curr_stim & psth_grp(:,3) == curr_k,:)')
    end
end
linkaxes(h.Children,'xy');

% Plot PSTH by learning day for stim
plot_stim = 90;
figure('Name','MUA by LD');
h = tiledlayout(n_k,max(psth_grp(:,1)),'TileSpacing','none');
for curr_k = unique(psth_grp(:,3))'
    for curr_day = unique(psth_grp(:,1))'
        nexttile; axis off;
        plot_data = psth_grp(:,1) == curr_day & ...
            psth_grp(:,2) == plot_stim & ...
            psth_grp(:,3) == curr_k;
        if ~any(plot_data)
            continue
        end
        ap.errorfill(psth_t,psth_avg(plot_data,:),psth_sem(plot_data,:), ...
            day_colormap(curr_day,:));
    end
end
linkaxes(h.Children,'xy');

% Plot PSTH max by learning day
figure;
h = tiledlayout(n_k,length(unique(stim_grp)),'TileIndexing','column');
for curr_stim = unique(stim_grp)'
    for curr_k = unique(kidx_grp)'
        nexttile; hold on; set(gca,'ColorOrder',day_colormap);
        plot_grps = psth_grp(:,2) == curr_stim & ...
            psth_grp(:,3) == curr_k;

        errorbar(psth_max(plot_grps),psth_max_sem(plot_grps),'k','linewidth',2);
        xline(find(plot_day_bins == 0));
    end
end
linkaxes(h.Children,'xy');

% Plot heatmaps sorted by reaction times (animals separately)
animals = unique(ephys.animal,'stable');
plot_k = 1;
figure; tiledlayout(max(plot_days_grp),max(animal_grp), ...
    'TileIndexing','columnmajor','TileSpacing','none');
for curr_animal = 1:max(animal_grp)
    for curr_day = 1:max(plot_days_grp)
        nexttile;
        curr_trials = find(animal_grp == curr_animal & stim_grp == 90 & ...
            plot_days_grp == curr_day & kidx_grp == plot_k);
        imagesc(psth_t,[],movmean(psth_norm_cat_smooth(curr_trials,:),[3,1]));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;
        if curr_day == plot_days(1)
            title(animals{curr_animal});
        end
    end
end



%% MUA (task)

[psth_sum,psth_sum_kidx] = cellfun(@(mua,kidx) ap.groupfun(@sum,mua,[],[],kidx), ...
    ephys.binned_spikes_stim_align,kidx_rec,'uni',false);

% (currently psth time is hard coded: save this somewhere)
psth_t = -0.5:0.001:1;
baseline_t = psth_t < -0.3;
softnorm = 10;
psth_norm = cell(size(psth_sum));
psth_norm(~cellfun(@isempty,psth_sum)) = cellfun(@(mua) ...
    (mua - mean(mua(:,baseline_t,:),[1,2])) ./ ...
    (mean(mua(:,baseline_t,:),[1,2]) + softnorm), ...
    psth_sum(~cellfun(@isempty,psth_sum)),'uni',false);

psth_norm_cat_smooth = smoothdata(cell2mat(cellfun(@(x) reshape(permute(x,[1,3,2]),[], ...
    length(psth_t)),psth_norm,'uni',false)),2,'gaussian',[30,0]);

[rxn_grp_sq,kidx_grp_sq] = cellfun(@(stim,kidx) ndgrid(stim,kidx), ...
    bhv.stim_to_move,psth_sum_kidx,'uni',false);

animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(grp2idx(ephys.animal)),rxn_grp_sq,'uni',false));

training_day = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));
td_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(training_day),rxn_grp_sq,'uni',false));

ld_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(bhv.days_from_learning),rxn_grp_sq,'uni',false));

rxn_grp = cell2mat(cellfun(@(x) reshape(x,[],1),rxn_grp_sq,'uni',false));
kidx_grp = cell2mat(cellfun(@(x) reshape(x,[],1),kidx_grp_sq,'uni',false));


day_grp = discretize(max(ld_grp,-inf),[-Inf,-2:1,Inf]);

% Plot heatmaps sorted by reaction times
figure; tiledlayout(n_k,max(day_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(day_grp)'
        nexttile;
        curr_trials = find(day_grp == curr_day & kidx_grp == curr_k);
       
        [~,sort_idx] = sort(rxn_grp(curr_trials));
        imagesc(psth_t,[],movmean(psth_norm_cat_smooth(curr_trials(sort_idx),:),[20,5]));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;
    end
end

% Plot heatmaps sorted by reaction times (animals separately)
% (currently plotted as training day)
n_learned_day = cell2mat(cellfun(@(x) min([NaN,find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'),'uni',false));

plot_days = 1:8;
animals = unique(ephys.animal,'stable');
plot_k = 1;
figure; tiledlayout(length(plot_days),max(animal_grp), ...
    'TileIndexing','columnmajor','TileSpacing','none');
for curr_animal = 1:max(animal_grp)
    for curr_day = plot_days
        nexttile;
        curr_trials = find(animal_grp == curr_animal & td_grp == curr_day & kidx_grp == plot_k);
        [~,sort_idx] = sort(rxn_grp(curr_trials));
        imagesc(psth_t,[],movmean(psth_norm_cat_smooth(curr_trials(sort_idx),:),[3,1]));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;
        if curr_day == plot_days(1)
            title(animals{curr_animal});
        end
        if curr_day >= n_learned_day(curr_animal)
            xline(0,'g');
        end
    end
end

% Choose split type
split_type = 'rxn_percentile';
switch split_type
    case 'percentile'
        n_split = 5;
    case 'rxn_percentile'
        n_split = 5;
    case 'rxn_bin'
%           rxn_bins = [0,0.2,0.4,0.6,0.8,1];
%         rxn_bins = [-Inf,0.05,0.15,0.25,0.5,1,Inf];
        rxn_bins = [0,0.5,2,Inf];
        n_split = length(rxn_bins)-1;
end

% Plot PSTH split
figure('color','w'); h = tiledlayout(n_k,max(day_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(day_grp)'
        nexttile;
        curr_trials = find(day_grp == curr_day & kidx_grp == curr_k);
        curr_data = psth_norm_cat_smooth(curr_trials,:);

        switch split_type
            case 'percentile'
                split_idx = cell2mat(arrayfun(@(x) ...
                    min(round(linspace(1,n_split+1,sum(animal_grp(curr_trials)==x))),n_split)', ...                
                unique(animal_grp(curr_trials)),'uni',false));
            case 'rxn_percentile'
%                 % (total)
%                 split_idx_unsorted = min(round(linspace(1,n_split+1,size(curr_data,1))),n_split)';
%                 [~,sort_idx] = sort(rxn_grp(curr_trials));
%                 [~,rxn_rank] = sort(sort_idx);
%                 split_idx = split_idx_unsorted(rxn_rank);

                % (by animal)
                split_idx = cell2mat(arrayfun(@(x) ...
                    discretize(tiedrank(rxn_grp(intersect(curr_trials,find(animal_grp == x)))), ...
                    linspace(1,length(intersect(curr_trials,find(animal_grp == x))),n_split+1)), ...
                    unique(animal_grp(curr_trials)),'uni',false));
            case 'rxn_bin'
                split_idx = discretize(rxn_grp(curr_trials),rxn_bins);
        end

        plot(psth_t,ap.nestgroupfun({@mean,@mean},curr_data,animal_grp(curr_trials),split_idx)');

    end
end
linkaxes(h.Children,'xy');
set(h.Children,'colororder',ap.colormap('KR',n_split)); 
axis(h.Children,'off')

% Plot PSTH max split
figure('color','w'); h = tiledlayout(n_k,max(day_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day = unique(day_grp)'
        nexttile;
        curr_trials = find(day_grp == curr_day & kidx_grp == curr_k);
        curr_data = psth_norm_cat_smooth(curr_trials,:);

        switch split_type
            case 'percentile'
                split_idx = cell2mat(arrayfun(@(x) ...
                    min(round(linspace(1,n_split+1,sum(animal_grp(curr_trials)==x))),n_split)', ...
                    unique(animal_grp(curr_trials)),'uni',false));
            case 'rxn_percentile'
%                 % (total)
%                 split_idx_unsorted = min(round(linspace(1,n_split+1,size(curr_data,1))),n_split)';
%                 [~,sort_idx] = sort(rxn_grp(curr_trials));
%                 [~,rxn_rank] = sort(sort_idx);
%                 split_idx = split_idx_unsorted(rxn_rank);

                % (by animal)
                split_idx = cell2mat(arrayfun(@(x) ...
                    discretize(tiedrank(rxn_grp(intersect(curr_trials,find(animal_grp == x)))), ...
                    linspace(1,length(intersect(curr_trials,find(animal_grp == x))),n_split+1)), ...
                    unique(animal_grp(curr_trials)),'uni',false));
            case 'rxn_bin'
                split_idx = discretize(rxn_grp(curr_trials),rxn_bins);
        end

        data_t = psth_t > 0.05 & psth_t < 0.15;

        % (mean PSTH, max t, mean across animals)
        [psth_mean,psth_mean_grp] = ap.nestgroupfun({@mean,@mean},curr_data, ...
            (1:size(curr_data,1))',[animal_grp(curr_trials),split_idx]);
        psth_max = max(psth_mean(:,data_t),[],2);
        psth_max_avg = ap.nestgroupfun({@mean,@mean},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
        psth_max_sem = ap.nestgroupfun({@mean,@AP_sem},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
        errorbar(psth_max_avg,psth_max_sem,'k','linewidth',2);

%         % (mean in time window)
%         curr_data_point = mean(curr_data(:,data_t),2);
%         data_avg = ap.nestgroupfun({@mean,@mean},curr_data_point,animal_grp(curr_trials),split_idx);
%         data_sem = ap.nestgroupfun({@mean,@AP_sem},curr_data_point,animal_grp(curr_trials),split_idx);
%         errorbar(data_avg,data_sem,'k','linewidth',2);

%         % (testing correcting for movement)
%         [psth_mean,psth_mean_grp] = ap.nestgroupfun({@mean,@mean},curr_data, ...
%             (1:size(curr_data,1))',[animal_grp(curr_trials),split_idx]);
%         psth_max = max(psth_mean(:,data_t),[],2);
% 
%         r_t = psth_t >= 0.4 & psth_t < 0.6;
%         r = max(psth_mean(:,r_t),[],2);
%         r_avg = ap.nestgroupfun({@mean,@mean},r,psth_mean_grp(:,1),psth_mean_grp(:,2));
% 
%         psth_max_avg = ap.nestgroupfun({@mean,@mean},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
%         psth_max_sem = ap.nestgroupfun({@mean,@AP_sem},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
%         errorbar(psth_max_avg - r_avg(3),psth_max_sem,'k','linewidth',2);

    end
end
linkaxes(h.Children,'xy');
set(h.Children,'colororder',ap.colormap('KR',n_split)); 


%% (testing from above - animal corr for fast rxn and next day)

% Choose split type
split_type = 'rxn_percentile';
switch split_type
    case 'percentile'
        n_split = 5;
    case 'rxn_percentile'
        n_split = 5;
    case 'rxn_bin'
%           rxn_bins = [0,0.2,0.4,0.6,0.8,1];
%         rxn_bins = [-Inf,0.05,0.15,0.25,0.5,1,Inf];
        rxn_bins = [0,0.5,2,Inf];
        n_split = length(rxn_bins)-1;
end


% Plot PSTH max split
td_grp
animal_grp
rxn_grp

split_idx = cell2mat(arrayfun(@(x) ...
    min(round(linspace(1,n_split+1,sum(animal_grp(curr_trials)==x))),n_split)', ...
    unique(animal_grp(curr_trials)),'uni',false));

% (mean PSTH, max t, mean across animals)
data_t = psth_t > 0.05 & psth_t < 0.15;

[psth_mean,psth_mean_grp] = ap.nestgroupfun({@mean,@mean},curr_data, ...
    (1:size(curr_data,1))',[animal_grp(curr_trials),split_idx]);
psth_max = max(psth_mean(:,data_t),[],2);
psth_max_avg = ap.nestgroupfun({@mean,@mean},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
psth_max_sem = ap.nestgroupfun({@mean,@AP_sem},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
errorbar(psth_max_avg,psth_max_sem,'k','linewidth',2);



figure('color','w'); h = tiledlayout(n_k,max(day_grp),'TileSpacing','compact');
curr_k = 1;
for curr_day = unique(day_grp)'
    nexttile;
    curr_trials = find(day_grp == curr_day & kidx_grp == curr_k);
    curr_data = psth_norm_cat_smooth(curr_trials,:);

    switch split_type
        case 'percentile'
            split_idx = cell2mat(arrayfun(@(x) ...
                min(round(linspace(1,n_split+1,sum(animal_grp(curr_trials)==x))),n_split)', ...
                unique(animal_grp(curr_trials)),'uni',false));
        case 'rxn_percentile'
            %                 % (total)
            %                 split_idx_unsorted = min(round(linspace(1,n_split+1,size(curr_data,1))),n_split)';
            %                 [~,sort_idx] = sort(rxn_grp(curr_trials));
            %                 [~,rxn_rank] = sort(sort_idx);
            %                 split_idx = split_idx_unsorted(rxn_rank);

            % (by animal)
            split_idx = cell2mat(arrayfun(@(x) ...
                discretize(tiedrank(rxn_grp(intersect(curr_trials,find(animal_grp == x)))), ...
                linspace(1,length(intersect(curr_trials,find(animal_grp == x))),n_split+1)), ...
                unique(animal_grp(curr_trials)),'uni',false));
        case 'rxn_bin'
            split_idx = discretize(rxn_grp(curr_trials),rxn_bins);
    end

    
    % (mean PSTH, max t, mean across animals)
    data_t = psth_t > 0.05 & psth_t < 0.15;

    [psth_mean,psth_mean_grp] = ap.nestgroupfun({@mean,@mean},curr_data, ...
        (1:size(curr_data,1))',[animal_grp(curr_trials),split_idx]);
    psth_max = max(psth_mean(:,data_t),[],2);
    psth_max_avg = ap.nestgroupfun({@mean,@mean},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
    psth_max_sem = ap.nestgroupfun({@mean,@AP_sem},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
    errorbar(psth_max_avg,psth_max_sem,'k','linewidth',2);

end

linkaxes(h.Children,'xy');
set(h.Children,'colororder',ap.colormap('KR',n_split));



%% Widefield

% Create group indicies
wf_t = wf.wf_stim_time{1};

wf_animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,length(grp),1),    ...
    num2cell(grp2idx(wf.animal)),wf.trial_stim_values,'uni',false));

wf_ld_grp = cell2mat(cellfun(@(ld,grp) repmat(ld,length(grp),1),    ...
    num2cell(bhv.days_from_learning),wf.trial_stim_values,'uni',false));

plot_day_bins = [-Inf,-2:1,Inf];
wf_plot_days_grp = discretize(max(wf_ld_grp,-inf),plot_day_bins);


% Average and plot widefield
[wf_avg,wf_avg_grp] = ap.nestgroupfun({@mean,@mean},cell2mat(wf.V_no_move_stim_align), ...
    wf_animal_grp,[wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

curr_stim = 90;

curr_data_idx = wf_avg_grp(:,2) == curr_stim;

curr_data_px = plab.wf.svd2px(U_master,permute(wf_avg(curr_data_idx,:,:),[3,2,1]));

ap.imscroll(curr_data_px - nanmean(curr_data_px(:,:,wf_t<0,:),3),wf_t)
colormap(ap.colormap('PWG',[],1.5));
clim(max(curr_data_px(:)).*[-1,1]*0.6);
axis image;
ap.wf_draw('ccf','k');

% Plot ROIs by striatum cluster
% % (weighted average)
% striatum_wf_roi = max(kmeans_cluster_mean,0)./max(kmeans_cluster_mean,[],[1,2]);
% (thresholded spot)
[~,m] = max(kmeans_cluster_mean,[],[1,2],'linear');
[mr,mc] = ind2sub(size(U_master,[1,2]),m);

striatum_wf_roi = kmeans_cluster_mean > std(kmeans_cluster_mean(:))*3;
striatum_wf_roi(:,round(size(striatum_wf_roi,2)/2):end,:) = false;
for k = 1:size(striatum_wf_roi,3)
    [~,m] = max(imgaussfilt(kmeans_cluster_mean(:,:,k),10).* ...
        round(linspace(1,0,size(kmeans_cluster_mean,2))),[],[1,2],'linear');
    [mx,my] = ind2sub(size(U_master,[1,2]),m);
    striatum_wf_roi(:,:,k) = bwselect(striatum_wf_roi(:,:,k),my,mx);
end

figure;tiledlayout(n_k,1,'tilespacing','none')
for curr_k = 1:n_k
    nexttile;
    imagesc(striatum_wf_roi(:,:,curr_k));
    axis image off; ap.wf_draw('ccf','r');
    colormap(ap.colormap('WK'));
end

wf_striatum_roi = permute(ap.wf_roi(U_master, ...
    permute(cell2mat(wf.V_no_move_stim_align),[3,2,1]),[],[],striatum_wf_roi),[3,2,1]);

[wf_striatum_roi_avg,wf_striatum_roi_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi, ...
    wf_animal_grp,[wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

wf_striatum_roi_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi, ...
    wf_animal_grp,[wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

plot_stim = 90;
figure;
h = tiledlayout(n_k,max(wf_plot_days_grp),'TileSpacing','none');
for curr_k = 1:size(striatum_wf_roi,3)
    for curr_ld = unique(wf_plot_days_grp)'
        nexttile; axis off;
        curr_data_idx = wf_striatum_roi_avg_grp(:,1) == curr_ld & ...
            wf_striatum_roi_avg_grp(:,2) == plot_stim;

        ap.errorfill(wf_t,wf_striatum_roi_avg(curr_data_idx,:,curr_k), ...
            wf_striatum_roi_sem(curr_data_idx,:,curr_k), ...
            day_colormap(curr_ld,:));
    end
end
linkaxes(h.Children,'xy');

% Plot average in ROI by day
stim_t = wf_t >= 0 & wf_t <= 0.25;
wf_striatum_roi_tavg = nanmean(wf_striatum_roi(:,stim_t,:),2);

[wf_striatum_roi_tavg_avg,wf_striatum_roi_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi_tavg, ...
    wf_animal_grp,[wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

wf_striatum_roi_tavg_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi_tavg, ...
    wf_animal_grp,[wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

plot_stim = 90;
figure;
h = tiledlayout(n_k,1,'TileSpacing','none');
for curr_k = 1:size(striatum_wf_roi,3)
        nexttile;
        curr_data_idx = wf_striatum_roi_avg_grp(:,2) == plot_stim;

        errorbar(wf_striatum_roi_tavg_avg(curr_data_idx,:,curr_k), ...
            wf_striatum_roi_tavg_sem(curr_data_idx,:,curr_k),'k','linewidth',2);
        ylabel('t avg')
        axis padded
end
linkaxes(h.Children,'xy');

% Plot max in ROI by day
stim_t = wf_t >= 0 & wf_t <= 0.25;

[wf_striatum_roi_dayavg,wf_striatum_roi_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi, ...
    (1:size(wf_striatum_roi,1))', ...
    [wf_animal_grp,wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

wf_striatum_roi_max = max(wf_striatum_roi_dayavg(:,stim_t,:),[],2);

[wf_striatum_roi_max_avg,wf_striatum_roi_max_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi_max, ...
    wf_striatum_roi_grp(:,1),wf_striatum_roi_grp(:,2:3));

wf_striatum_roi_max_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi_max, ...
    wf_striatum_roi_grp(:,1),wf_striatum_roi_grp(:,2:3));

plot_stim = 90;
figure;
h = tiledlayout(n_k,1,'TileSpacing','none');
for curr_k = 1:size(striatum_wf_roi,3)
        nexttile;
        curr_data_idx = wf_striatum_roi_max_avg_grp(:,2) == plot_stim;

        errorbar(wf_striatum_roi_max_avg(curr_data_idx,:,curr_k), ...
            wf_striatum_roi_max_sem(curr_data_idx,:,curr_k),'k','linewidth',2);
        ylabel('t max');
        axis padded
end

linkaxes(h.Children,'xy');


















