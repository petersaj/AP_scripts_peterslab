%% Notes

% After AM data packaging


%% Load ephys

% AM-packaged: 
data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');

load_workflow = 'passive';
switch load_workflow
    case 'passive'
        load(fullfile(data_path,'ephys'));
    case 'task'
        load(fullfile(data_path,'task_ephys'));
        ephys = task_ephys;
        clear task_ephys;
    case 'task_outcome'
        load(fullfile(data_path,'task_ephys_outcome'));
        ephys = task_ephys;
        clear task_ephys;
end


%% Load ephys (no-stim)

data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data\AM_packaged\nostim\';

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


%% Load widefield

data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');
U_master = plab.wf.load_master_U;

load_workflow = 'passive';
switch load_workflow
    case 'passive'
        load(fullfile(data_path,'ctx_wf'));
    case 'task'
        load(fullfile(data_path,'task_ctx_wf'));
end

%% Behavior

% Load behavior
data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');
load(fullfile(data_path,'swr_bhv'));

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

%% Behavior v2
% (this one includes other stats - choose how to define learning)

use_stat = 'firstmove_mean';
learn_stat_p = 0.05;

% Load behavior
data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');
load(fullfile(data_path,'swr_bhv_v2'));

% Set (overwrite) "learned_days" and "days_from_learning"
bhv.learned_days = cellfun(@(x) x < learn_stat_p,bhv.(['stimwheel_pval_',use_stat]));
for curr_animal = unique(bhv.animal)'
    curr_rows = strcmp(bhv.animal,curr_animal);
    bhv.days_from_learning(curr_rows) = ...
        (1:sum(curr_rows))' - ...
        max([NaN,find(bhv.learned_days(curr_rows),1)]);
end

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

% Load maps
data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');
load(fullfile(data_path,'ctx_maps_to_str'));
U_size = size(all_ctx_maps_to_str.cortex_kernel_px{1},[1,2]);

% %%%%%%%%%%%%%%%%%%%
% %%% Blur maps
% all_ctx_maps_to_str.cortex_kernel_px = ...
%     cellfun(@(x) imgaussfilt(x,5),all_ctx_maps_to_str.cortex_kernel_px,'uni',false);
% %%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%
% %%% Rectify maps
% all_ctx_maps_to_str.cortex_kernel_px = ...
%     cellfun(@(x) max(x,0,'IncludeMissing'),all_ctx_maps_to_str.cortex_kernel_px,'uni',false);
% %%%%%%%%%%%%%%%%%%%

% K-means cluster maps
maps_flat_cat = reshape(cat(3,all_ctx_maps_to_str.cortex_kernel_px{:}),prod(U_size),[]);
expl_var_cat = vertcat(all_ctx_maps_to_str.explained_var{:});

n_k = 6;

% % (kmeans starting as index)
% kmeans_starting = nanmean(cell2mat(permute(cellfun(@(x) ...
%     x(:,:,round(linspace(1,size(x,3),n_k))), ...
%     all_ctx_maps_to_str.cortex_kernel_px(~cellfun(@isempty, ...
%     all_ctx_maps_to_str.cortex_kernel_px)),'uni',false),[2,3,4,1])),4);

% % (kmeans depth avg using discretize)
% kmeans_starting = nanmean(cell2mat(permute(cellfun(@(x) ap.groupfun(@mean,x,[],[],discretize(1:size(x,3),n_k)), ...
%     all_ctx_maps_to_str.cortex_kernel_px(cellfun(@(x) ...
%     size(x,3) >= n_k,all_ctx_maps_to_str.cortex_kernel_px)),'uni',false),[2,3,4,1])),4);

% (kmeans depth avg using quantiles)
kmeans_starting = nanmean(cell2mat(permute(cellfun(@(x) ...
    ap.groupfun(@mean,x,[],[],ap.quantile_bin(size(x,3),n_k)), ...
    all_ctx_maps_to_str.cortex_kernel_px(cellfun(@(x) ...
    size(x,3) >= n_k,all_ctx_maps_to_str.cortex_kernel_px)),'uni',false),[2,3,4,1])),4);

% (testing: depth by day)
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


% %%%%% BINARY KMEANS
% [~,m] = max(kmeans_starting,[],[1,2],'linear');
% [mr,mc] = ind2sub(size(U_master,[1,2]),m);
% kmeans_starting_binary = kmeans_starting > max(kmeans_starting,[],[1,2])*0.5;
% kmeans_starting_binary(:,round(size(kmeans_starting_binary,2)/2):end,:) = false;
% for k = 1:size(kmeans_starting_binary,3)
%     [~,m] = max(imgaussfilt(kmeans_cluster_mean(:,:,k),10).* ...
%         round(linspace(1,0,size(kmeans_cluster_mean,2))),[],[1,2],'linear');
%     [mx,my] = ind2sub(size(U_master,[1,2]),m);
%     kmeans_starting_binary(:,:,k) = bwselect(kmeans_starting_binary(:,:,k),my,mx);
% end
% 
% [domain_idx,kmeans_centroid_flat,~,kmeans_dist] = kmeans(...
%     (maps_flat_cat> prctile(maps_flat_cat,50,1))',n_k, ...
%     'Distance','hamming','start',reshape(+kmeans_starting_binary,[],n_k)');
% kmeans_centroid = reshape(kmeans_centroid_flat',size(kmeans_starting,1),size(kmeans_starting,2),[]);
% %%%%%%

[domain_idx,kmeans_centroid_flat,~,kmeans_dist] = kmeans(...
    maps_flat_cat',n_k, ...
    'Distance','correlation','start',reshape(kmeans_starting,[],n_k)');
kmeans_centroid = reshape(kmeans_centroid_flat',size(kmeans_starting,1),size(kmeans_starting,2),[]);

% %%%%%%%%%%%%%%%%%%%
% % ADJUSTMENTS FOR N_K=4
% % looks like value of ~0.75 works, just picked empirical way to get there
% domain_idx_raw = domain_idx;
% 
% dist_cutoff_k1 = prctile(kmeans_dist(domain_idx_raw==1,1),95);
% dist_cutoff_k2 = prctile(kmeans_dist(domain_idx_raw==2,2),95);
% 
% domain_idx = domain_idx_raw;
% 
% % (to fold 2's into cluster 1)
% domain_idx(domain_idx_raw==2 & kmeans_dist(:,1) < dist_cutoff_k2) = 1;
% 
% % % (make 2's into new cluster)
% % domain_idx(domain_idx_raw==2 & kmeans_dist(:,1) < dist_cutoff) = 5;
% % [~,domain_idx] = ismember(domain_idx,[1,5,2,3,4]); domain_idx(domain_idx==0) = NaN;% re-assign to order depths
% % n_k = 5;
% % % domain_idx(domain_idx==3) = 2; domain_idx(domain_idx==4) = 3;n_k = 3; % (combine pfc 2+3)
% 
% % % (make 1's and 2's into new cluster)
% % domain_idx(domain_idx_raw==1 & kmeans_dist(:,2) < dist_cutoff_k1) = 5;
% % domain_idx(domain_idx_raw==2 & kmeans_dist(:,1) < dist_cutoff_k2) = 5;
% % [~,domain_idx] = ismember(domain_idx,[1,5,2,3,4]); domain_idx(domain_idx==0) = NaN;% re-assign to order depths
% % n_k = 5;
% 
% %%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%
if n_k == 6
    % % ADJUSTMENTS FOR N_K=6
    % Turn 6 into 3 (vis, frontal, lateral)
    domain_idx(ismember(domain_idx,2:4)) = 2;
    domain_idx(ismember(domain_idx,5:6)) = 3;
    n_k = 3;
end
% %%%%%%%%%%%%%%%%%%%

domain_idx_rec = mat2cell(domain_idx,cellfun(@(x) size(x,3).*(size(x,1)>0), ...
    all_ctx_maps_to_str.cortex_kernel_px));

% %%%%%%%%%%%%%%% 
% %%% Moving median domain_idx to reduce random flipping
% domain_idx_rec = cellfun(@(x) floor(movmedian(x,3)),domain_idx_rec,'uni',false);
% domain_idx = cell2mat(domain_idx_rec);
% %%%%%%%%%%%%%%

kmeans_cluster_mean = reshape(ap.groupfun(@nanmean,maps_flat_cat,[],domain_idx),[U_size,n_k]);

figure;
h = tiledlayout(n_k,1,'tilespacing','none');
for curr_k = 1:n_k
    % nexttile;
    % imagesc(kmeans_centroid(:,:,curr_k));
    % axis image off
    % clim(max(abs(clim)).*[-1,1]);
    % colormap(ap.colormap('PWG'));
    % ap.wf_draw('ccf','k');

    nexttile;
    imagesc(kmeans_cluster_mean(:,:,curr_k));
    axis image off
    clim(max(abs(clim)).*[-1,1]);
    colormap(ap.colormap('PWG'));
    ap.wf_draw('ccf','k');
end

% (to draw points at max position - unused at the moment)
% figure; hold on; axis image; set(gca,'ydir','reverse')
% max_map = imgaussfilt(maps_cat,10).*round(linspace(1,0,size(maps_cat,2)));
% [~,kernel_max_idx_full] = max(max_map,[],[1,2],'linear');
% [kernel_max_y,kernel_max_x,~] = ind2sub(size(maps_cat),kernel_max_idx_full);
% for plot_domain_idx = 1:n_k
%     plot(squeeze(kernel_max_x(domain_idx==plot_domain_idx)),squeeze(kernel_max_y(domain_idx==plot_domain_idx)),'x','linewidth',2.5);
% end
% ap.wf_draw('ccf','k');
% legend(string(1:n_k))


%% Map difference across days

% Maps and explained variance by day for each cluster

map_ld_grp = cell2mat(cellfun(@(x,maps) repelem(x,size(maps,3)*(size(maps,1)>0),1), ...
    num2cell(bhv.days_from_learning), ...
    all_ctx_maps_to_str.cortex_kernel_px,'uni',false));

map_animal_grp = cell2mat(cellfun(@(x,maps) repelem(x,size(maps,3)*(size(maps,1)>0),1), ...
    num2cell(grp2idx(bhv.animal)), ...
    all_ctx_maps_to_str.cortex_kernel_px,'uni',false));

plotday_bins = [-Inf,-2:2,Inf];
map_plotday_grp = discretize(map_ld_grp,plotday_bins);


% Plot maps by day
[map_avg,map_grp] = ap.nestgroupfun({@nanmean,@nanmean}, maps_flat_cat',...
    map_animal_grp,[domain_idx,map_plotday_grp]);

figure;
h = tiledlayout(n_k,length(plotday_bins)-1,'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = 1:length(plotday_bins)-1
        nexttile;
        curr_idx = map_grp(:,1) == curr_k & map_grp(:,2) == curr_day;
        imagesc(reshape(map_avg(curr_idx,:),U_size));
        ap.wf_draw('ccf','k');

        axis off image;
        colormap(ap.colormap('PWG'));
        clim(max(abs(map_avg),[],'all').*[-1,1].*0.6);
    end
end

% Plot explained variance by day
expl_var_avg = ap.nestgroupfun({@nanmean,@nanmean}, expl_var_cat,...
    map_animal_grp,[domain_idx,map_plotday_grp]);
expl_var_sem = ap.nestgroupfun({@nanmean,@AP_sem}, expl_var_cat,...
    map_animal_grp,[domain_idx,map_plotday_grp]);

figure; hold on;
for curr_k = 1:n_k
    curr_idx = map_grp(:,1) == curr_k;
    errorbar(expl_var_avg(curr_idx),expl_var_sem(curr_idx),'linewidth',2);
    ylabel('Explained variance');
end



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

% (smooth and normalize)
psth_t = -0.5:0.001:1;
t_baseline = psth_t < 0;
softnorm = 1;
smooth_window = 100;

unit_psth_cat = cell2mat(cellfun(@(x) smoothdata(...
    (x - nanmean(x(:,t_baseline,:),[2,3]))./(nanmean(x(:,t_baseline,:),[2,3]) + softnorm),...
    2,'gaussian',[smooth_window,0]), ...
    ephys.unit_event_psths,'uni',false,'ErrorHandler',@(varargin) []));

% Make index of recording per unit
unit_animal = cell2mat(cellfun(@(grp,units) repmat(grp,length(units),1), ...
    num2cell(grp2idx(bhv.animal)),ephys.single_unit_idx,'uni',false));

% Make learning day per unit
unit_ld = cell2mat(cellfun(@(grp,units) repmat(grp,length(units),1), ...
    num2cell(bhv.days_from_learning),ephys.single_unit_idx,'uni',false));

plot_day_bins = [-Inf,-2:2,Inf];
unit_plot_days_grp = discretize(max(unit_ld,-inf),plot_day_bins);

% Make kmeans cluster per unit
unit_domain_idx_subset = cell2mat(cellfun(@(depth,domain_idx) domain_idx(depth(~isnan(depth))), ...
    ephys.unit_depth_group,domain_idx_rec,'uni',false));
unit_domain_idx = nan(size(unit_rec_idx));
unit_domain_idx(~isnan(cell2mat(ephys.unit_depth_group))) = unit_domain_idx_subset;


% Group data and plot
group_labels = [unit_rec_idx];
split_labels = [unit_plot_days_grp,unit_domain_idx];

% (frac responsive cells)
% use_units = true(size(unit_single_cat));
% use_units = unit_single_cat;
% use_units = logical(vertcat(ephys.str_msn_idx{:}));
% use_units = logical(vertcat(ephys.str_msn_idx{:})) & unit_resp_cat(:,2);
use_units = logical(vertcat(ephys.str_msn_idx{:})) & unit_single_cat;

figure;
h = tiledlayout(n_k,1);
for curr_k = 1:n_k

    nexttile; hold on;
    set(gca,'ColorOrder',ap.colormap('BKR',3));

    for curr_stim = 1:3
        [unit_group_mean,groups] = ap.nestgroupfun({@mean,@mean}, ...
            +unit_resp_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

        unit_group_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
            +unit_resp_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

        [n_units_ld,n_units_ld_group] = ap.nestgroupfun({@sum,@mean},+use_units,group_labels,split_labels);

        plot_data = groups(:,2) == curr_k;
        errorbar(groups(plot_data,1),unit_group_mean(plot_data,1), ...
            unit_group_sem(plot_data,1),'linewidth',2);
        ylabel('Frac responsive cells');

        % nexttile;
        % plot_data = n_units_ld_group(:,2) == curr_k;
        % plot(n_units_ld_group(plot_data,1),n_units_ld(plot_data),'k');
        % ylabel('N units');
    end
end
linkaxes(h.Children,'xy');

% (psth: average of single units)
curr_stim = 3;
% use_units = true(size(unit_domain_idx));
% use_units = unit_single_cat;
use_units = unit_resp_cat(:,3);
% use_units = logical(vertcat(ephys.str_msn_idx{:}));

[unit_group_mean,groups] = ap.nestgroupfun({@mean,@mean}, ...
    unit_psth_cat(use_units,:,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

unit_group_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    unit_psth_cat(use_units,:,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

figure; 
h = tiledlayout(n_k,max(unit_plot_days_grp),'TileSpacing','none');
for curr_k = unique(groups(:,2))'
    for curr_day = unique(unit_plot_days_grp)'
        nexttile; axis off;
        ap.errorfill([], ...
            unit_group_mean(groups(:,1) == curr_day & groups(:,2) == curr_k,:)', ...
            unit_group_sem(groups(:,1) == curr_day & groups(:,2) == curr_k,:)','k');
    end
end
linkaxes(h.Children,'xy');




%% (from above) Cell type heatmap

% plot_celltype = vertcat(ephys.str_msn_idx{:});% & vertcat(ephys.single_unit_idx{:});
% plot_celltype = vertcat(ephys.str_fsi_idx{:});% & vertcat(ephys.single_unit_idx{:});
plot_celltype = vertcat(ephys.str_tan_idx{:});% & vertcat(ephys.single_unit_idx{:});

plot_domain_idx = [1:2];
plot_stim_idx = 2;

stim_t = psth_t > 0 & psth_t < 0.2;

% Plot grouped days

% unit_ld_prepost = unit_ld >= 0;
plot_day_bins = [-Inf,-2:2,Inf];
unit_plot_day_grp = discretize(unit_ld,plot_day_bins);

figure;
colormap(ap.colormap('BWR',[],2));
h = tiledlayout(2,max(unit_plot_day_grp),'TileIndexing','column','TileSpacing','compact');
for curr_ld = unique(unit_plot_day_grp(~isnan(unit_plot_day_grp)))'
    nexttile;

    curr_data = unit_psth_cat(ismember(unit_domain_idx,plot_domain_idx) & unit_plot_day_grp == curr_ld & ...
        plot_celltype,:,plot_stim_idx);
    [~,sort_idx] = sort(max(curr_data(:,stim_t),[],2),'descend');
    imagesc(psth_t,[],curr_data(sort_idx,:));
    clim([-1,1])
    title(plot_day_bins(curr_ld));

    nexttile;
    curr_data = unit_psth_cat(ismember(unit_domain_idx,plot_domain_idx) & unit_plot_day_grp == curr_ld & ...
        plot_celltype,:,plot_stim_idx);
    ap.errorfill(psth_t,mean(curr_data,1),AP_sem(curr_data,1),'k');
end
linkaxes(h.Children(1:2:end),'xy');
linkaxes(h.Children(2:2:end),'xy'); 


%% (from above) Unit heatmap across stimuli

plot_day_bins = [-Inf,-2:1,Inf];
unit_plot_day_grp = discretize(unit_ld,plot_day_bins);

figure;
h = tiledlayout(1,length(plot_day_bins)-1);

curr_domain_idx = 1;
for curr_day_grp = 1:length(plot_day_bins)-1

    curr_units = ismember([unit_domain_idx,unit_plot_day_grp],[curr_domain_idx,curr_day_grp],'rows') & ...
        logical(vertcat(ephys.str_msn_idx{:}));

    curr_data = unit_psth_cat(curr_units,:,:);

    use_t = psth_t > 0 & psth_t < 0.2;
    [~,sort_idx] = sort(max(curr_data(:,use_t,2),[],2) - max(curr_data(:,use_t,3),[],2));

    nexttile;
    imagesc(movmean(reshape(permute(curr_data(sort_idx,:,:),[2,3,1]),[],length(sort_idx))',[20,1]));

    clim([-3,3])
    colormap(ap.colormap('BWR'));

end





%% (from above) Unit response overlap

% Get max response in window 
use_t = psth_t > 0 & psth_t < 0.2;
unit_psth_max = permute(max(unit_psth_cat(:,use_t,:),[],2),[1,3,2]);

plot_day_bins = [-Inf,-2:1,Inf];
unit_plot_day_grp = discretize(unit_ld,plot_day_bins);

% Scatter (raw)
figure;
h = tiledlayout(1,length(plot_day_bins)-1);

curr_domain_idx = 1;
for curr_day_grp = 1:length(plot_day_bins)-1
    nexttile; axis equal
    curr_units = ismember([unit_domain_idx,unit_plot_day_grp],[curr_domain_idx,curr_day_grp],'rows') & ...
        logical(vertcat(ephys.str_msn_idx{:}));
    plot(unit_psth_max(curr_units,2),unit_psth_max(curr_units,3),'.k');
end
linkaxes(h.Children,'xy');

% Scatter (normalized)
unit_response_norm = unit_psth_max./sum(abs(unit_psth_max),2);

figure;
h = tiledlayout(1,length(plot_day_bins)-1);

curr_domain_idx = 1;
for curr_day_grp = 1:length(plot_day_bins)-1
    nexttile; axis equal
    curr_units = ismember([unit_domain_idx,unit_plot_day_grp],[curr_domain_idx,curr_day_grp],'rows') & ...
        logical(vertcat(ephys.str_msn_idx{:}));
    plot3(unit_response_norm(curr_units,1),unit_response_norm(curr_units,2),unit_response_norm(curr_units,3),'.k');
end
linkaxes(h.Children,'xy');



% Trying out snowflake plot
unit_response_norm = unit_psth_max;
unit_response_relative = abs(unit_response_norm)./max(abs(unit_response_norm),[],2);
% unit_response_relative = min(1,abs(unit_response_norm)./prctile(abs(unit_response_norm),80,'all'));

ternary_tform = [0,1;sqrt(3/4),-0.5;-sqrt(3/4),-0.5];
unit_response_norm_ternary = (ternary_tform\unit_response_relative')';

unit_max_response = max(abs(unit_response_norm),[],2);

response_limits = reshape([eye(3),eye(3)+circshift(eye(3),[0,1])]',3,[])';
ternary_limits = (ternary_tform\response_limits')';

figure;
h = tiledlayout(1,length(plot_day_bins)-1);

curr_domain_idx = 1;
for curr_day_grp = 1:length(plot_day_bins)-1
    curr_data_idx = ismember([unit_domain_idx,unit_plot_day_grp],[curr_domain_idx,curr_day_grp],'rows') & ...
        logical(vertcat(ephys.str_msn_idx{:}));

    % Plot each unit
    size_norm_denom = prctile(unit_max_response(unit_domain_idx == curr_domain_idx),90);
    size_norm = min(unit_max_response(curr_data_idx)./size_norm_denom,1);
    dot_size = size_norm.*20+1;

    response_sign = sign(ap.signed_max(unit_response_norm(curr_data_idx,:),2));
    dot_color = [response_sign == 1,zeros(size(response_sign)),response_sign == -1];

    nexttile; hold on; axis image off;
    patch(ternary_limits(:,1),ternary_limits(:,2),'w','linewidth',1);
    scatter(unit_response_norm_ternary(curr_data_idx,1), ...
        unit_response_norm_ternary(curr_data_idx,2), ...
        dot_size,dot_color,'filled');
    % ap.heatscatter(unit_response_norm_ternary(curr_data_idx,1), ...
    %     unit_response_norm_ternary(curr_data_idx,2), ...
    %     10,1);
end




%% MUA (passive)

animal_groupfun = @mean;

[psth_sum,psth_sum_domain_idx] = cellfun(@(mua,domain_idx) ap.groupfun(@sum,mua,[],[],domain_idx), ...
    ephys.binned_msn_spikes_stim_align,domain_idx_rec,'uni',false);

% (currently psth time is hard coded: save this somewhetre)
psth_t = -0.5:0.001:1;
baseline_t = psth_t < 0;
softnorm = 10;
psth_baseline = cellfun(@(mua) ...
    mean(mua(:,baseline_t,:),[1,2]), ...
    psth_sum,'uni',false,'ErrorHandler',@(varargin) NaN);

psth_norm = cellfun(@(psth,baseline) ...
    (psth-nanmean(psth(:,baseline_t,:),2))./(baseline+softnorm),psth_sum,psth_baseline,'uni',false, ...
    'ErrorHandler',@(varargin) []);

psth_norm_cat_smooth = smoothdata(cell2mat(cellfun(@(x) reshape(permute(x,[1,3,2]),[], ...
    length(psth_t)),psth_norm,'uni',false)),2,'gaussian',[100,0]);

[stim_grp_sq,domain_idx_grp_sq] = cellfun(@(stim,domain_idx) ndgrid(stim,domain_idx), ...
    ephys.trial_stim_values,psth_sum_domain_idx,'uni',false);

animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(grp2idx(ephys.animal)),stim_grp_sq,'uni',false));

ld_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(bhv.days_from_learning),stim_grp_sq,'uni',false));

plot_day_bins = [-Inf,-2:2,Inf];
% plot_day_bins = [-3:3];
plot_days_grp = discretize(max(ld_grp,-inf),plot_day_bins);

stim_grp = cell2mat(cellfun(@(x) reshape(x,[],1),stim_grp_sq,'uni',false));
domain_idx_grp = cell2mat(cellfun(@(x) reshape(x,[],1),domain_idx_grp_sq,'uni',false));

[psth_avg,psth_grp] = ap.nestgroupfun({@mean,animal_groupfun},psth_norm_cat_smooth,animal_grp, ...
    [plot_days_grp,stim_grp,domain_idx_grp]);
psth_sem = ap.nestgroupfun({@mean,@AP_sem},psth_norm_cat_smooth,animal_grp, ...
    [plot_days_grp,stim_grp,domain_idx_grp]);

% (quick and dirty: plot n's)
[psth_n,psth_n_grp] = ap.nestgroupfun({@mean,@length}, ...
    psth_norm_cat_smooth,animal_grp, ...
    [ld_grp,stim_grp,domain_idx_grp]);
figure; hold on;
for curr_k = 1:n_k
    curr_plot_idx = psth_n_grp(:,2)==90 & psth_n_grp(:,3)==curr_k;
    plot(psth_n_grp(curr_plot_idx,1),psth_n(curr_plot_idx));
end
yline(5)
xlabel('Learning day');
ylabel('N');

% (get max: avg within day, max in time, avg agross animals)
[psth_dayavg,psth_dayavg_grp] = ap.groupfun(@mean, ...
    psth_norm_cat_smooth, ...
    [animal_grp,plot_days_grp,stim_grp,domain_idx_grp]);

max_t = psth_t > 0 & psth_t < 0.2;
[psth_max,psth_max_grp] = ap.nestgroupfun({@mean,animal_groupfun}, ...
    max(psth_dayavg(:,max_t),[],2),psth_dayavg_grp(:,1), ...
    psth_dayavg_grp(:,2:end));
psth_max_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    max(psth_dayavg(:,max_t),[],2),psth_dayavg_grp(:,1), ...
    psth_dayavg_grp(:,2:end));

day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

% Plot overlaid PSTH by learning day for stim
figure;
h = tiledlayout(n_k,length(unique(stim_grp)),'TileIndexing','column');
for curr_stim = unique(stim_grp)'
    for curr_k = unique(domain_idx_grp)'
        nexttile; hold on; set(gca,'ColorOrder',day_colormap);
        plot(psth_avg(psth_grp(:,2) == curr_stim & psth_grp(:,3) == curr_k,:)')
    end
end
linkaxes(h.Children,'xy');

% Plot PSTH by learning day for stim
% (color by day)
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

% Plot PSTH by learning day for all stim
stim_plot = [-90,0,90];
stim_color = [0,0,1;0,0,0;1,0,0];
figure('Name','MUA by LD');
h = tiledlayout(n_k,max(psth_grp(:,1)),'TileSpacing','none');
for curr_k = unique(psth_grp(:,3))'
    for curr_day = unique(psth_grp(:,1))'
        nexttile; axis off;
        for curr_stim_idx = 1:length(stim_plot)            
            plot_data = psth_grp(:,1) == curr_day & ...
                psth_grp(:,2) == stim_plot(curr_stim_idx) & ...
                psth_grp(:,3) == curr_k;
            if ~any(plot_data)
                continue
            end
            ap.errorfill(psth_t,psth_avg(plot_data,:),psth_sem(plot_data,:), ...
                stim_color(curr_stim_idx,:));
        end
    end
end
linkaxes(h.Children,'xy');

% Plot PSTH max by learning day
figure;
h = tiledlayout(n_k,1);
for curr_k = unique(domain_idx_grp)'
    nexttile; hold on; set(gca,'ColorOrder',ap.colormap('BKR',3));
    for curr_stim = unique(stim_grp)'
        plot_grps = psth_max_grp(:,2) == curr_stim & ...
            psth_max_grp(:,3) == curr_k;

        errorbar(psth_max(plot_grps),psth_max_sem(plot_grps),'linewidth',2);
    end
    xline(find(plot_day_bins == 0));
    axis padded
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
            plot_days_grp == curr_day & domain_idx_grp == plot_k);
        imagesc(psth_t,[],movmean(psth_norm_cat_smooth(curr_trials,:),[3,1]));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;
        if curr_day == min(plot_days_grp)
            title(animals{curr_animal});
        end
    end
end

%% (from above - yesterday behaviors vs increase)

% Get fraction of fast reaction times each day
n_fast_trials = cellfun(@(x) sum(x <= 0.25),bhv.stim_to_move);

% Get max of avg psth
[psth_avg,psth_grp] = ap.nestgroupfun({@mean,@mean}, ...
    psth_norm_cat_smooth,(1:length(animal_grp))', ...
    [animal_grp,ld_grp,stim_grp,domain_idx_grp]);

max_t = psth_t > 0 & psth_t < 0.2;
psth_avg_tmax = max(psth_avg(:,max_t),[],2);

x = [];
y = [];

figure; hold on; set(gca,'ColorOrder',ap.colormap('tube',length(animals)));
for curr_animal = 1:length(animals)

    use_psth_idx = psth_grp(:,1) == curr_animal & psth_grp(:,3) == 90 & psth_grp(:,4) == 1 & psth_grp(:,2) <= 0;
    use_bhv_idx = grp2idx(bhv.animal) == curr_animal & ismember(bhv.days_from_learning,psth_grp(use_psth_idx,2));

    curr_ld = psth_grp(use_psth_idx,2);
    curr_psth = psth_avg_tmax(use_psth_idx);
    curr_fastrxn = n_fast_trials(use_bhv_idx);

    % psth_diff = diff(curr_psth)./(curr_psth(2:end)+1);
%     psth_diff = curr_psth(2:end);
    psth_diff = curr_psth(2:end)./(curr_psth(1:end-1));
    use_psth_diff = diff(curr_ld) == 1;

    x = vertcat(x,curr_fastrxn(use_psth_diff));
    y = vertcat(y,psth_diff(use_psth_diff));
    
   % plot(curr_fastrxn(1:end-1),curr_psth(2:end),'.','MarkerSize',20)
   % plot(curr_fastrxn(1:end-1),diff(curr_psth),'.','MarkerSize',20)
    plot(curr_fastrxn(use_psth_diff),psth_diff(use_psth_diff),'.','MarkerSize',20)

end
xlabel('% fast reaction');
ylabel('Passive str increase next day');
axis padded






 %% MUA (task)

[psth_sum,psth_sum_domain_idx] = cellfun(@(mua,domain_idx) ap.groupfun(@sum,mua,[],[],domain_idx), ...
    ephys.binned_msn_spikes_stim_align,domain_idx_rec,'uni',false);

% (currently psth time is hard coded: save this somewhere)
psth_t = -0.5:0.001:1;
baseline_t = psth_t < 0;
softnorm = 10;
psth_baseline = cellfun(@(mua) ...
    mean(mua(:,baseline_t,:),[1,2]) + softnorm, ...
    psth_sum,'uni',false,'ErrorHandler',@(varargin) NaN);

% %%%%%%%%%%%%%%%% MOVE ALIGN
% [psth_sum,psth_sum_domain_idx] = cellfun(@(mua,domain_idx) ap.groupfun(@sum,mua,[],[],domain_idx), ...
%     ephys.binned_msn_spikes_move_align,domain_idx_rec,'uni',false);
% %%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%% OUTCOME ALIGN
% [psth_sum,psth_sum_domain_idx] = cellfun(@(mua,domain_idx) ap.groupfun(@sum,mua,[],[],domain_idx), ...
%     ephys.binned_msn_spikes_outcome_align,domain_idx_rec,'uni',false);
% %%%%%%%%%%%%%%%%

psth_norm = cellfun(@(psth,baseline) ...
    (psth-baseline)./baseline,psth_sum,psth_baseline,'uni',false);

psth_norm_cat_smooth = smoothdata(cell2mat(cellfun(@(x) reshape(permute(x,[1,3,2]),[], ...
    length(psth_t)),psth_norm,'uni',false)),2,'gaussian',[100,0]);

[rxn_grp_sq,domain_idx_grp_sq] = cellfun(@(stim,domain_idx) ndgrid(stim,domain_idx), ...
    bhv.stim_to_move,psth_sum_domain_idx,'uni',false);

animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(grp2idx(ephys.animal)),rxn_grp_sq,'uni',false));

training_day = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));
td_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(training_day),rxn_grp_sq,'uni',false));

ld_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(bhv.days_from_learning),rxn_grp_sq,'uni',false));

rxn_grp = cell2mat(cellfun(@(x) reshape(x,[],1),rxn_grp_sq,'uni',false));
domain_idx_grp = cell2mat(cellfun(@(x) reshape(x,[],1),domain_idx_grp_sq,'uni',false));


day_grp = discretize(max(ld_grp,-inf),[-Inf,-3:2,Inf]);

% Plot heatmaps sorted by reaction times
figure; tiledlayout(n_k,max(day_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(day_grp)'
        nexttile;
        curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k);
       
        [~,sort_idx] = sort(rxn_grp(curr_trials));

        imagesc(psth_t,[],movmean(psth_norm_cat_smooth(curr_trials(sort_idx),:),[20,5]));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;
    end
end

% Plot heatmaps sorted by reaction times (animals separately)
animals = unique(ephys.animal,'stable');
plot_k = 1;
figure; tiledlayout(max(day_grp),max(animal_grp), ...
    'TileIndexing','columnmajor','TileSpacing','none');
for curr_animal = 1:max(animal_grp)
    for curr_day = unique(day_grp)'
        nexttile;
        curr_trials = find(animal_grp == curr_animal & day_grp == curr_day & domain_idx_grp == plot_k);
        [~,sort_idx] = sort(rxn_grp(curr_trials));
        imagesc(psth_t,[],movmean(psth_norm_cat_smooth(curr_trials(sort_idx),:),[3,1]));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;
        if curr_day == min(day_grp)
            title(animals{curr_animal});
        end
    end
end

% Choose split type
split_type = 'rxn_bin';
switch split_type
    case 'trial_percentile'
        n_split = 3;
    case 'rxn_percentile'
        n_split = 4;
    case 'rxn_bin'
%         rxn_bins = [0,0.2,0.4,0.6,0.8,1];
%         rxn_bins = [-Inf,0:0.2:0.6,Inf];
%         rxn_bins = [-Inf,0.05,0.15,0.25,0.5,1,Inf];
        rxn_bins = [-Inf,0.25,Inf];
        n_split = length(rxn_bins)-1;
end


% Plot PSTH split
figure('color','w'); h = tiledlayout(n_k,max(day_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(day_grp)'
        nexttile;
        curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k);
        curr_data = psth_norm_cat_smooth(curr_trials,:);

        switch split_type
            case 'trial_percentile'
                split_idx = cell2mat(arrayfun(@(x) ...
                    min(floor(linspace(1,n_split+1,sum(animal_grp(curr_trials)==x))),n_split)', ...
                    unique(animal_grp(curr_trials)),'uni',false));
            case 'rxn_percentile'
%                 % (total)
%                 split_idx_unsorted = min(floor(linspace(1,n_split+1,size(curr_data,1))),n_split)';
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
        curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k);
        curr_data = psth_norm_cat_smooth(curr_trials,:);

        switch split_type
            case 'trial_percentile'
                split_idx = cell2mat(arrayfun(@(x) ...
                    min(floor(linspace(1,n_split+1,sum(animal_grp(curr_trials)==x))),n_split)', ...
                    unique(animal_grp(curr_trials)),'uni',false));
            case 'rxn_percentile'
%                 % (total)
%                 split_idx_unsorted = min(floor(linspace(1,n_split+1,size(curr_data,1))),n_split)';
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

        data_t = psth_t > 0 & psth_t < 0.2;

        % (mean PSTH, max t, mean across animals)
        [psth_mean,psth_mean_grp] = ap.groupfun(@mean,curr_data, ...
            [animal_grp(curr_trials),split_idx]);
        psth_max = max(psth_mean(:,data_t),[],2);
        psth_max_avg = ap.nestgroupfun({@mean,@mean},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
        psth_max_sem = ap.nestgroupfun({@mean,@AP_sem},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
        errorbar(psth_max_avg,psth_max_sem,'k','linewidth',2);

%         % (mean in time window)
%         curr_data_point = mean(curr_data(:,data_t),2);
%         data_avg = ap.nestgroupfun({@mean,@mean},curr_data_point,animal_grp(curr_trials),split_idx);
%         data_sem = ap.nestgroupfun({@mean,@AP_sem},curr_data_point,animal_grp(curr_trials),split_idx);
%         errorbar(data_avg,data_sem,'k','linewidth',2);

        axis padded
    end
end
linkaxes(h.Children,'xy');
set(h.Children,'colororder',ap.colormap('KR',n_split)); 

%% (from above - to trial percentile by reaction group)

day_bins = [-Inf,-2:2,Inf];
day_grp = discretize(max(ld_grp,-inf),day_bins);

% Split by trial percentile within day
n_split = 3;
rxn_bins = [0.3,Inf];

rxn_bin_grp = discretize(rxn_grp,rxn_bins);
[~,~,rec_idx] = unique([animal_grp,day_grp,domain_idx_grp],'rows');

split_idx = nan(size(rec_idx));
for curr_rec = 1:max(rec_idx)
    split_idx(rec_idx == curr_rec) = ...
        min(floor(linspace(1,n_split+1,sum(rec_idx == curr_rec))),n_split)';
end


% % Plot PSTH split for rxn bin 2
% figure('color','w'); h = tiledlayout(n_k,max(day_grp),'TileSpacing','none');
% for curr_k = 1:n_k
%     for curr_day = unique(day_grp)'
%         nexttile;
%         curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k & rxn_bin_grp == 2);
%         plot(psth_t,ap.nestgroupfun({@mean,@mean},psth_norm_cat_smooth(curr_trials,:), ...
%             animal_grp(curr_trials),split_idx(curr_trials))');
%     end
% end
% linkaxes(h.Children,'xy');
% set(h.Children,'colororder',ap.colormap('KR',n_split)); 
% axis(h.Children,'off')



% Plot PSTH max split
figure('color','w'); h = tiledlayout(n_k,max(day_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day = 1:length(day_bins)-1
        nexttile; hold on;
        curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k);
        curr_data = psth_norm_cat_smooth(curr_trials,:);

        data_t = psth_t > 0 & psth_t < 0.2;

        % (mean PSTH, max t, mean across animals)
        [psth_mean,psth_mean_grp] = ap.groupfun(@mean,curr_data, ...
            [animal_grp(curr_trials),split_idx(curr_trials),rxn_bin_grp(curr_trials)]);
        psth_max = max(psth_mean(:,data_t),[],2);
        [psth_max_avg,psth_max_grp] = ap.nestgroupfun({@mean,@mean},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2:3));
        psth_max_sem = ap.nestgroupfun({@mean,@AP_sem},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2:3));

        arrayfun(@(x) errorbar(psth_max_avg(psth_max_grp(:,2) == x), ...
            psth_max_sem(psth_max_grp(:,2) == x),'linewidth',2), ...
            1:max(rxn_bin_grp));

        axis padded off
    end
end
linkaxes(h.Children,'xy');





%% (from above - animal corr for fast rxn and next day)

plot_day_bins = [-Inf,-2:1,Inf];
plot_day_grp = discretize(ld_grp,plot_day_bins);

% Plot PSTH max split

% % (split = reaction percentile in day)
% n_split = 5;
% [~,~,rec_grp] = unique([animal_grp,td_grp],'rows','stable');
% split_idx = cell2mat(arrayfun(@(x) ...
%     discretize(tiedrank(rxn_grp(rec_grp==x)), ...
%     linspace(1,sum(rec_grp==x),n_split+1)), ...
%     unique(rec_grp),'uni',false));

% (split = reaction bins)
rxn_bins = [0,0.25,Inf];
split_idx = discretize(rxn_grp,rxn_bins);

% (mean PSTH, max t, mean across animals)
data_t = psth_t > 0 & psth_t < 0.2;

[psth_mean,psth_mean_grp] = ap.groupfun(@mean,psth_norm_cat_smooth, ...
    [animal_grp,plot_day_grp,split_idx,domain_idx_grp]);

psth_tgrp = mean(psth_mean(:,data_t),2);

plot_k = 1;
bin_x = 1;
bin_y = 2;

figure;
h = tiledlayout(1,max(plot_day_grp)-1);
for curr_day_grp = 1:max(plot_day_grp)-1
    nexttile; hold on;
    arrayfun(@(animal) plot(...
        psth_tgrp(ismember(psth_mean_grp,[animal,curr_day_grp,  bin_x,plot_k],'rows')), ...
        psth_tgrp(ismember(psth_mean_grp,[animal,curr_day_grp+1,bin_y,plot_k],'rows')),'.','MarkerSize',20), ...
        1:length(animals),'ErrorHandler',@(varargin) [],'uni',false);
    xlim(ylim);
    axis square;
    refline(1,0);
end

figure; hold on;
col = copper(max(plot_day_grp)-1);
for curr_day_grp = 1:max(plot_day_grp)-1
    arrayfun(@(animal) plot(...
        psth_tgrp(ismember(psth_mean_grp,[animal,curr_day_grp,  bin_x,plot_k],'rows')), ...
        psth_tgrp(ismember(psth_mean_grp,[animal,curr_day_grp+1,bin_y,plot_k],'rows')),'.', ...
        'MarkerSize',20,'color',col(curr_day_grp,:)), ...
        1:length(animals),'ErrorHandler',@(varargin) [],'uni',false);
end
xlim(ylim);
axis square;
refline(1,0);



% pre-learning bin 1 vs bin 2
% (was packagking data into matrix here)
x = nan(length(animals),max(plot_day_grp),max(split_idx),n_k);
x_ind = sub2ind(size(x),psth_mean_grp(:,1),psth_mean_grp(:,2),psth_mean_grp(:,3),psth_mean_grp(:,4));
x(x_ind) = psth_tgrp;

figure; h = tiledlayout(4,1);

m = x(:,:,:,1);
m(:,:,3) = NaN;
m2 = reshape(permute(m,[3,2,1]),[],length(animals));
nexttile;plot(m2);
title('within-day')
set(gca,'colororder',jet(length(animals)));

mz = cat(3,m(:,1:end-1,1),m(:,2:end,2),m(:,2:end,3));
mz2 = reshape(permute(mz,[3,2,1]),[],length(animals));
nexttile;plot(mz2);
title('bin1 day 1 vs bin2 day 2')
set(gca,'colororder',jet(length(animals)));

mr = cat(3,m(:,1:end-1,1),m(:,2:end,1),m(:,2:end,3));
mr2 = reshape(permute(mr,[3,2,1]),[],length(animals));
nexttile;plot(mr2);
title('bin1 day 1 vs bin1 day 2')
set(gca,'colororder',jet(length(animals)));

ma = cat(3,m(:,1:end-1,2),m(:,2:end,2),m(:,2:end,3));
ma2 = reshape(permute(ma,[3,2,1]),[],length(animals));
nexttile;plot(mr2);
title('bin2 day 1 vs bin2 day 2')
set(gca,'colororder',jet(length(animals)));

linkaxes(h.Children,'xy');

% (as above, but as differences)
m = x(:,:,:,1);

figure; h = tiledlayout(4,1);

nexttile; hold on;
m2 = diff(m,[],3);
plot(m2');
errorbar(nanmean(m2,1),AP_sem(m2,1),'k','linewidth',2);
yline(0);

nexttile; hold on;
m2 = (m(:,2:end,2) - m(:,1:end-1,1));
plot(m2');
errorbar(nanmean(m2,1),AP_sem(m2,1),'k','linewidth',2);
yline(0);

nexttile; hold on;
m2 = (m(:,2:end,1) - m(:,1:end-1,1));
plot(m2');
errorbar(nanmean(m2,1),AP_sem(m2,1),'k','linewidth',2);
yline(0);

nexttile; hold on;
m2 = (m(:,2:end,2) - m(:,1:end-1,2));
plot(m2');
errorbar(nanmean(m2,1),AP_sem(m2,1),'k','linewidth',2);
yline(0);

linkaxes(h.Children,'xy');
axis padded

%% MUA (task, raw vs move-removed)

[psth_sum,psth_sum_domain_idx] = cellfun(@(mua,domain_idx) ap.groupfun(@sum,mua,[],[],domain_idx), ...
    ephys.binned_msn_spikes_stim_align,domain_idx_rec,'uni',false);

% (currently psth time is hard coded: save this somewhere)
psth_t = -0.5:0.001:1;
baseline_t = psth_t < 0;
softnorm = 10;
psth_baseline = cellfun(@(mua) ...
    mean(mua(:,baseline_t,:),[1,2]) + softnorm, ...
    psth_sum,'uni',false,'ErrorHandler',@(varargin) NaN);

% %%%%%%%%%%%%%%%% MOVE ALIGN
% [psth_sum,psth_sum_domain_idx] = cellfun(@(mua,domain_idx) ap.groupfun(@sum,mua,[],[],domain_idx), ...
%     ephys.binned_msn_spikes_move_align,domain_idx_rec,'uni',false);
% %%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%% OUTCOME ALIGN
% [psth_sum,psth_sum_domain_idx] = cellfun(@(mua,domain_idx) ap.groupfun(@sum,mua,[],[],domain_idx), ...
%     ephys.binned_msn_spikes_outcome_align,domain_idx_rec,'uni',false);
% %%%%%%%%%%%%%%%%

psth_norm = cellfun(@(psth,baseline) ...
    (psth-baseline)./baseline,psth_sum,psth_baseline,'uni',false);

psth_norm_cat_smooth = smoothdata(cell2mat(cellfun(@(x) reshape(permute(x,[1,3,2]),[], ...
    length(psth_t)),psth_norm,'uni',false)),2,'gaussian',[100,0]);

[rxn_grp_sq,domain_idx_grp_sq] = cellfun(@(stim,domain_idx) ndgrid(stim,domain_idx), ...
    bhv.stim_to_move,psth_sum_domain_idx,'uni',false);

animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(grp2idx(ephys.animal)),rxn_grp_sq,'uni',false));

training_day = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));
td_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(training_day),rxn_grp_sq,'uni',false));

ld_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(bhv.days_from_learning),rxn_grp_sq,'uni',false));

rxn_grp = cell2mat(cellfun(@(x) reshape(x,[],1),rxn_grp_sq,'uni',false));
domain_idx_grp = cell2mat(cellfun(@(x) reshape(x,[],1),domain_idx_grp_sq,'uni',false));

day_grp = discretize(max(ld_grp,-inf),[-Inf,-2:2,Inf]);


% Plot heatmaps sorted by reaction times
figure; tiledlayout(n_k,max(day_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(day_grp)'
        nexttile;
        curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k);
       
        [~,sort_idx] = sort(rxn_grp(curr_trials));

        imagesc(psth_t,[],movmean(psth_norm_cat_smooth(curr_trials(sort_idx),:),[20,5]));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;
    end
end

% Plot heatmaps sorted by reaction times (NaN-out move times)
figure; tiledlayout(n_k,max(day_grp),'TileSpacing','none');
move_leeway = 0.1;
for curr_k = 1:n_k
    for curr_day = unique(day_grp)'
        nexttile;
        curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k);      
        [~,sort_idx] = sort(rxn_grp(curr_trials));

        curr_data = psth_norm_cat_smooth(curr_trials(sort_idx),:).* ...
            AP_nanout(psth_t > rxn_grp(curr_trials(sort_idx))-move_leeway);

        imagesc(psth_t,[],movmean(curr_data,[20,5]));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;
    end
end

% Plot average activity in trial (NaN-out movement times)
figure; h = tiledlayout(n_k,max(day_grp),'TileSpacing','none');
move_leeway = 0.1;
for curr_k = 1:n_k
    for curr_day = unique(day_grp)'
        nexttile; hold on;
        curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k);      
        
        % Average all data
        curr_data = psth_norm_cat_smooth(curr_trials,:);
        curr_data_no_move = psth_norm_cat_smooth(curr_trials,:).* ...
            AP_nanout(psth_t > rxn_grp(curr_trials)-move_leeway);

        curr_data_mean = mean(ap.groupfun(@nanmean,curr_data,animal_grp(curr_trials)),1);
        curr_data_sem = AP_sem(ap.groupfun(@nanmean,curr_data,animal_grp(curr_trials)),1);

        % Average data with movement removede
        % (only plot no-move timepoints with at least N trials and N animals)
        min_nomove_trials = 10;
        min_noanimals_trials = 5;

        curr_data_no_move_n = ap.groupfun(@(x) sum(~isnan(x)),curr_data_no_move,animal_grp(curr_trials));
        curr_data_no_move_animal_mean = ap.groupfun(@nanmean,curr_data_no_move,animal_grp(curr_trials));
        curr_data_no_move_animal_mean_mindata = curr_data_no_move_animal_mean.* ...
            AP_nanout(curr_data_no_move_n < min_nomove_trials).* ...
            AP_nanout(sum(curr_data_no_move_n >= min_nomove_trials) < min_noanimals_trials);
        curr_data_no_move_mean = nanmean(curr_data_no_move_animal_mean_mindata,1);
        curr_data_no_move_sem = AP_sem(curr_data_no_move_animal_mean_mindata,1);

        % Plot
        ap.errorfill(psth_t,curr_data_mean,curr_data_sem,'k',0.3,false,2);
        % ap.errorfill(psth_t,curr_data_no_move_mean,curr_data_no_move_sem,'k',[],[],2);
        plot(psth_t,curr_data_no_move_mean,'k','linewidth',2);
        axis off;
    end
end
linkaxes(h.Children,'xy');

%% MUA (task pre/post, by depth)

n_depths = 6;

depth_rec = cellfun(@(x) discretize(1:size(x,3),n_depths)', ...
    all_ctx_maps_to_str.cortex_kernel_px,'uni',false);

[psth_stim_sum,psth_sum_domain_idx] = cellfun(@(mua,mua_group) ap.groupfun(@sum,mua,[],[],mua_group), ...
    ephys.binned_msn_spikes_stim_align,depth_rec,'uni',false);

psth_move_sum = cellfun(@(mua,mua_group) ap.groupfun(@sum,mua,[],[],mua_group), ...
    ephys.binned_msn_spikes_move_align,depth_rec,'uni',false);

psth_outcome_sum = cellfun(@(mua,mua_group) ap.groupfun(@sum,mua,[],[],mua_group), ...
    ephys.binned_msn_spikes_outcome_align,depth_rec,'uni',false);

% (currently psth time is hard coded: save this somewhere)
psth_t = -0.5:0.001:1;
baseline_t = psth_t < 0;
softnorm = 10;
psth_baseline = cellfun(@(mua) ...
    mean(mua(:,baseline_t,:),[1,2]) + softnorm, ...
    psth_stim_sum,'uni',false,'ErrorHandler',@(varargin) NaN);

psth_norm_smooth_reshape_fcn = @(psth,baseline) ...
    reshape(permute(smoothdata((psth-baseline)./baseline,2, ...
    'gaussian',[100,0]),[1,3,2]),[],length(psth_t));

psth_stim = cell2mat(cellfun(psth_norm_smooth_reshape_fcn,psth_stim_sum,psth_baseline,'uni',false));
psth_move = cell2mat(cellfun(psth_norm_smooth_reshape_fcn,psth_move_sum,psth_baseline,'uni',false));
psth_outcome = cell2mat(cellfun(psth_norm_smooth_reshape_fcn,psth_outcome_sum,psth_baseline,'uni',false));

[rxn_grp_sq,domain_idx_grp_sq] = cellfun(@(stim,domain_idx) ndgrid(stim,domain_idx), ...
    bhv.stim_to_move,psth_sum_domain_idx,'uni',false);

animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(grp2idx(ephys.animal)),rxn_grp_sq,'uni',false));

training_day = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));
td_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(training_day),rxn_grp_sq,'uni',false));

ld_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(bhv.days_from_learning),rxn_grp_sq,'uni',false));

rxn_grp = cell2mat(cellfun(@(x) reshape(x,[],1),rxn_grp_sq,'uni',false));
domain_idx_grp = cell2mat(cellfun(@(x) reshape(x,[],1),domain_idx_grp_sq,'uni',false));

day_grp = discretize(max(ld_grp,-inf),[-Inf,-0,Inf]);

% Plot average activity in trial (pre/post learning)
stim_x = [-0.2,0.3];
move_x = [0,0.4];
outcome_x = [-0.1,0.5];

figure; h = tiledlayout(n_depths,3,'TileSpacing','tight');
for curr_depth = 1:n_depths

    curr_trials = domain_idx_grp == curr_depth;

    % Stim
    nexttile; hold on; axis off;
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean},psth_stim(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem},psth_stim(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(stim_x);

    % Move
    nexttile; hold on; axis off;
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean},psth_move(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem},psth_move(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(move_x);

    % Outcome
    nexttile; hold on; axis off;
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean},psth_outcome(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem},psth_outcome(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(outcome_x);

end
% (link all y, and x of same-alignment)
linkaxes(h.Children,'y');
for ax = 1:3
    linkaxes(h.Children(ax:3:end),'x');
end
% (set same data aspect to have same x-size)
[h.Children.DataAspectRatio] = deal(min(vertcat(h.Children.DataAspectRatio),[],1));

%% MUA (task pre/post, by cluster)

[psth_stim_sum,psth_sum_domain_idx] = cellfun(@(mua,mua_group) ap.groupfun(@sum,mua,[],[],mua_group), ...
    ephys.binned_msn_spikes_stim_align,domain_idx_rec,'uni',false);

psth_move_sum = cellfun(@(mua,mua_group) ap.groupfun(@sum,mua,[],[],mua_group), ...
    ephys.binned_msn_spikes_move_align,domain_idx_rec,'uni',false);

psth_outcome_sum = cellfun(@(mua,mua_group) ap.groupfun(@sum,mua,[],[],mua_group), ...
    ephys.binned_msn_spikes_outcome_align,domain_idx_rec,'uni',false);

% (currently psth time is hard coded: save this somewhere)
psth_t = -0.5:0.001:1;
baseline_t = psth_t < 0;
softnorm = 10;
psth_baseline = cellfun(@(mua) ...
    mean(mua(:,baseline_t,:),[1,2]) + softnorm, ...
    psth_stim_sum,'uni',false,'ErrorHandler',@(varargin) NaN);

psth_norm_smooth_reshape_fcn = @(psth,baseline) ...
    reshape(permute(smoothdata((psth-baseline)./baseline,2, ...
    'gaussian',[100,0]),[1,3,2]),[],length(psth_t));

psth_stim = cell2mat(cellfun(psth_norm_smooth_reshape_fcn,psth_stim_sum,psth_baseline,'uni',false));
psth_move = cell2mat(cellfun(psth_norm_smooth_reshape_fcn,psth_move_sum,psth_baseline,'uni',false));
psth_outcome = cell2mat(cellfun(psth_norm_smooth_reshape_fcn,psth_outcome_sum,psth_baseline,'uni',false));

[rxn_grp_sq,domain_idx_grp_sq] = cellfun(@(stim,domain_idx) ndgrid(stim,domain_idx), ...
    bhv.stim_to_move,psth_sum_domain_idx,'uni',false);

animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(grp2idx(ephys.animal)),rxn_grp_sq,'uni',false));

training_day = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));
td_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(training_day),rxn_grp_sq,'uni',false));

ld_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(bhv.days_from_learning),rxn_grp_sq,'uni',false));

rxn_grp = cell2mat(cellfun(@(x) reshape(x,[],1),rxn_grp_sq,'uni',false));
domain_idx_grp = cell2mat(cellfun(@(x) reshape(x,[],1),domain_idx_grp_sq,'uni',false));

day_grp = discretize(max(ld_grp,-inf),[-Inf,-0,Inf]);

% Plot average activity in trial (pre/post learning)
stim_x = [-0.2,0.3];
move_x = [0,0.4];
outcome_x = [-0.1,0.5];

figure; h = tiledlayout(n_k,3,'TileSpacing','tight');
for curr_depth = 1:n_k

    curr_trials = domain_idx_grp == curr_depth;

    % Stim
    nexttile; hold on; axis off;
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean},psth_stim(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem},psth_stim(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(stim_x);

    % Move
    nexttile; hold on; axis off;
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean},psth_move(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem},psth_move(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(move_x);

    % Outcome
    nexttile; hold on; axis off;
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean},psth_outcome(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem},psth_outcome(curr_trials,:),animal_grp(curr_trials),day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(outcome_x);

end
% (link all y, and x of same-alignment)
linkaxes(h.Children,'y');
for ax = 1:3
    linkaxes(h.Children(ax:3:end),'x');
end
% (set same data aspect to have same x-size)
[h.Children.DataAspectRatio] = deal(min(vertcat(h.Children.DataAspectRatio),[],1));



%% MUA (task, iti movements)

[psth_sum,psth_sum_domain_idx] = cellfun(@(mua,domain_idx) ap.groupfun(@sum,mua,[],[],domain_idx), ...
    ephys.binned_msn_spikes_itimove_align,domain_idx_rec,'uni',false);

% (currently psth time is hard coded: save this somewhere)
psth_t = -0.5:0.001:1;
baseline_t = psth_t < 0;
softnorm = 10;
psth_baseline = cellfun(@(mua) ...=
    mean(mua(:,baseline_t,:),[1,2]) + softnorm, ...
    psth_sum,'uni',false,'ErrorHandler',@(varargin) NaN);

psth_norm = cellfun(@(psth,baseline) ...
    (psth-baseline)./baseline,psth_sum,psth_baseline,'uni',false);

psth_norm_cat_smooth = smoothdata(cell2mat(cellfun(@(x) reshape(permute(x,[1,3,2]),[], ...
    length(psth_t)),psth_norm,'uni',false)),2,'gaussian',[100,0]);

% (for iti movements)
% (ephys.iti_fastmove_maxvel or ephys.iti_fastmove_tduration)
[rxn_grp_sq,domain_idx_grp_sq] = cellfun(@(move_stat,domain_idx) ndgrid(move_stat,domain_idx), ...
    ephys.iti_fastmove_maxvel,psth_sum_domain_idx,'uni',false);

animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(grp2idx(ephys.animal)),rxn_grp_sq,'uni',false));

training_day = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));
td_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(training_day),rxn_grp_sq,'uni',false));

ld_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(bhv.days_from_learning),rxn_grp_sq,'uni',false));

rxn_grp = cell2mat(cellfun(@(x) reshape(x,[],1),rxn_grp_sq,'uni',false));
domain_idx_grp = cell2mat(cellfun(@(x) reshape(x,[],1),domain_idx_grp_sq,'uni',false));


day_grp = discretize(max(ld_grp,-inf),[-Inf,-3:2,Inf]);

% Plot heatmaps sorted by reaction times
figure; tiledlayout(n_k,max(day_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(day_grp)'
        nexttile;
        curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k);
       
        [~,sort_idx] = sort(rxn_grp(curr_trials));

        imagesc(psth_t,[],movmean(psth_norm_cat_smooth(curr_trials(sort_idx),:),[20,5]));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;
    end
end



% Choose split type
split_type = 'iti_stat_percentile';
switch split_type
    case 'iti_stat_percentile'
        n_split = 1;
end


% Plot PSTH split
figure('color','w'); h = tiledlayout(n_k,max(day_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(day_grp)'
        nexttile;
        curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k);
        curr_data = psth_norm_cat_smooth(curr_trials,:);

        switch split_type
            case 'iti_stat_percentile'
%                 % (total)
%                 split_idx_unsorted = min(floor(linspace(1,n_split+1,size(curr_data,1))),n_split)';
%                 [~,sort_idx] = sort(rxn_grp(curr_trials));
%                 [~,rxn_rank] = sort(sort_idx);
%                 split_idx = split_idx_unsorted(rxn_rank);

                % (by animal)
                split_idx = cell2mat(arrayfun(@(x) ...
                    discretize(tiedrank(rxn_grp(intersect(curr_trials,find(animal_grp == x)))), ...
                    linspace(1,length(intersect(curr_trials,find(animal_grp == x))),n_split+1)), ...
                    unique(animal_grp(curr_trials)),'uni',false));
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
        curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k);
        curr_data = psth_norm_cat_smooth(curr_trials,:);

        switch split_type
            case 'iti_stat_percentile'
%                 % (total)
%                 split_idx_unsorted = min(floor(linspace(1,n_split+1,size(curr_data,1))),n_split)';
%                 [~,sort_idx] = sort(rxn_grp(curr_trials));
%                 [~,rxn_rank] = sort(sort_idx);
%                 split_idx = split_idx_unsorted(rxn_rank);

                % (by animal)
                split_idx = cell2mat(arrayfun(@(x) ...
                    discretize(tiedrank(rxn_grp(intersect(curr_trials,find(animal_grp == x)))), ...
                    linspace(1,length(intersect(curr_trials,find(animal_grp == x))),n_split+1)), ...
                    unique(animal_grp(curr_trials)),'uni',false));
        end

        data_t = psth_t > 0 & psth_t < 0.2;

        % (mean PSTH, max t, mean across animals)
        [psth_mean,psth_mean_grp] = ap.groupfun(@mean,curr_data, ...
            [animal_grp(curr_trials),split_idx]);
        psth_max = max(psth_mean(:,data_t),[],2);
        psth_max_avg = ap.nestgroupfun({@mean,@mean},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
        psth_max_sem = ap.nestgroupfun({@mean,@AP_sem},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2));
        errorbar(psth_max_avg,psth_max_sem,'k','linewidth',2);

%         % (mean in time window)
%         curr_data_point = mean(curr_data(:,data_t),2);
%         data_avg = ap.nestgroupfun({@mean,@mean},curr_data_point,animal_grp(curr_trials),split_idx);
%         data_sem = ap.nestgroupfun({@mean,@AP_sem},curr_data_point,animal_grp(curr_trials),split_idx);
%         errorbar(data_avg,data_sem,'k','linewidth',2);

        axis padded
    end
end
linkaxes(h.Children,'xy');
set(h.Children,'colororder',ap.colormap('KR',n_split)); 


%% MUA (task, no-stim)

[psth_sum,psth_sum_domain_idx] = cellfun(@(mua,domain_idx) ap.groupfun(@sum,mua,[],[],domain_idx), ...
    ephys.binned_spikes_stim_align,domain_idx_rec,'uni',false);

% (currently psth time is hard coded: save this somewhere)
psth_t = -0.5:0.001:1;
baseline_t = psth_t < -0.3;
softnorm = 10;
psth_baseline = cellfun(@(mua) ...
    mean(mua(:,baseline_t,:),[1,2]) + softnorm, ...
    psth_sum,'uni',false,'ErrorHandler',@(varargin) NaN);

% %%%%%%%%%%%%%%%% MOVE ALIGN
% [psth_sum,psth_sum_domain_idx] = cellfun(@(mua,domain_idx) ap.groupfun(@sum,mua,[],[],domain_idx), ...
%     ephys.binned_spikes_move_align,domain_idx_rec,'uni',false);
% %%%%%%%%%%%%%%%%

psth_norm = cellfun(@(psth,baseline) ...
    (psth-baseline)./baseline,psth_sum,psth_baseline,'uni',false);

psth_norm_cat_smooth = smoothdata(cell2mat(cellfun(@(x) reshape(permute(x,[1,3,2]),[], ...
    length(psth_t)),psth_norm,'uni',false)),2,'gaussian',[30,0]);
[rxn_grp_sq,domain_idx_grp_sq] = cellfun(@(stim,domain_idx) ndgrid(stim,domain_idx), ...
    bhv.stim_to_move,psth_sum_domain_idx,'uni',false);

opacity_grp_sq = cellfun(@(stim,domain_idx) ndgrid(stim,domain_idx), ...
    ephys.trial_opacity,psth_sum_domain_idx,'uni',false);

animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(grp2idx(ephys.animal)),rxn_grp_sq,'uni',false));

training_day = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));
td_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(training_day),rxn_grp_sq,'uni',false));

ld_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,numel(grp),1),    ...
    num2cell(bhv.days_from_learning),rxn_grp_sq,'uni',false));

rxn_grp = cell2mat(cellfun(@(x) reshape(x,[],1),rxn_grp_sq,'uni',false));
domain_idx_grp = cell2mat(cellfun(@(x) reshape(x,[],1),domain_idx_grp_sq,'uni',false));

opacity_grp = cell2mat(cellfun(@(x) reshape(x,[],1),opacity_grp_sq,'uni',false));

day_grp = discretize(max(ld_grp,-inf),[-Inf,-2:1,Inf]);

% Plot heatmaps sorted by reaction times
figure; tiledlayout(n_k,2,'TileSpacing','none','TileIndexing','column');
for curr_opacity = 0:1
    for curr_k = 1:n_k
        for curr_day = unique(day_grp)'
            nexttile;
            curr_trials = find(td_grp >= 2 & domain_idx_grp == curr_k & opacity_grp == curr_opacity);

            [~,sort_idx] = sort(rxn_grp(curr_trials));
            imagesc(psth_t,[],movmean(psth_norm_cat_smooth(curr_trials(sort_idx),:),[10,5]));
            colormap(ap.colormap('WK'));
            clim([0,3]);
            axis off;
        end
    end
end


% Choose split type
split_type = 'rxn_bin';
switch split_type
    case 'trial_percentile'
        n_split = 4;
    case 'rxn_percentile'
        n_split = 4;
    case 'rxn_bin'
%         rxn_bins = [0,0.2,0.4,0.6,0.8,1];
%         rxn_bins = [-Inf,0:0.2:0.6,Inf];
%         rxn_bins = [-Inf,0.05,0.15,0.25,0.5,1,Inf];
        rxn_bins = [0,0.5,Inf];
        n_split = length(rxn_bins)-1;
end


% Plot PSTH split
for curr_opacity = 0:1
    figure('color','w'); h = tiledlayout(n_k,max(day_grp),'TileSpacing','none');
    for curr_k = 1:n_k
        for curr_day = unique(day_grp)'
            nexttile;
            curr_trials = find(day_grp == curr_day & domain_idx_grp == curr_k & opacity_grp == curr_opacity);
            curr_data = psth_norm_cat_smooth(curr_trials,:);

            switch split_type
                case 'trial_percentile'
                    split_idx = cell2mat(arrayfun(@(x) ...
                        min(floor(linspace(1,n_split+1,sum(animal_grp(curr_trials)==x))),n_split)', ...
                        unique(animal_grp(curr_trials)),'uni',false));
                case 'rxn_percentile'
                    %                 % (total)
                    %                 split_idx_unsorted = min(floor(linspace(1,n_split+1,size(curr_data,1))),n_split)';
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
end




%% Widefield (passive)

plot_stim = 90;

% Create group indicies
wf_t = wf.wf_stim_time{1};

wf_animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,length(grp),1),    ...
    num2cell(grp2idx(wf.animal)),wf.trial_stim_values,'uni',false));

wf_ld_grp = cell2mat(cellfun(@(ld,grp) repmat(ld,length(grp),1),    ...
    num2cell(bhv.days_from_learning),wf.trial_stim_values,'uni',false));

plot_day_bins = [-Inf,-2:2,Inf];
wf_plot_days_grp = discretize(max(wf_ld_grp,-inf),plot_day_bins);


% Average and plot widefield
[wf_avg,wf_avg_grp] = ap.nestgroupfun({@mean,@mean},cell2mat(wf.V_no_move_stim_align), ...
    wf_animal_grp,[wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

curr_data_idx = wf_avg_grp(:,2) == plot_stim;

curr_data_px = plab.wf.svd2px(U_master,permute(wf_avg(curr_data_idx,:,:),[3,2,1]));

ap.imscroll(curr_data_px - nanmean(curr_data_px(:,:,wf_t<0,:),3),wf_t);
colormap(ap.colormap('PWG',[],1.5));
clim(max(curr_data_px(:)).*[-1,1]*0.6);
axis image;
ap.wf_draw('ccf','k');

% Widefield max
stim_t = wf_t > 0 & wf_t < 0.2;
[wf_avg,wf_avg_grp] = ap.groupfun(@mean,cell2mat(wf.V_no_move_stim_align), ...
    [wf_animal_grp,wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

wf_max_px = permute(max(plab.wf.svd2px(U_master,permute(wf_avg(:,stim_t,:),[3,2,1])),[],3),[1,2,4,3]);
[wf_max_px_avg,wf_max_px_avg_grp] = ap.nestgroupfun({@mean,@mean},reshape(wf_max_px,[],size(wf_max_px,3))', ...
    wf_avg_grp(:,1),wf_avg_grp(:,2:3));

figure;
h = tiledlayout(3,max(wf_plot_days_grp),'TileIndexing','ColumnMajor','TileSpacing','none');
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


% Plot ROIs by striatum cluster
% % (weighted average)
% striatum_wf_roi = max(kmeans_cluster_mean,0)./max(kmeans_cluster_mean,[],[1,2]);
% (thresholded spot)

kmeans_cluster_mean_gauss = imgaussfilt(kmeans_cluster_mean,10);
[~,m] = max(kmeans_cluster_mean_gauss,[],[1,2],'linear');
[mr,mc] = ind2sub(size(U_master,[1,2]),m);
striatum_wf_roi = kmeans_cluster_mean_gauss > prctile(kmeans_cluster_mean_gauss,100,[1,2])*0.50;
striatum_wf_roi(:,round(size(striatum_wf_roi,2)/2):end,:) = false;
for k = 1:size(striatum_wf_roi,3)
    [~,m] = max(imgaussfilt(kmeans_cluster_mean_gauss(:,:,k),10).* ...
        round(linspace(1,0,size(kmeans_cluster_mean_gauss,2))),[],[1,2],'linear');
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

% Get ROI activity
wf_striatum_roi = permute(ap.wf_roi(U_master, ...
    permute(cell2mat(wf.V_no_move_stim_align),[3,2,1]),[],[],striatum_wf_roi),[3,2,1]);

% (baseline-subtract)
baseline_t = wf_t >= -0.1 & wf_t <= 0.05;
wf_striatum_roi = wf_striatum_roi - nanmean(wf_striatum_roi(:,baseline_t,:),2);

[wf_striatum_roi_avg,wf_striatum_roi_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi, ...
    wf_animal_grp,[wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

wf_striatum_roi_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi, ...
    wf_animal_grp,[wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

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
linkaxes(h.Children,'xy');t

% Plot average in ROI by day
stim_t = wf_t >= 0 & wf_t <= 0.2;
wf_striatum_roi_tavg = nanmean(wf_striatum_roi(:,stim_t,:),2);

[wf_striatum_roi_tavg_avg,wf_striatum_roi_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi_tavg, ...
    wf_animal_grp,[wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

wf_striatum_roi_tavg_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi_tavg, ...
    wf_animal_grp,[wf_plot_days_grp,cell2mat(wf.trial_stim_values)]);

figure;
h = tiledlayout(n_k,1,'TileSpacing','none');
for curr_k = 1:size(striatum_wf_roi,3)
    nexttile; hold on;
    set(gca,'ColorOrder',ap.colormap('BKR',3));
    for curr_stim = unique(unique(wf_striatum_roi_avg_grp(:,2)))'
        curr_data_idx = wf_striatum_roi_avg_grp(:,2) == curr_stim;

        errorbar(wf_striatum_roi_tavg_avg(curr_data_idx,:,curr_k), ...
            wf_striatum_roi_tavg_sem(curr_data_idx,:,curr_k),'linewidth',2);
        ylabel('t avg');
        axis padded
    end
end
linkaxes(h.Children,'xy');


% Plot max in ROI by day
stim_t = wf_t >= 0 & wf_t <= 0.20;

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

figure;
h = tiledlayout(n_k,1,'TileSpacing','none');
for curr_k = 1:size(striatum_wf_roi,3)
    nexttile; hold on;
    set(gca,'ColorOrder',ap.colormap('BKR',3));
    for curr_stim = unique(unique(wf_striatum_roi_max_avg_grp(:,2)))'
        curr_data_idx = wf_striatum_roi_max_avg_grp(:,2) == curr_stim;

        errorbar(wf_striatum_roi_max_avg(curr_data_idx,:,curr_k), ...
            wf_striatum_roi_max_sem(curr_data_idx,:,curr_k),'linewidth',2);
        ylabel('t max');
        axis padded
    end
end
linkaxes(h.Children,'xy');


%% Widefield (task)

% Create group indicies
wf_t = wf.wf_stim_time{1};

wf_animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,size(grp,1),1),    ...
    num2cell(grp2idx(wf.animal)),wf.V_stim_align,'uni',false));

wf_ld_grp = cell2mat(cellfun(@(ld,grp) repmat(ld,size(grp,1),1),    ...
    num2cell(bhv.days_from_learning),wf.V_stim_align,'uni',false));

plot_day_bins = [-Inf,-2:2,Inf];
wf_plot_days_grp = discretize(max(wf_ld_grp,-inf),plot_day_bins);

wf_rxn_grp = cell2mat(bhv.stim_to_move);


% Make ROIs from striatum maps
% % (weighted average)
% striatum_wf_roi = max(kmeans_cluster_mean,0)./max(kmeans_cluster_mean,[],[1,2]);
% (thresholded spot)
[~,m] = max(kmeans_cluster_mean,[],[1,2],'linear');
[mr,mc] = ind2sub(size(U_master,[1,2]),m);
striatum_wf_roi = kmeans_cluster_mean > max(kmeans_cluster_mean,[],[1,2])*0.5;
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
    axis image off; ap.wf_draw('ccf','k');
    colormap(ap.colormap('WG'));
end

% Get ROI activity
wf_striatum_roi = permute(ap.wf_roi(U_master, ...
    permute(cell2mat(wf.V_stim_align),[3,2,1]),[],[],striatum_wf_roi),[3,2,1]);

% (baseline-subtract)
wf_striatum_roi = wf_striatum_roi - nanmean(wf_striatum_roi(:,wf_t < -0.3,:),2);

% Plot heatmaps sorted by reaction times
figure; tiledlayout(n_k,max(wf_plot_days_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(wf_plot_days_grp)'
        nexttile;
        curr_trials = find(wf_plot_days_grp == curr_day);
       
        [~,sort_idx] = sort(wf_rxn_grp(curr_trials));
        imagesc(wf_t,[],movmean(wf_striatum_roi(curr_trials(sort_idx),:,curr_k),[20,5]));
        colormap(ap.colormap('PWG'));
        clim(max(abs(clim)).*[-1,1]);
        axis off;
    end
end

% Choose split type
split_type = 'rxn_bin';
switch split_type
    case 'percentile'
        n_split = 5;
    case 'rxn_percentile'
        n_split = 5;
    case 'rxn_bin'
%         rxn_bins = [0,0.2,0.4,0.6,0.8,1];
%         rxn_bins = [-Inf,0.05,0.15,0.25,0.5,1,Inf];
        rxn_bins = [0,0.25,Inf];
        n_split = length(rxn_bins)-1;
end

% Plot ROI split
figure('color','w'); h = tiledlayout(n_k,max(wf_plot_days_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(wf_plot_days_grp)'
        nexttile;
        curr_trials = find(wf_plot_days_grp == curr_day);
        curr_data = wf_striatum_roi(curr_trials,:,curr_k);

        switch split_type
            case 'percentile'
                split_idx = cell2mat(arrayfun(@(x) ...
                    min(floor(linspace(1,n_split+1,sum(wf_animal_grp(curr_trials)==x))),n_split)', ...                
                unique(wf_animal_grp(curr_trials)),'uni',false));
            case 'rxn_percentile'
%                 % (total)
%                 split_idx_unsorted = min(floor(linspace(1,n_split+1,size(curr_data,1))),n_split)';
%                 [~,sort_idx] = sort(rxn_grp(curr_trials));
%                 [~,rxn_rank] = sort(sort_idx);
%                 split_idx = split_idx_unsorted(rxn_rank);

                % (by animal)
                split_idx = cell2mat(arrayfun(@(x) ...
                    discretize(tiedrank(wf_rxn_grp(intersect(curr_trials,find(wf_animal_grp == x)))), ...
                    linspace(1,length(intersect(curr_trials,find(wf_animal_grp == x))),n_split+1)), ...
                    unique(wf_animal_grp(curr_trials)),'uni',false));
            case 'rxn_bin'
                split_idx = discretize(wf_rxn_grp(curr_trials),rxn_bins);
        end

        curr_trial_group_data = ap.nestgroupfun({@mean,@mean},curr_data,wf_animal_grp(curr_trials),split_idx);
        plot(wf_t,curr_trial_group_data');

    end
end
linkaxes(h.Children,'xy');
set(h.Children,'colororder',ap.colormap('KR',n_split)); 
axis(h.Children,'off');


% Plot ROI max split
figure('color','w'); h = tiledlayout(n_k,max(wf_plot_days_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day = unique(wf_plot_days_grp)'
        nexttile;
        curr_trials = find(wf_plot_days_grp == curr_day);
        curr_data = wf_striatum_roi(curr_trials,:,curr_k);

        switch split_type
            case 'percentile'
                split_idx = cell2mat(arrayfun(@(x) ...
                    min(floor(linspace(1,n_split+1,sum(wf_animal_grp(curr_trials)==x))),n_split)', ...                
                unique(wf_animal_grp(curr_trials)),'uni',false));
            case 'rxn_percentile'
%                 % (total)
%                 split_idx_unsorted = min(floor(linspace(1,n_split+1,size(curr_data,1))),n_split)';
%                 [~,sort_idx] = sort(rxn_grp(curr_trials));
%                 [~,rxn_rank] = sort(sort_idx);
%                 split_idx = split_idx_unsorted(rxn_rank);

                % (by animal)
                split_idx = cell2mat(arrayfun(@(x) ...
                    discretize(tiedrank(wf_rxn_grp(intersect(curr_trials,find(wf_animal_grp == x)))), ...
                    linspace(1,length(intersect(curr_trials,find(wf_animal_grp == x))),n_split+1)), ...
                    unique(wf_animal_grp(curr_trials)),'uni',false));
            case 'rxn_bin'
                split_idx = discretize(wf_rxn_grp(curr_trials),rxn_bins);
        end

        data_t = wf_t > 0 & wf_t < 0.2;

        % (mean PSTH, max t, mean across animals)
        [wf_mean,wf_mean_grp] = ap.groupfun(@mean,curr_data,[wf_animal_grp(curr_trials),split_idx]);
        wf_max = max(wf_mean(:,data_t),[],2);
        wf_max_avg = ap.nestgroupfun({@mean,@mean},wf_max,wf_mean_grp(:,1),wf_mean_grp(:,2));
        wf_max_sem = ap.nestgroupfun({@mean,@AP_sem},wf_max,wf_mean_grp(:,1),wf_mean_grp(:,2));
        errorbar(wf_max_avg,wf_max_sem,'color',[0,0.7,0],'linewidth',2);

        % % (mean in time window)
        % curr_data_point = mean(curr_data(:,data_t),2);
        % data_avg = ap.nestgroupfun({@mean,@mean},curr_data_point,wf_animal_grp(curr_trials),split_idx);
        % data_sem = ap.nestgroupfun({@mean,@AP_sem},curr_data_point,wf_animal_grp(curr_trials),split_idx);
        % errorbar(data_avg,data_sem,'k','linewidth',2);

        axis padded off
    end
end
linkaxes(h.Children,'xy');
set(h.Children,'colororder',ap.colormap('KR',n_split)); 

%% (from above - to trial percentile by reaction group)
% (copied from striatum and variables not renamed yet)

plot_day_bins = [-Inf,-2:2,Inf];
wf_plot_days_grp = discretize(max(wf_ld_grp,-inf),plot_day_bins);

% Split by trial percentile within day
n_split = 3;
rxn_bins = [0.3,Inf];

rxn_bin_grp = discretize(wf_rxn_grp,rxn_bins);
[~,~,rec_idx] = unique([wf_animal_grp,wf_plot_days_grp],'rows');

split_idx = nan(size(rec_idx));
for curr_rec = 1:max(rec_idx)
    split_idx(rec_idx == curr_rec) = ...
        min(floor(linspace(1,n_split+1,sum(rec_idx == curr_rec))),n_split)';
end


% % Plot PSTH split in rxn bin 2
% figure('color','w'); h = tiledlayout(n_k,max(wf_plot_days_grp),'TileSpacing','none');
% for curr_k = 1:n_k
%     for curr_day = unique(wf_plot_days_grp)'
%         nexttile;
%         curr_trials = find(wf_plot_days_grp == curr_day & rxn_bin_grp == 2);
%         plot(wf_t,ap.nestgroupfun({@mean,@mean},wf_striatum_roi(curr_trials,:,curr_k), ...
%             wf_animal_grp(curr_trials),split_idx(curr_trials))');
%     end
% end
% linkaxes(h.Children,'xy');
% set(h.Children,'colororder',ap.colormap('KR',n_split)); 
% axis(h.Children,'off')



% Plot PSTH max split
figure('color','w'); h = tiledlayout(n_k,max(wf_plot_days_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day = unique(wf_plot_days_grp)'
        nexttile; hold on;
        curr_trials = find(wf_plot_days_grp == curr_day);
        curr_data = wf_striatum_roi(curr_trials,:,curr_k);

        data_t = wf_t > 0 & wf_t < 0.2;

        % (mean PSTH, max t, mean across animals)
        [psth_mean,psth_mean_grp] = ap.groupfun(@mean,curr_data, ...
            [wf_animal_grp(curr_trials),split_idx(curr_trials),rxn_bin_grp(curr_trials)]);
        psth_max = max(psth_mean(:,data_t),[],2);
        [psth_max_avg,psth_max_grp] = ap.nestgroupfun({@mean,@mean},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2:3));
        psth_max_sem = ap.nestgroupfun({@mean,@AP_sem},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2:3));

        arrayfun(@(x) errorbar(psth_max_avg(psth_max_grp(:,2) == x), ...
            psth_max_sem(psth_max_grp(:,2) == x),'linewidth',2), ...
            1:max(rxn_bin_grp));

        axis padded off
    end
end
linkaxes(h.Children,'xy');


%% Widefield (task, raw vs move-removed)

% Create group indicies
wf_t = wf.wf_stim_time{1};

wf_animal_grp = cell2mat(cellfun(@(animal,grp) repmat(animal,size(grp,1),1),    ...
    num2cell(grp2idx(wf.animal)),wf.V_stim_align,'uni',false));

wf_ld_grp = cell2mat(cellfun(@(ld,grp) repmat(ld,size(grp,1),1),    ...
    num2cell(bhv.days_from_learning),wf.V_stim_align,'uni',false));

plot_day_bins = [-Inf,-2:2,Inf];
wf_plot_days_grp = discretize(max(wf_ld_grp,-inf),plot_day_bins);

wf_rxn_grp = cell2mat(bhv.stim_to_move);

% Make ROIs from striatum maps
% % (weighted average)
% striatum_wf_roi = max(kmeans_cluster_mean,0)./max(kmeans_cluster_mean,[],[1,2]);
% (thresholded spot)
[~,m] = max(kmeans_cluster_mean,[],[1,2],'linear');
[mr,mc] = ind2sub(size(U_master,[1,2]),m);
striatum_wf_roi = kmeans_cluster_mean > max(kmeans_cluster_mean,[],[1,2])*0.5;
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
    axis image off; ap.wf_draw('ccf','k');
    colormap(ap.colormap('WG'));
end

% Get ROI activity
wf_striatum_roi = permute(ap.wf_roi(U_master, ...
    permute(cell2mat(wf.V_stim_align),[3,2,1]),[],[],striatum_wf_roi),[3,2,1]);

% (baseline-subtract)
wf_striatum_roi = wf_striatum_roi - nanmean(wf_striatum_roi(:,wf_t < -0.3,:),2);

% Plot average activity in trial (NaN-out movement times)
figure; h = tiledlayout(n_k,max(wf_plot_days_grp),'TileSpacing','none');
move_leeway = 0.1;
for curr_k = 1:n_k
    for curr_day = unique(wf_plot_days_grp)'
        nexttile; hold on;
        curr_trials = find(wf_plot_days_grp == curr_day);      
        
        % Average all data
        curr_data = wf_striatum_roi(curr_trials,:,curr_k);

        curr_data_mean = mean(ap.groupfun(@nanmean,curr_data,wf_animal_grp(curr_trials)),1);
        curr_data_sem = AP_sem(ap.groupfun(@nanmean,curr_data,wf_animal_grp(curr_trials)),1);

        % Average data with movement removede
        % (only plot no-move timepoints with at least N trials and N animals)
        curr_data_no_move = wf_striatum_roi(curr_trials,:,curr_k).* ...
            AP_nanout(wf_t > wf_rxn_grp(curr_trials)-move_leeway);

        min_nomove_trials = 10;
        min_noanimals_trials = 5;

        curr_data_no_move_n = ap.groupfun(@(x) sum(~isnan(x)),curr_data_no_move,wf_animal_grp(curr_trials));
        curr_data_no_move_animal_mean = ap.groupfun(@nanmean,curr_data_no_move,wf_animal_grp(curr_trials));
        curr_data_no_move_animal_mean_mindata = curr_data_no_move_animal_mean.* ...
            AP_nanout(curr_data_no_move_n < min_nomove_trials).* ...
            AP_nanout(sum(curr_data_no_move_n >= min_nomove_trials) < min_noanimals_trials);
        curr_data_no_move_mean = nanmean(curr_data_no_move_animal_mean_mindata,1);
        curr_data_no_move_sem = AP_sem(curr_data_no_move_animal_mean_mindata,1);

        % Plot
        ap.errorfill(wf_t,curr_data_mean,curr_data_sem,'k',[],[],2);
        ap.errorfill(wf_t,curr_data_no_move_mean,curr_data_no_move_sem,'r',[],[],2);
        axis off;
    end
end
linkaxes(h.Children,'xy');

plot_movie = false;
if plot_movie
    % Plot movie of non-move trials
    V_no_move = nan(size(U_master,3),length(wf_t),max(wf_plot_days_grp),'single');
    for curr_day = unique(wf_plot_days_grp)'

        curr_trials = wf_plot_days_grp == curr_day;
        curr_trials_rec = mat2cell(curr_trials,cellfun(@(x) size(x,1),wf.V_stim_align),1);

        curr_data_no_move = cell2mat(cellfun(@(V,trials) V(trials,:,:),wf.V_stim_align,curr_trials_rec,'uni',false)).* ...
            AP_nanout(wf_t > wf_rxn_grp(curr_trials)-move_leeway);

        min_nomove_trials = 10;
        min_noanimals_trials = 5;

        curr_data_no_move_n = ap.groupfun(@(x) sum(~isnan(x)),curr_data_no_move(:,:,1),wf_animal_grp(curr_trials));
        curr_data_no_move_animal_mean = ap.groupfun(@nanmean,curr_data_no_move,wf_animal_grp(curr_trials));
        curr_data_no_move_animal_mean_mindata = curr_data_no_move_animal_mean.* ...
            AP_nanout(curr_data_no_move_n < min_nomove_trials).* ...
            AP_nanout(sum(curr_data_no_move_n >= min_nomove_trials) < min_noanimals_trials);
        curr_data_no_move_mean = nanmean(curr_data_no_move_animal_mean_mindata,1);

        V_no_move(:,:,curr_day) = permute(curr_data_no_move_mean,[3,2,1]);
    end

    V_no_move_baselinesub = V_no_move - nanmean(V_no_move(:,wf_t < -0.3,:),2);
    ap.imscroll(plab.wf.svd2px(U_master,V_no_move_baselinesub));
    colormap(ap.colormap('PWG',[],1.5));
    clim(max(abs(clim)).*[-1,1]);
    axis image;
    ap.wf_draw('ccf','k');
end





