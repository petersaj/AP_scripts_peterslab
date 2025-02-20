
%% Load AM-packaged data

data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');

load(fullfile(data_path,'swr_bhv'));
load(fullfile(data_path,'ctx_maps_to_str'));
% load(fullfile(data_path,'ctx_wf'));
load(fullfile(data_path,'ephys'));
U_master = plab.wf.load_master_U;

% K-means cluster maps
% (experiment 2 is bad - AM will re-run) 
maps_cat = cat(3,all_ctx_maps_to_str.cortex_kernel_px{:});
expl_var_cat = vertcat(all_ctx_maps_to_str.explained_var{:});

n_k = 4;
kidx = nan(size(maps_cat,3),1);
kidx_use_maps = expl_var_cat > 0;
[kidx(kidx_use_maps),kmeans_map] = kmeans(...
    reshape(maps_cat(:,:,kidx_use_maps),prod(size(U_master,[1,2])),[])',n_k, ...
    'Distance','Correlation','Replicates',5);
kmeans_map = reshape(kmeans_map',size(U_master,1),size(U_master,2),[]);

figure;
imagesc(reshape(permute(kmeans_map,[1,3,2]),[],size(kmeans_map,2)))
axis image off
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));

kidx_rec = mat2cell(kidx,cellfun(@(x) size(x,3).*(size(x,1)>0), ...
    all_ctx_maps_to_str.cortex_kernel_px));

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

[unit_group_mean,groups] = ap.nestgroupfun({@mean,@mean}, ...
    +unit_resp_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

unit_group_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    +unit_resp_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

figure; hold on;
h = arrayfun(@(x) ap.errorfill(groups(groups(:,2) == x,1),unit_group_mean(groups(:,2) == x,1),unit_group_sem(groups(:,2) == x,1)),unique(groups(:,2)));
legend(h,string(num2cell(unique(groups(:,2)))));


% (response amplitude)
curr_stim = 3;
use_units = unit_single_cat;

[unit_group_mean,groups] = ap.nestgroupfun({@mean,@mean}, ...
    unit_mean_post_stim_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

unit_group_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    unit_mean_post_stim_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

figure; hold on;
h = arrayfun(@(x) ap.errorfill(groups(groups(:,2) == x,1),unit_group_mean(groups(:,2) == x,1),unit_group_sem(groups(:,2) == x,1)),unique(groups(:,2)));
legend(h,string(num2cell(unique(groups(:,2)))));


% (psth)
curr_stim = 3;
% use_units = true(size(unit_kidx));
% use_units = unit_single_cat;
use_units = unit_resp_cat(:,3);

tic
[unit_group_mean,groups] = ap.nestgroupfun({@mean,@mean}, ...
    unit_psth_cat{curr_stim}(use_units,:),group_labels(use_units,:),split_labels(use_units,:));
toc

unit_group_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    unit_psth_cat{curr_stim}(use_units,:),group_labels(use_units,:),split_labels(use_units,:));

plot_days = -3:2;

day_colormap_sym = ap.colormap('BKR',max(abs(plot_days))*2+1);
day_colormap = day_colormap_sym(ismember( ...
    -max(abs(plot_days)):max(abs(plot_days)),plot_days),:);

figure; 
h = tiledlayout(n_k,1);
for curr_k = 1:n_k
    nexttile; set(gca,'ColorOrder',day_colormap);
    arrayfun(@(x) ...
        ap.errorfill([], ...
        unit_group_mean(groups(:,1) == x & groups(:,2) == curr_k,:)', ...
        unit_group_sem(groups(:,1) == x & groups(:,2) == curr_k,:)'), ...
        plot_days);
end
linkaxes(h.Children,'xy');




[unit_group_mean,groups] = ap.nestgroupfun({@mean,@mean}, ...
    +unit_resp_cat(use_units,curr_stim),group_labels(use_units,:),split_labels(use_units,:));

















