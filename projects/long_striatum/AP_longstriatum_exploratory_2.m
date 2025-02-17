
%% Load AM-packaged data

data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');

load(fullfile(data_path,'swr_bhv'));
load(fullfile(data_path,'ctx_maps_to_str'));
% load(fullfile(data_path,'ctx_wf'));
load(fullfile(data_path,'ephys'));
U_master = plab.wf.load_master_U;

% K-means cluster maps
% (why is expt 2 crazy?)
maps_cat = cat(3,all_ctx_maps_to_str.cortex_kernel_px{:});
expl_var_cat = vertcat(all_ctx_maps_to_str.explained_var{:});

n_k = 4;
kidx = nan(size(maps_cat,3),1);
kidx_use_maps = expl_var_cat > 0;
[kidx(kidx_use_maps),kmeans_map] = kmeans(...
    reshape(maps_cat(:,:,kidx_use_maps),prod(size(U_master,[1,2])),[])',n_k);
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

% Make index of kidx per unit
unit_kidx_subset = cell2mat(cellfun(@(depth,kidx) kidx(depth(~isnan(depth))), ...
    ephys.unit_depth_group,kidx_rec,'uni',false));
unit_kidx = nan(size(unit_rec_idx));
unit_kidx(~isnan(cell2mat(ephys.unit_depth_group))) = unit_kidx_subset;

% Concatenate unit data (responsive, single, psth, mean)
unit_resp_cat = cell2mat(cellfun(@(x) cell2mat(x)', ...
    horzcat(ephys.unit_resp_p_value{:})','uni',false));

unit_single_cat = cell2mat(cellfun(@logical,ephys.single_unit_idx,'uni',false));

unit_psth_cat = arrayfun(@(stim) cell2mat(cellfun(@(x) x{stim}, ...
    ephys.unit_smooth_event_psths,'ErrorHandler',@(varargin) [],'uni',false)), ...
    1:size(unit_resp_cat,2),'uni',false);

unit_mean_post_stim_cat = cell2mat(vertcat(ephys.mean_post_stim{:}));

% Make index of recording per unit

unit_animal = cell2mat(cellfun(@(grp,units) repmat(grp,length(units),1), ...
    num2cell(grp2idx(bhv.animal)),ephys.single_unit_idx,'uni',false));

unit_ld = cell2mat(cellfun(@(grp,units) repmat(grp,length(units),1), ...
    num2cell(bhv.days_from_learning),ephys.single_unit_idx,'uni',false));

unit_kidx_subset = cell2mat(cellfun(@(depth,kidx) kidx(depth(~isnan(depth))), ...
    ephys.unit_depth_group,kidx_rec,'uni',false));
unit_kidx = nan(size(unit_rec_idx));
unit_kidx(~isnan(cell2mat(ephys.unit_depth_group))) = unit_kidx_subset;

% Group data and plot
group = [unit_kidx,unit_ld,unit_animal];
group_use = ~any(isnan(group),2);

unit_group_idx = nan(size(group,1),1);
[group_unique,~,unit_group_idx(group_use)] = unique(use_grp(group_use,:),'rows');

unit_kidx_rec_grouped = ap.groupfun(@mean,unit_mean_post_stim_cat(:,curr_stim),unit_group_idx,[]);

figure;
for curr_k = 1:n_k
    nexttile
end
























