% Load data for longitudinal striatum+cortex project
%
% Required variable: 
% load_dataset - dataset to load ('task' or 'passive');
% 
% Optional: 
% load_dataset_overwrite - overwrite if data already loaded (true/false)

%% Check options

% Check for valid dataset selection
if ~exist('load_dataset','var')
    error('No load_dataset specified')
elseif ~ismember(load_dataset,{'passive','task'})
    error('Unknown load_dataset: %s',load_dataset);
end

% If no overwrite flag set, turn off by default
if ~exist('load_dataset_overwrite','var')
    load_dataset_overwrite = false;
end

% If selected dataset loaded and no overwrite, skip
if ~load_dataset_overwrite && ...
        exist('loaded_dataset','var') && strcmp(loaded_dataset,load_dataset)
    return
end

% Before loading, clear workspace
clearvars -except load_dataset
fprintf('Loading dataset: %s...\n',load_dataset);


%% Behavior

% Set stat and p-value to define learning day
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


%% K-means on maps 

% Load maps
data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');
load(fullfile(data_path,'ctx_maps_to_str'));
U_size = size(all_ctx_maps_to_str.cortex_kernel_px{1},[1,2]);

% K-means cluster maps
n_k = 3;

% (kmeans starting as maps averaged by depth)
kmeans_starting = nanmean(cell2mat(permute(cellfun(@(x) ap.groupfun(@mean,x,[],[],ap.quantile_bin(size(x,3),n_k)), ...
    all_ctx_maps_to_str.cortex_kernel_px(cellfun(@(x) ...
    size(x,3) >= n_k,all_ctx_maps_to_str.cortex_kernel_px)),'uni',false),[2,3,4,1])),4);

[kidx,kmeans_centroid_flat,~,kmeans_dist] = kmeans(...
    reshape(cat(3,all_ctx_maps_to_str.cortex_kernel_px{:}),prod(U_size),[])',n_k, ...
    'Distance','correlation','start',reshape(kmeans_starting,[],n_k)');
kmeans_centroid = reshape(kmeans_centroid_flat',size(kmeans_starting,1),size(kmeans_starting,2),[]);

% Package kidx by recording (used for some groupings)
kidx_rec = mat2cell(kidx,cellfun(@(x) size(x,3).*(size(x,1)>0), ...
    all_ctx_maps_to_str.cortex_kernel_px));


%% Widefield and ephys

% Set data path
data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');

% Load master U
U_master = plab.wf.load_master_U;

switch load_dataset

    case 'passive'

        % Load ephys
        load(fullfile(data_path,'ephys'));

        % Load widefield
        load(fullfile(data_path,'ctx_wf'));
        % (TO DO: change saved field to wf.V_no_move_stim_align>V_stim_align
        wf = renamevars(wf,'V_no_move_stim_align','V_stim_align');

    case 'task'

        % Load ephys
        load(fullfile(data_path,'task_ephys_outcome'));
        % (TO DO: change this saved variable to be 'ephys', not 'task_ephys'
        ephys = task_ephys;
        clear task_ephys;

        % Load widefield
        load(fullfile(data_path,'task_ctx_wf'));

end


%% Ephys: create cluster MUA and data indices

% Sum spikes from depth into cluster multiunit
[striatum_mua_sum,striatum_kidx] = cellfun(@(mua,kidx) ...
    ap.groupfun(@sum,mua,[],[],kidx), ...
    ephys.binned_msn_spikes_stim_align,kidx_rec,'uni',false);

% Normalize and smooth multiunit
% (currently psth time is hard coded: save this somewhere)
psth_t = -0.5:0.001:1;
baseline_t = psth_t < 0;
softnorm = 10;
mua_baseline = cellfun(@(mua) ...
    mean(mua(:,baseline_t,:),[1,2]) + softnorm, ...
    striatum_mua_sum,'uni',false,'ErrorHandler',@(varargin) NaN);

spikes_norm_smooth_reshape_fcn = @(spikes,baseline) ...
    reshape(permute(smoothdata((spikes-baseline)./baseline,2, ...
    'gaussian',[100,0]),[1,3,2]),[],length(psth_t));

striatum_mua = cell2mat(cellfun(spikes_norm_smooth_reshape_fcn, ...
    striatum_mua_sum,mua_baseline,'uni',false));

% Create grouping indicies for striatum MUA
striatum_mua_grp = struct;
striatum_trials_rec_n = cellfun(@(x) size(x,1),striatum_mua_sum);
striatum_mua_rec_n = cellfun(@length,striatum_kidx);

striatum_mua_grp.kidx = cell2mat(cellfun(@(grp,n_trials) reshape(repmat(grp',n_trials,1),[],1), ...
    striatum_kidx,num2cell(striatum_trials_rec_n),'uni',false));

striatum_mua_grp.animal = cell2mat(cellfun(@(grp,n_trials,n_mua) repmat(grp,n_trials*n_mua,1), ...
    num2cell(grp2idx(ephys.animal)),num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

striatum_mua_grp.ld = cell2mat(cellfun(@(grp,n_trials,n_mua) repmat(grp,n_trials*n_mua,1), ...
    num2cell(bhv.days_from_learning),num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

switch load_dataset
    case 'task'
        % Reaction times index
        striatum_mua_grp.rxn = cell2mat(cellfun(@(rxn,n_mua) reshape(repmat(rxn,1,n_mua),[],1), ...
            bhv.stim_to_move,num2cell(striatum_mua_rec_n),'uni',false));

        % Create move-aligned MUA
        striatum_move_psth_sum = cellfun(@(mua,mua_group) ap.groupfun(@sum,mua,[],[],mua_group), ...
            ephys.binned_msn_spikes_move_align,kidx_rec,'uni',false);
        striatum_move_psth = cell2mat(cellfun(spikes_norm_smooth_reshape_fcn, ...
            striatum_move_psth_sum,mua_baseline,'uni',false));

        % Create move-aligned MUA
        striatum_outcome_psth_sum = cellfun(@(mua,mua_group) ap.groupfun(@sum,mua,[],[],mua_group), ...
            ephys.binned_msn_spikes_outcome_align,kidx_rec,'uni',false);
        striatum_outcome_psth = cell2mat(cellfun(spikes_norm_smooth_reshape_fcn, ...
            striatum_outcome_psth_sum,mua_baseline,'uni',false));

    case 'passive'
        % Stim index
        striatum_mua_grp.stim = cell2mat(cellfun(@(stim,n_mua) reshape(repmat(stim,1,n_mua),[],1), ...
            ephys.trial_stim_values,num2cell(striatum_mua_rec_n),'uni',false));
end


%% Ephys: create indicies for single units

% Normalize and smooth single
psth_t = -0.5:0.001:1;
baseline_t = psth_t < 0;
softnorm = 1;
sua_baseline = cellfun(@(sua) ...
    mean(sua(:,baseline_t,:),[2,3]) + softnorm, ...
    ephys.unit_event_psths,'uni',false,'ErrorHandler',@(varargin) NaN);

spikes_norm_smooth_reshape_fcn = @(spikes,baseline) ...
    smoothdata((spikes-baseline)./baseline,2, ...
    'gaussian',[100,0]);

% (units outside the striatum have ephys.unit_depth_group = NaN)
striatum_units = cellfun(@(x) ~isnan(x),ephys.unit_depth_group,'uni',false);

striatum_sua = cell2mat(cellfun(@(data,baseline,striatum_units) ...
    spikes_norm_smooth_reshape_fcn(data(striatum_units,:,:),baseline(striatum_units)), ...
    ephys.unit_event_psths,sua_baseline,striatum_units,'uni',false));

% Create grouping indicies for striatum SUA
striatum_sua_grp = struct;
striatum_units_n_rec = cellfun(@sum,striatum_units);

striatum_sua_grp.animal = cell2mat(cellfun(@(grp,n_units) repmat(grp,n_units,1), ...
    num2cell(grp2idx(bhv.animal)),num2cell(striatum_units_n_rec),'uni',false));

striatum_sua_grp.ld = cell2mat(cellfun(@(grp,n_units) repmat(grp,n_units,1), ...
    num2cell(bhv.days_from_learning),num2cell(striatum_units_n_rec),'uni',false));

striatum_sua_grp.kidx = cell2mat(cellfun(@(depth,kidx,striatum_units) kidx(depth(striatum_units)), ...
    ephys.unit_depth_group,kidx_rec,striatum_units,'uni',false));



%% Widefield: create indices

% Grab aligned data time
wf_t = wf.wf_stim_time{1};

% Create group indicies
wf_grp = struct;

wf_grp.animal = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
    num2cell(grp2idx(wf.animal)),wf.V_stim_align,'uni',false));

wf_grp.ld = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
    num2cell(bhv.days_from_learning),wf.V_stim_align,'uni',false));

switch load_dataset
    case 'task'
        % Reaction time index
        wf_grp.rxn = cell2mat(bhv.stim_to_move);
    case 'passive'
        % Stim index
        wf_grp.stim = cell2mat(wf.trial_stim_values);
end

cortex_trials_rec_n = cellfun(@(x) size(x,1),wf.V_stim_align);


%% Widefield ROIs by corticostriatal maps

% Create ROIs by striatum cluster maps
kmeans_centroid_blur = imgaussfilt(kmeans_centroid,10);
[~,m] = max(kmeans_centroid_blur,[],[1,2],'linear');
[mr,mc] = ind2sub(size(U_master,[1,2]),m);
striatum_wf_roi = kmeans_centroid_blur > prctile(kmeans_centroid_blur,100,[1,2])*0.50;
striatum_wf_roi(:,round(size(striatum_wf_roi,2)/2):end,:) = false;
for k = 1:size(striatum_wf_roi,3)
    [~,m] = max(imgaussfilt(kmeans_centroid_blur(:,:,k),10).* ...
        round(linspace(1,0,size(kmeans_centroid_blur,2))),[],[1,2],'linear');
    [mx,my] = ind2sub(size(U_master,[1,2]),m);
    striatum_wf_roi(:,:,k) = bwselect(striatum_wf_roi(:,:,k),my,mx);
end

% Get ROI activity
wf_striatum_roi = permute(ap.wf_roi(U_master, ...
    permute(cell2mat(wf.V_stim_align),[3,2,1]),[],[],striatum_wf_roi),[3,2,1]);

% (baseline-subtract)
baseline_t = wf_t < 0;
wf_striatum_roi = wf_striatum_roi - nanmean(wf_striatum_roi(:,baseline_t,:),2);


%% Set flag for loaded data

loaded_dataset = load_dataset;
fprintf('Finished loading dataset: %s\n',loaded_dataset);

