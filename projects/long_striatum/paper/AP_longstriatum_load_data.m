% Load data for longitudinal striatum+cortex project
%
% Required variable:
% load_dataset - dataset to load ('task','passive','noact');
%
% Optional:
% load_dataset_overwrite - overwrite if data already loaded (true/false)

%% Check options

% Check for valid dataset selection
if ~exist('load_dataset','var')
    error('No load_dataset specified')
elseif ~ismember(load_dataset,{'passive','task','noact'})
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

% Before loading, clear workspace except for fig settings
clearvars -except load_dataset fig_save_flag stat_fid print_stat save_figs
fprintf('Loading dataset: %s...\n',load_dataset);

% For printing stats: if no stat file, print to command line
if ~exist('stat_fid','var') || ~stat_fid
    print_stat = @(varargin) fprintf(varargin{:});
end


%% Set data path

data_path = fullfile(plab.locations.server_path,'Lab','Papers','Marica_2025','data');


%% Behavior

% Set stat and p-value to define learning day
use_stat = 'firstmove_mean';
learn_stat_p = 0.05;

% Load behavior
load(fullfile(data_path,'bhv'));

% Set "learned_days" and "days_from_learning"
bhv.learned_days = cellfun(@(x) x < learn_stat_p,bhv.(['stimwheel_pval_',use_stat]));
for curr_animal = unique(bhv.animal)'
    curr_rows = strcmp(bhv.animal,curr_animal);
    bhv.days_from_learning(curr_rows) = ...
        (1:sum(curr_rows))' - ...
        max([NaN,find(bhv.learned_days(curr_rows),1)]);
end


%% K-means on maps

% Load maps
load(fullfile(data_path,'ctx_str_maps'));
% load(fullfile(data_path,'ctx_str_maps_Umaster'));
U_size = size(ctx_str_maps.cortex_striatum_map{1},[1,2]);

% K-means cluster maps
n_domains = 3;

% (kmeans starting as maps averaged by depth)
kmeans_starting = nanmean(cell2mat(permute(cellfun(@(x) ...
    ap.groupfun(@mean,x,[],[],ap.quantile_bin(size(x,3),n_domains)), ...
    ctx_str_maps.cortex_striatum_map(cellfun(@(x) ...
    size(x,3) >= n_domains,ctx_str_maps.cortex_striatum_map)),'uni',false),[2,3,4,1])),4);

[domain_idx,kmeans_centroid_flat] = kmeans(...
    reshape(cat(3,ctx_str_maps.cortex_striatum_map{:}),prod(U_size),[])',n_domains, ...
    'Distance','correlation','start',reshape(kmeans_starting,[],n_domains)');

kmeans_centroid = reshape(kmeans_centroid_flat', ...
    size(kmeans_starting,1),size(kmeans_starting,2),[]);

kmeans_cluster_mean = reshape(ap.groupfun(@nanmean, ...
    reshape(cat(3,ctx_str_maps.cortex_striatum_map{:}), ...
    prod(U_size),[]),[],domain_idx),[U_size,n_domains]);

% Package domain_idx by recording
domain_idx_rec = mat2cell(domain_idx,cellfun(@(x) size(x,3).*(size(x,1)>0), ...
    ctx_str_maps.cortex_striatum_map));

% %%%% TESTING
% % cumulative max
% domain_idx_rec = cellfun(@cummax,domain_idx_rec,'uni',false);
% domain_idx = cell2mat(domain_idx_rec);
% 
% % moving median
% domain_idx_rec = cellfun(@(x) medfilt1(x,3),domain_idx_rec,'uni',false);
% domain_idx = cell2mat(domain_idx_rec);
% %%%%%%


%% Widefield and ephys

% Load master U
U_master = plab.wf.load_master_U;

switch load_dataset
    case 'passive'
        load(fullfile(data_path,'ephys_passive'));
        load(fullfile(data_path,'wf_passive'));
    case 'task'
        load(fullfile(data_path,'ephys_task'));
        load(fullfile(data_path,'wf_task'));
end


%% Ephys: create cluster MUA and data indices

if ~strcmp(load_dataset,'noact')

    % Sum spikes from depth into cluster multiunit
    [striatum_mua_sum,striatum_domain_idx] = cellfun(@(mua,domain_idx) ...
        ap.groupfun(@sum,mua,[],[],domain_idx), ...
        ephys.binned_msn_spikes_event_align,domain_idx_rec,'uni',false);
    % (standardize dimensionality of empty data)
    striatum_mua_sum(cellfun(@isempty,striatum_mua_sum)) = ...
        {zeros([0,size(striatum_mua_sum{find(~cellfun(@isempty,striatum_mua_sum),1)},[2,3,4])])};

    % Normalize and smooth multiunit
    % (currently psth time is hard coded: save this somewhere)
    psth_t = -0.5:0.001:1;
    baseline_t = psth_t < 0;
    softnorm = 10;

    % (to use day baseline)
    % mua_baseline = cellfun(@(mua) ...
    %     repmat(mean(mua(:,baseline_t,:,1),[1,2]),1,1,1,size(mua,4)), ...
    %     striatum_mua_sum,'uni',false,'ErrorHandler',@(varargin) NaN);
    % OR
    % (to use trial baseline)
    mua_baseline = cellfun(@(mua) ...
        repmat(mean(mua(:,baseline_t,:,1),[2]),1,1,1,size(mua,4)), ...
        striatum_mua_sum,'uni',false,'ErrorHandler',@(varargin) NaN);

    spikes_norm_smooth_reshape_fcn = @(spikes,baseline) ...
        reshape(permute(smoothdata((spikes-baseline)./(baseline+softnorm),2, ...
        'gaussian',[100,0]),[1,3,2,4]),[],length(psth_t),size(spikes,4));

    striatum_mua = cell2mat(cellfun(spikes_norm_smooth_reshape_fcn, ...
        striatum_mua_sum,mua_baseline,'uni',false));

    % Create grouping indicies for striatum MUA
    striatum_mua_grp = struct;
    striatum_trials_rec_n = cellfun(@(x) size(x,1),striatum_mua_sum);
    striatum_mua_rec_n = cellfun(@length,striatum_domain_idx);

    striatum_mua_grp.domain_idx = cell2mat(cellfun(@(grp,n_trials) reshape(repmat(grp',n_trials,1),[],1), ...
        striatum_domain_idx,num2cell(striatum_trials_rec_n),'uni',false));

    striatum_mua_grp.animal = cell2mat(cellfun(@(grp,n_trials,n_mua) repmat(grp,n_trials*n_mua,1), ...
        num2cell(grp2idx(ephys.animal)),num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

    striatum_mua_grp.ld = cell2mat(cellfun(@(grp,n_trials,n_mua) repmat(grp,n_trials*n_mua,1), ...
        num2cell(bhv.days_from_learning),num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

    switch load_dataset
        case 'task'
            % Reaction times index
            striatum_mua_grp.rxn = cell2mat(cellfun(@(rxn,n_mua) reshape(repmat(rxn,1,n_mua),[],1), ...
                bhv.stim_to_move,num2cell(striatum_mua_rec_n),'uni',false));
        case 'passive'
            % Stim index
            striatum_mua_grp.stim = cell2mat(cellfun(@(stim,n_mua) reshape(repmat(stim,1,n_mua),[],1), ...
                ephys.trial_stim_values,num2cell(striatum_mua_rec_n),'uni',false));
    end

end


%% Ephys: concatenate single units and create indicies

if ~strcmp(load_dataset,'noact')

    % Normalize and smooth single unit PSTHs
    psth_t = -0.5:0.001:1;
    baseline_t = psth_t < 0;
    softnorm = 1;
    switch load_dataset
        case 'task'
            % (task baseline is only align 1 = stim-align)
            sua_baseline = cellfun(@(sua) ...
                mean(sua(:,baseline_t,1),[2,3]), ...
                ephys.unit_event_psths,'uni',false,'ErrorHandler',@(varargin) NaN);
        case 'passive'
            % (passive baseline is average of all stim align)
            sua_baseline = cellfun(@(sua) ...
                mean(sua(:,baseline_t,:),[2,3]), ...
                ephys.unit_event_psths,'uni',false,'ErrorHandler',@(varargin) NaN);
    end

    spikes_norm_smooth_reshape_fcn = @(spikes,baseline) ...
        smoothdata((spikes-baseline)./(baseline+softnorm),2, ...
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

    striatum_sua_grp.domain_idx = cell2mat(cellfun(@(depth,domain_idx,striatum_units) domain_idx(depth(striatum_units)), ...
        ephys.unit_depth_group,domain_idx_rec,striatum_units,'uni',false));

    % (indicies to find specific cell in recording)
    striatum_sua_grp.rec = cell2mat(cellfun(@(grp,n_units) repmat(grp,n_units,1), ...
        num2cell(1:length(striatum_units_n_rec))',num2cell(striatum_units_n_rec),'uni',false));

    striatum_sua_grp.unit_id = cell2mat(cellfun(@(grp) find(grp), ...
        striatum_units,'uni',false));

    % (bombcell cell types)
    striatum_sua_grp.msn = cell2mat(cellfun(@(grp,str) logical(grp(str)), ...
        ephys.str_msn_idx,striatum_units,'uni',false));

    striatum_sua_grp.fsi = cell2mat(cellfun(@(grp,str) logical(grp(str)), ...
        ephys.str_fsi_idx,striatum_units,'uni',false));

    striatum_sua_grp.tan = cell2mat(cellfun(@(grp,str) logical(grp(str)), ...
        ephys.str_tan_idx,striatum_units,'uni',false));

    striatum_sua_grp.good_unit = cell2mat(cellfun(@(grp,str) logical(grp(str)), ...
        ephys.single_unit_idx,striatum_units,'uni',false));

end


%% Widefield: create indices

if ~strcmp(load_dataset,'noact')

    % Grab aligned data time
    wf_t = wf.wf_stim_time{1};

    % Create group indicies
    wf_grp = struct;

    wf_grp.animal = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
        num2cell(grp2idx(wf.animal)),wf.V_event_align,'uni',false));

    wf_grp.ld = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
        num2cell(bhv.days_from_learning),wf.V_event_align,'uni',false));

    switch load_dataset
        case 'task'
            % Reaction time index
            wf_grp.rxn = cell2mat(bhv.stim_to_move);
        case 'passive'
            % Stim index
            wf_grp.stim = cell2mat(wf.trial_stim_values);
    end

    cortex_trials_rec_n = cellfun(@(x) size(x,1),wf.V_event_align);

end


%% Widefield ROIs by corticostriatal maps

% Create ROIs by striatum cluster maps
kmeans_centroid_blur = imgaussfilt(kmeans_cluster_mean,10);
striatum_wf_roi = kmeans_centroid_blur > prctile(kmeans_centroid_blur,100,[1,2])*0.75;

% (not needed anymore after increasing thresold?)
% striatum_wf_roi(:,round(size(striatum_wf_roi,2)/2):end,:) = false;
% for k = 1:size(striatum_wf_roi,3)
%     [~,m] = max(imgaussfilt(kmeans_centroid_blur(:,:,k),10).* ...
%         round(linspace(1,0,size(kmeans_centroid_blur,2))),[],[1,2],'linear');
%     [mx,my] = ind2sub(size(U_master,[1,2]),m);
%     striatum_wf_roi(:,:,k) = bwselect(striatum_wf_roi(:,:,k),my,mx);
% end

if ~strcmp(load_dataset,'noact')

    % Get ROI activity

    wf_striatum_roi = cell2mat(cellfun(@(x) ...
        permute(ap.wf_roi(U_master,permute(x,[3,2,1,4]),[],[],striatum_wf_roi), ...
        [3,2,1,4]),wf.V_event_align,'uni',false));

    % (baseline-subtract)
    baseline_t = wf_t < 0;
    wf_striatum_roi = wf_striatum_roi - nanmean(wf_striatum_roi(:,baseline_t,:,:),2);

end


%% Set flag for loaded data

loaded_dataset = load_dataset;
fprintf('Finished loading dataset: %s\n',loaded_dataset);

