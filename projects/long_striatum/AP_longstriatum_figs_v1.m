%% Notes

% Figures for longitudinal striatum+cortex project

% Data packaged by AM save scripts
% Analysis taken fom AP_longstriautm_exploratory_2

%% ------- LOAD/PREP DATA ------------------

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


%% Passive widefield + ephys

dataset = 'passive';

% Set data path
data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');

% Load ephys
load(fullfile(data_path,'ephys'));

% Load widefield
load(fullfile(data_path,'ctx_wf'));
% (TO DO: change saved field to wf.V_no_move_stim_align>V_stim_align
wf = renamevars(wf,'V_no_move_stim_align','V_stim_align');

% Load master U
U_master = plab.wf.load_master_U;


%% Task widefield + ephys

dataset = 'task';

% Set data path
data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');

% Load ephys
load(fullfile(data_path,'task_ephys_outcome'));
% (TO DO: change this saved variable to be 'ephys', not 'task_ephys'
ephys = task_ephys;
clear task_ephys;

% Load widefield
load(fullfile(data_path,'task_ctx_wf'));

% Load master U
U_master = plab.wf.load_master_U;


%% Striatum: create cluster MUA and data indices

% Sum spikes from depth into cluster multiunit
[striatum_psth_sum,striatum_kidx] = cellfun(@(mua,kidx) ...
    ap.groupfun(@sum,mua,[],[],kidx), ...
    ephys.binned_msn_spikes_stim_align,kidx_rec,'uni',false);

% Normalize and smooth multiunit
% (currently psth time is hard coded: save this somewhere)
psth_t = -0.5:0.001:1;
baseline_t = psth_t < 0;
softnorm = 10;
psth_baseline = cellfun(@(mua) ...
    mean(mua(:,baseline_t,:),[1,2]) + softnorm, ...
    striatum_psth_sum,'uni',false,'ErrorHandler',@(varargin) NaN);

psth_norm_smooth_reshape_fcn = @(psth,baseline) ...
    reshape(permute(smoothdata((psth-baseline)./baseline,2, ...
    'gaussian',[100,0]),[1,3,2]),[],length(psth_t));

striatum_psth = cell2mat(cellfun(psth_norm_smooth_reshape_fcn, ...
    striatum_psth_sum,psth_baseline,'uni',false));

% Create grouping indicies for all striatum data
striatum_psth_grp = struct;

n_trials_rec = cellfun(@(x) size(x,1),striatum_psth_sum);
n_mua_rec = cellfun(@length,striatum_kidx);

striatum_psth_grp.kidx = cell2mat(cellfun(@(kidx,n_trials) reshape(repmat(kidx',n_trials,1),[],1), ...
    striatum_kidx,num2cell(n_trials_rec),'uni',false));

striatum_psth_grp.animal = cell2mat(cellfun(@(animal,n_trials,n_mua) repmat(animal,n_trials*n_mua,1), ...
    num2cell(grp2idx(ephys.animal)),num2cell(n_trials_rec),num2cell(n_mua_rec),'uni',false));

striatum_psth_grp.ld = cell2mat(cellfun(@(animal,n_trials,n_mua) repmat(animal,n_trials*n_mua,1), ...
    num2cell(bhv.days_from_learning),num2cell(n_trials_rec),num2cell(n_mua_rec),'uni',false));

switch dataset
    case 'task'
        striatum_psth_grp.rxn = cell2mat(cellfun(@(rxn,n_mua) reshape(repmat(rxn,1,n_mua),[],1), ...
            bhv.stim_to_move,num2cell(n_mua_rec),'uni',false));
    case 'passive'
        striatum_psth_grp.stim = cell2mat(cellfun(@(stim,n_mua) reshape(repmat(stim,1,n_mua),[],1), ...
            ephys.trial_stim_values,num2cell(n_mua_rec),'uni',false));
end



%% Widefield: create indices

% Grab aligned data time
wf_t = wf.wf_stim_time{1};

% Create group indicies
wf_grp = struct;

wf_grp.animal = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
    num2cell(grp2idx(wf.animal)),wf.V_stim_align,'uni',false));

wf_grp.ld = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
    num2cell(bhv.days_from_learning),wf.V_stim_align,'uni',false));

switch dataset
    case 'task'
        wf_grp.rxn = cell2mat(bhv.stim_to_move);
    case 'passive'
        wf_grp.stim = cell2mat(wf.trial_stim_values);
end



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


%% ------- GENERATE FIGURES ----------------


%% [Fig 1X] Reaction stat and learning histogram

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
            use_plot_day_grp = training_day;
            plot_x = 1:7;
        case 2
            use_plot_day_grp = bhv.days_from_learning;
            plot_x = -3:2;
    end

    [rxn_idx_mean,rxn_group_x] = ap.groupfun(@mean,rxn_idx_cat,use_plot_day_grp);
    rxn_idx_sem = ap.groupfun(@AP_sem,rxn_idx_cat,use_plot_day_grp);

    rxn_mean = ap.groupfun(@mean,rxn_cat,use_plot_day_grp);
    rxn_sem = ap.groupfun(@AP_sem,rxn_cat,use_plot_day_grp);

    rxn_null_mean = ap.groupfun(@mean,rxn_null_cat,use_plot_day_grp);
    rxn_null_sem = ap.groupfun(@AP_sem,rxn_null_cat,use_plot_day_grp);

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
n_learned_day = cellfun(@(x) max([0, ...
    find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'));
figure;
histogram(n_learned_day,-0.5:max(n_learned_day)+0.5)
ylabel('Number of mice');
xlabel('Days to learn');


%% [Fig 1X] Cortex map depths pre/post learning

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


%% [Fig 2X] K-means cluster centroids

figure;
h = tiledlayout(n_k,1,'tilespacing','none');
for curr_k = 1:n_k
    nexttile;
    imagesc(kmeans_centroid(:,:,curr_k));
    axis image off
    clim(max(abs(clim)).*[-1,1]);
    colormap(ap.colormap('PWG'));
    ap.wf_draw('ccf',[0.5,0.5,0.5]);
end

ap.prettyfig;

%% [Fig 2X] Striatum task trial heatmap (reaction-sorted)

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(striatum_psth_grp.ld,-inf),plot_day_bins);

heatmap_smooth = [20,1]; % ([trials,time] to smooth for graphics)
figure; tiledlayout(n_k,max(plot_day_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(plot_day_grp)'
        nexttile;
        curr_trials = find(plot_day_grp == curr_day & striatum_psth_grp.kidx == curr_k);
       
        [~,sort_idx] = sort(striatum_psth_grp.rxn(curr_trials));      
        imagesc(psth_t,[],movmean(striatum_psth(curr_trials(sort_idx),:),heatmap_smooth));
        colormap(ap.colormap('WK'));
        clim([0,3]);
        axis off;
        if curr_k == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end

ap.prettyfig;

%% [Fig 2X] Striatum task average PSTH w/ and w/o movement

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(striatum_psth_grp.ld,-inf),plot_day_bins);

% Plot average activity in trial (NaN-out movement times)
figure; h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','none');
move_leeway = 0.1; % time pre-movement to exclude
for curr_k = 1:n_k
    for curr_day = unique(plot_day_grp)'
        nexttile; hold on;
        curr_trials = find(plot_day_grp == curr_day & striatum_psth_grp.kidx == curr_k);      
        
        % Average all data
        curr_data = striatum_psth(curr_trials,:);
        curr_data_no_move = striatum_psth(curr_trials,:).* ...
            ap.nanout(psth_t > striatum_psth_grp.rxn(curr_trials)-move_leeway);

        curr_data_mean = mean(ap.groupfun(@nanmean,curr_data,striatum_psth_grp.animal(curr_trials)),1);
        curr_data_sem = AP_sem(ap.groupfun(@nanmean,curr_data,striatum_psth_grp.animal(curr_trials)),1);

        % Average data with movement removed
        % (only plot no-move timepoints with at least N trials and N animals)
        min_nomove_trials = 10;
        min_noanimals_trials = 5;

        curr_data_no_move_n = ap.groupfun(@(x) sum(~isnan(x)),curr_data_no_move,striatum_psth_grp.animal(curr_trials));
        curr_data_no_move_animal_mean = ap.groupfun(@nanmean,curr_data_no_move,striatum_psth_grp.animal(curr_trials));
        curr_data_no_move_animal_mean_mindata = curr_data_no_move_animal_mean.* ...
            ap.nanout(curr_data_no_move_n < min_nomove_trials).* ...
            ap.nanout(sum(curr_data_no_move_n >= min_nomove_trials) < min_noanimals_trials);
        curr_data_no_move_mean = nanmean(curr_data_no_move_animal_mean_mindata,1);
        curr_data_no_move_sem = AP_sem(curr_data_no_move_animal_mean_mindata,1);

        % Plot w/ movement (shaded) and w/o movement (line)
        ap.errorfill(psth_t,curr_data_mean,curr_data_sem,'k',0.3,false,2);
        plot(psth_t,curr_data_no_move_mean,'k','linewidth',2);
        axis off;
        if curr_k == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end
linkaxes(h.Children,'xy');

ap.prettyfig;


%% [Fig 2X] Striatum task max stim activity split aross session

rxn_cutoff = 0.3; % only plot trials with slow reaction times

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(ld_grp,-inf),plot_day_bins);

% Split by trial percentile within day
n_split = 3;
split_idx = cell2mat(cellfun(@(n_trials,n_mua) ...
    reshape(repmat(ap.quantile_bin(n_trials,n_split),1,n_mua),[],1), ...
    num2cell(n_trials_rec),num2cell(n_mua_rec),'uni',false));

% (mean PSTH, max t, mean across animals)
use_trials = striatum_psth_grp.rxn > rxn_cutoff;
[psth_mean,psth_mean_grp] = ap.groupfun(@mean,striatum_psth(use_trials,:), ...
    [striatum_psth_grp.animal(use_trials),plot_day_grp(use_trials), ...
    split_idx(use_trials),striatum_psth_grp.kidx(use_trials)]);

max_t = psth_t > 0 & psth_t < 0.2;
psth_max = max(psth_mean(:,max_t),[],2);
[psth_max_avg,psth_max_grp] = ap.nestgroupfun({@mean,@mean},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2:end));
psth_max_sem = ap.nestgroupfun({@mean,@AP_sem},psth_max,psth_mean_grp(:,1),psth_mean_grp(:,2:end));

figure;
h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day = 1:length(plot_day_bins)-1
        nexttile; hold on;
        curr_data_idx = psth_max_grp(:,1) == curr_day & ...
            psth_max_grp(:,3) == curr_k;
        errorbar(psth_max_avg(curr_data_idx),psth_max_sem(curr_data_idx),'k','linewidth',2);
        axis padded
        if curr_k == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Fig 2X] Cortex task trial heatmap (reaction-sorted)

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% Plot heatmaps sorted by reaction times
figure; tiledlayout(n_k,max(plot_day_grp),'TileSpacing','none');
for curr_k = 1:n_k
    for curr_day = unique(plot_day_grp)'
        nexttile;
        curr_trials = find(plot_day_grp == curr_day);
       
        [~,sort_idx] = sort(wf_grp.rxn(curr_trials));
        imagesc(wf_t,[],movmean(wf_striatum_roi(curr_trials(sort_idx),:,curr_k),[20,5]));
        colormap(ap.colormap('PWG'));
        clim(max(abs(clim)).*[-1,1]);
        axis off;
        if curr_k == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end

ap.prettyfig;


%% [Fig 2X] Cortex task average PSTH w/ and w/o movement

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% Plot average activity in trial (NaN-out movement times)
figure; h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','none');
move_leeway = 0.1;
for curr_k = 1:n_k
    for curr_day = unique(plot_day_grp)'
        nexttile; hold on;
        curr_trials = find(plot_day_grp == curr_day);      
        
        % Average all data
        curr_data = wf_striatum_roi(curr_trials,:,curr_k);

        curr_data_mean = mean(ap.groupfun(@nanmean,curr_data,wf_grp.animal(curr_trials)),1);
        curr_data_sem = AP_sem(ap.groupfun(@nanmean,curr_data,wf_grp.animal(curr_trials)),1);

        % Average data with movement removede
        % (only plot no-move timepoints with at least N trials and N animals)
        curr_data_no_move = wf_striatum_roi(curr_trials,:,curr_k).* ...
            AP_nanout(wf_t > wf_grp.rxn(curr_trials)-move_leeway);

        min_nomove_trials = 10;
        min_noanimals_trials = 5;

        curr_data_no_move_n = ap.groupfun(@(x) sum(~isnan(x)),curr_data_no_move,wf_grp.animal(curr_trials));
        curr_data_no_move_animal_mean = ap.groupfun(@nanmean,curr_data_no_move,wf_grp.animal(curr_trials));
        curr_data_no_move_animal_mean_mindata = curr_data_no_move_animal_mean.* ...
            AP_nanout(curr_data_no_move_n < min_nomove_trials).* ...
            AP_nanout(sum(curr_data_no_move_n >= min_nomove_trials) < min_noanimals_trials);
        curr_data_no_move_mean = nanmean(curr_data_no_move_animal_mean_mindata,1);
        curr_data_no_move_sem = AP_sem(curr_data_no_move_animal_mean_mindata,1);

        % Plot
        ap.errorfill(wf_t,curr_data_mean,curr_data_sem,'k',0.3,false,2);
        plot(wf_t,curr_data_no_move_mean,'k','linewidth',2);
        axis off;
        if curr_k == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Fig 2X] Cortex task max stim activity split aross session




%% [Fig 3X] Striatum PSTH passive

plot_day_bins = [-Inf,-1:1,Inf];
plot_days_grp = discretize(max(striatum_psth_grp.ld,-inf),plot_day_bins);
day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

[striatum_psth_avg,striatum_psth_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean}, ...
    striatum_psth,striatum_psth_grp.animal, ...
    [plot_days_grp,striatum_psth_grp.stim,striatum_psth_grp.kidx]);
striatum_psth_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    striatum_psth,striatum_psth_grp.animal, ...
    [plot_days_grp,striatum_psth_grp.stim,striatum_psth_grp.kidx]);

unique_stim = unique(striatum_psth_grp.stim);
figure;
h = tiledlayout(n_k,max(plot_days_grp)*length(unique_stim),'TileSpacing','tight');
for curr_k = unique(striatum_psth_avg_grp(:,3))'
    for curr_day = unique(striatum_psth_avg_grp(:,1))'
        for plot_stim = unique_stim'
            nexttile; axis off;
            plot_data = striatum_psth_avg_grp(:,1) == curr_day & ...
                striatum_psth_avg_grp(:,2) == plot_stim & ...
                striatum_psth_avg_grp(:,3) == curr_k;
            if ~any(plot_data)
                continue
            end
            ap.errorfill(psth_t,striatum_psth_avg(plot_data,:),striatum_psth_sem(plot_data,:), ...
                day_colormap(curr_day,:));
        end
    end
end
linkaxes(h.Children,'xy');
xlim([-0.2,0.8]);
ap.prettyfig;


%% [Fig 3X] Striatum PSTH max passive

plot_day_bins = [-Inf,-2:2,Inf];
plot_days_grp = discretize(max(striatum_psth_grp.ld,-inf),plot_day_bins);
day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

[striatum_psth_dayavg,striatum_psth_dayavg_grp] = ...
    ap.groupfun(@mean,striatum_psth, ...
    [striatum_psth_grp.animal,plot_days_grp,striatum_psth_grp.stim,striatum_psth_grp.kidx]);

max_t = psth_t > 0 & psth_t < 0.2;
[striatum_psth_max,striatum_psth_max_grp] = ap.nestgroupfun({@mean,animal_groupfun}, ...
    max(striatum_psth_dayavg(:,max_t),[],2),striatum_psth_dayavg_grp(:,1), ...
    striatum_psth_dayavg_grp(:,2:end));
striatum_psth_max_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    max(striatum_psth_dayavg(:,max_t),[],2),striatum_psth_dayavg_grp(:,1), ...
    striatum_psth_dayavg_grp(:,2:end));

figure;
h = tiledlayout(n_k,1);
binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');
for curr_k = unique(striatum_psth_grp.kidx)'
    nexttile; hold on; set(gca,'ColorOrder',ap.colormap('BKR',3));
    for curr_stim = unique(striatum_psth_grp.stim)'
        plot_grps = striatum_psth_max_grp(:,2) == curr_stim & ...
            striatum_psth_max_grp(:,3) == curr_k;

        errorbar(binned_days_x,striatum_psth_max(plot_grps), ...
            striatum_psth_max_sem(plot_grps),'linewidth',2);
    end
    axis padded
    xline(0);
end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Fig 3X] Cortex passive max

plot_day_bins = [-Inf,-1:1,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

stim_t = wf_t > 0 & wf_t < 0.2;
[wf_avg,wf_avg_grp] = ap.groupfun(@mean,cell2mat(wf.V_no_move_stim_align), ...
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
        clim([-1,1]*1e-3*2);
        axis image off;
        ap.wf_draw('ccf',[0.5,0.5,0.5]);
    end
end
ap.prettyfig;

%% [Fig 3X] Cortex ROIs

figure;tiledlayout(n_k,1,'tilespacing','none')
for curr_k = 1:n_k
    nexttile;
    imagesc(striatum_wf_roi(:,:,curr_k));
    axis image off; ap.wf_draw('ccf',[0.5,0.5,0.5]);
    colormap(ap.colormap('WG'));
end
ap.prettyfig;


%% [Fig 3X] Cortex passive ROI average

plot_day_bins = [-Inf,-1:1,Inf];
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
h = tiledlayout(n_k,max(plot_day_grp)*length(unique_stim),'TileSpacing','tight');
for curr_k = 1:size(striatum_wf_roi,3)
    for curr_ld = unique(plot_day_grp)'
        for curr_stim = unique_stim'
            nexttile; axis off;
            curr_data_idx = wf_striatum_roi_avg_grp(:,1) == curr_ld & ...
                wf_striatum_roi_avg_grp(:,2) == curr_stim;

            ap.errorfill(wf_t,wf_striatum_roi_avg(curr_data_idx,:,curr_k), ...
                wf_striatum_roi_sem(curr_data_idx,:,curr_k), ...
                day_colormap(curr_ld,:));
        end
    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Fig 3X] Cortex passive ROI max

plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

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
h = tiledlayout(n_k,1,'TileSpacing','tight');
binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');
for curr_k = 1:size(striatum_wf_roi,3)
    nexttile; hold on;
    set(gca,'ColorOrder',ap.colormap('BKR',3));
    for curr_stim = unique(unique(wf_striatum_roi_max_avg_grp(:,2)))'
        curr_data_idx = wf_striatum_roi_max_avg_grp(:,2) == curr_stim;

        errorbar(binned_days_x,wf_striatum_roi_max_avg(curr_data_idx,:,curr_k), ...
            wf_striatum_roi_max_sem(curr_data_idx,:,curr_k),'linewidth',2);
        ylabel('\DeltaF/F_0 max');
        axis padded
    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;





