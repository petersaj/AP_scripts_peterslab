%% Notes

% Putting together paper figs

%% Behavior v2
% (this one includes other stats - choose how to define learning)

use_stat = 'firstmove_mad';
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

% K-means cluster maps
maps_flat_cat = reshape(cat(3,all_ctx_maps_to_str.cortex_kernel_px{:}),prod(U_size),[]);
expl_var_cat = vertcat(all_ctx_maps_to_str.explained_var{:});

n_k = 6;
kmeans_starting = nanmean(cell2mat(permute(cellfun(@(x) ...
    x(:,:,round(linspace(1,size(x,3),n_k))), ...
    all_ctx_maps_to_str.cortex_kernel_px(~cellfun(@isempty, ...
    all_ctx_maps_to_str.cortex_kernel_px)),'uni',false),[2,3,4,1])),4);

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
% [kidx,kmeans_centroid_flat,~,kmeans_dist] = kmeans(...
%     (maps_flat_cat> prctile(maps_flat_cat,50,1))',n_k, ...
%     'Distance','hamming','start',reshape(+kmeans_starting_binary,[],n_k)');
% kmeans_centroid = reshape(kmeans_centroid_flat',size(kmeans_starting,1),size(kmeans_starting,2),[]);
% %%%%%%

[kidx,kmeans_centroid_flat,~,kmeans_dist] = kmeans(...
    maps_flat_cat',n_k, ...
    'Distance','cosine','start',reshape(kmeans_starting,[],n_k)');
kmeans_centroid = reshape(kmeans_centroid_flat',size(kmeans_starting,1),size(kmeans_starting,2),[]);

% %%%%%%%%%%%%%%%%%%%
% % ADJUSTMENTS FOR N_K=4
% % looks like value of ~0.75 works, just picked empirical way to get there
% kidx_raw = kidx;
% 
% dist_cutoff_k1 = prctile(kmeans_dist(kidx_raw==1,1),95);
% dist_cutoff_k2 = prctile(kmeans_dist(kidx_raw==2,2),95);
% 
% kidx = kidx_raw;
% 
% % (to fold 2's into cluster 1)
% kidx(kidx_raw==2 & kmeans_dist(:,1) < dist_cutoff_k2) = 1;
% 
% % % (make 2's into new cluster)
% % kidx(kidx_raw==2 & kmeans_dist(:,1) < dist_cutoff) = 5;
% % [~,kidx] = ismember(kidx,[1,5,2,3,4]); kidx(kidx==0) = NaN;% re-assign to order depths
% % n_k = 5;
% % % kidx(kidx==3) = 2; kidx(kidx==4) = 3;n_k = 3; % (combine pfc 2+3)
% 
% % % (make 1's and 2's into new cluster)
% % kidx(kidx_raw==1 & kmeans_dist(:,2) < dist_cutoff_k1) = 5;
% % kidx(kidx_raw==2 & kmeans_dist(:,1) < dist_cutoff_k2) = 5;
% % [~,kidx] = ismember(kidx,[1,5,2,3,4]); kidx(kidx==0) = NaN;% re-assign to order depths
% % n_k = 5;
% 
% %%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%
% % ADJUSTMENTS FOR N_K=6
% Turn 6 into 3 (vis, frontal, lateral)
kidx(ismember(kidx,2:4)) = 2;
kidx(ismember(kidx,5:6)) = 3;
n_k = 3;
% %%%%%%%%%%%%%%%%%%%


kmeans_cluster_mean = reshape(ap.groupfun(@nanmean,maps_flat_cat,[],kidx),[U_size,n_k]);

kidx_rec = mat2cell(kidx,cellfun(@(x) size(x,3).*(size(x,1)>0), ...
    all_ctx_maps_to_str.cortex_kernel_px));

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
% for plot_kidx = 1:n_k
%     plot(squeeze(kernel_max_x(kidx==plot_kidx)),squeeze(kernel_max_y(kidx==plot_kidx)),'x','linewidth',2.5);
% end
% ap.wf_draw('ccf','k');
% legend(string(1:n_k))


%% Load ephys

% AM-packaged: 
data_path = fullfile(plab.locations.server_path,'Users','Andrada-Maria_Marica','long_str_ctx_data');

% % Added move-align: 
% data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data\AM_packaged\stim\';

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