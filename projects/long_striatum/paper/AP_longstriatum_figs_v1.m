%% Figures for longitudinal striatum+cortex project
%
% Data packaged by AM save scripts
% Analysis taken fom `AP_longstriatum_exploratory_2`
%
% Uses `AP_longstriatum_load_data` to load/prep relevant data at the start
% of each figure (doesn't re-load if already overloaded - this lets each
% panel be run independently, or for the entire script to be run straight
% through)

tic

%% [Fig 1X] Reaction stat and learning histogram

%%% Load data for figure
load_dataset = 'noact';
AP_longstriatum_load_data;
%%%

% Plot reaction time and association index, split within day
n_daysplit = 3;
use_rxn = cellfun(@(x) ~isnan(x),bhv.stim_to_move_nullmean,'uni',false);

rxn_mean_daysplit = cell2mat(cellfun(@(x,idx) ap.groupfun(@mean,x(idx), ...
    ap.quantile_bin(sum(idx),n_daysplit)),bhv.stim_to_move,use_rxn,'uni',false)')';

rxn_null_mean_daysplit = cell2mat(cellfun(@(x,idx) ap.groupfun(@mean,x(idx), ...
    ap.quantile_bin(sum(idx),n_daysplit)),bhv.stim_to_move_nullmean,use_rxn,'uni',false)')';

rxn_idx_daysplit = (rxn_null_mean_daysplit-rxn_mean_daysplit)./ ...
    (rxn_null_mean_daysplit+rxn_mean_daysplit);

[rxn_daysplit_mean,rxn_group_x] = ap.groupfun(@mean,rxn_mean_daysplit, ...
    bhv.days_from_learning);
rxn_null_daysplit_mean = ap.groupfun(@mean,rxn_null_mean_daysplit, ...
    bhv.days_from_learning);

rxn_idx_daysplit_mean = ap.groupfun(@mean,rxn_idx_daysplit, ...
    bhv.days_from_learning);
rxn_idx_daysplit_sem = ap.groupfun(@AP_sem,rxn_idx_daysplit, ...
    bhv.days_from_learning);

plot_days = -3:2;
plot_day_idx = ismember(rxn_group_x,plot_days);

figure; tiledlayout(2,1);
rxn_group_x_daysplit = rxn_group_x+(0:n_daysplit)./n_daysplit;

nexttile; hold on; set(gca,'YScale','log');
plot(reshape(rxn_group_x_daysplit(plot_day_idx,:)',[],1), ...
    reshape(padarray(rxn_daysplit_mean(plot_day_idx,:),[0,1],nan,'post')',[],1),'k','linewidth',2);
plot(reshape(rxn_group_x_daysplit(plot_day_idx,:)',[],1), ...
    reshape(padarray(rxn_null_daysplit_mean(plot_day_idx,:),[0,1],nan,'post')',[],1),'r','linewidth',2);
xline(0,'r');
ylabel('Reaction time');
xlabel('Day from learning');

nexttile;
errorbar(reshape(rxn_group_x_daysplit(plot_day_idx,:)',[],1), ...
    reshape(padarray(rxn_idx_daysplit_mean(plot_day_idx,:),[0,1],nan,'post')',[],1), ...
    reshape(padarray(rxn_idx_daysplit_sem(plot_day_idx,:),[0,1],nan,'post')',[],1),'k','linewidth',2);
xline(0,'r');
ylabel('Association index');
xlabel('Day from learning');
ap.prettyfig;


% Plot histogram of learning days
n_learned_day = cellfun(@(x) max([0, ...
    find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'));
figure;
histogram(n_learned_day,-0.5:max(n_learned_day)+0.5,'FaceColor','k','EdgeColor','none');
ylabel('Number of mice');
xlabel('Days to learn');
ap.prettyfig;


%% [Fig 1X] Example striatal units and corticostriatal maps

%%% Load data for figure
load_dataset = 'noact';
AP_longstriatum_load_data;
%%%

% Choose animal and day to plot
use_animal = 'AM026';

animal_days = find(strcmp(bhv.animal,use_animal));
plot_days = [1,length(animal_days)];

figure;
h = tiledlayout(2,2);
title(h,sprintf('%s, days %d,%d',use_animal,plot_days));
for curr_rec_idx = plot_days

    % Load data from example recording
    use_rec = animal_days(curr_rec_idx);
    animal = bhv.animal{use_rec};
    rec_day = bhv.rec_day{use_rec};
    recordings = plab.find_recordings(animal,rec_day,'stim_wheel*');
    rec_time = recordings.recording{end};
    load_parts.ephys = true;
    ap.load_recording;

    % Plot widefield average (last day recording only)
    if curr_rec_idx == plot_days(end)
        wf_day_path = plab.locations.filename('server',animal,rec_day,[],'widefield');
        mean_image_fn = fullfile(wf_day_path,sprintf('meanImage_blue.npy'));
        wf_avg = plab.wf.wf_align(readNPY(mean_image_fn),animal,rec_day);
        figure;imagesc(wf_avg);
        ap.wf_draw('cortex','y');
        axis image off;
        ap.prettyfig;
    end

    % Plot units and overlay clustering
    domain_color = {'R','G','B'};
    domain_color_rgb = [1,0,0;0,1,0;0,0,1];

    ax = nexttile(h); hold on;

    domain_im = permute(domain_color_rgb(domain_idx_rec{use_rec},:),[1,3,2]);
    imagesc(ax,[],ctx_str_maps.depth_group_edges{use_rec},domain_im);
    ax.YDir = 'reverse';

    ap.plot_unit_depthrate(spike_times_timelite,spike_templates,template_depths,[],ax)
    yline(ctx_str_maps.depth_group_edges{use_rec},'linewidth',2,'color',[0.5,0.5,0.5]);

    % Plot average domain map colored and combined
    domain_avg = ap.groupfun(@mean,ctx_str_maps.cortex_striatum_map{use_rec},[],[],domain_idx_rec{use_rec});

    col_lim = [0,0.01];
    domain_colored = nan([size(domain_avg),3]);
    for curr_domain = 1:n_domains
        curr_colormap = ap.colormap(['W',domain_color{curr_domain}],[],2);
        curr_map_gray = 1+round(mat2gray(domain_avg(:,:,curr_domain),col_lim).*(size(curr_colormap,1)-1));
        domain_colored(:,:,curr_domain,:) = reshape(curr_colormap(curr_map_gray,:),size(curr_map_gray,1),size(curr_map_gray,2),3,[]);
    end

    domain_colored_combined = squeeze(min(domain_colored,[],3));
    
    nexttile(h);
    image(domain_colored_combined)
    axis image off;
    ap.wf_draw('cortex',[0.5,0.5,0.5]);

end

ap.prettyfig;


%% [Fig 1X] Corticostriatal map pre/post learning

%%% Load data for figure
load_dataset = 'noact';
AP_longstriatum_load_data;
%%%

% Make map indicies (not done in load data)
wf_map_grp = struct;
wf_map_grp.ld = cell2mat(cellfun(@(idx,maps) repmat(idx,size(maps,3).*~isempty(maps),1), ...
    num2cell(bhv.days_from_learning),ctx_str_maps.cortex_striatum_map,'uni',false));
wf_map_grp.animal = cell2mat(cellfun(@(idx,maps) repmat(idx,size(maps,3).*~isempty(maps),1), ...
    num2cell(grp2idx(bhv.animal)),ctx_str_maps.cortex_striatum_map,'uni',false));

% Group pre/post learn days
plot_day_bins = [-Inf,-0,Inf];
plot_day_grp = discretize(max(wf_map_grp.ld,-inf),plot_day_bins);

% Get and plot average maps
[wf_map_avg,wf_map_avg_grp] = ap.nestgroupfun({@mean,@mean}, ...
    reshape(cat(3,ctx_str_maps.cortex_striatum_map{:}),prod(U_size),[])', ...
    wf_map_grp.animal,[plot_day_grp,domain_idx]);

figure;
domain_color = {'R','G','B'};
h = tiledlayout(n_domains,length(plot_day_bins)-1,'TileSpacing','none');
for curr_domain = 1:n_domains
    for curr_day = 1:length(plot_day_bins)-1
        curr_data_idx = ismember(wf_map_avg_grp,[curr_day,curr_domain],'rows');
        
        nexttile;
        imagesc(reshape(wf_map_avg(curr_data_idx,:),size(U_master,[1,2])));
        axis image off;
        colormap(gca,ap.colormap(['W',domain_color{curr_domain}],[],2));
        clim([0,0.01]);
        ap.wf_draw('cortex',[0.5,0.5,0.5]);
    end
end

ap.prettyfig;


%% [Fig 1X] Striatal MUA pre/post learning

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,0,Inf];
plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
learn_colormap = ap.colormap('BKR',3);
prepost_colormap = max(0,learn_colormap([1,end],:)-0.2);

% Plot average task activity (pre/post learning)
stim_x = [-0.2,0.15];
move_x = [-0.05,0.4];
outcome_x = [-0.1,0.5];

figure; h = tiledlayout(n_domains,3,'TileSpacing','tight');
for curr_depth = 1:n_domains

    curr_trials = striatum_mua_grp.domain_idx == curr_depth;

    % Stim
    nexttile; hold on; set(gca,'ColorOrder',prepost_colormap);
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean}, ...
        striatum_mua(curr_trials,:,1),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem}, ...
        striatum_mua(curr_trials,:,1),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(stim_x);
    xline(0);

    % Move
    nexttile; hold on; set(gca,'ColorOrder',prepost_colormap);
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean}, ...
        striatum_mua(curr_trials,:,2),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem}, ...
        striatum_mua(curr_trials,:,2),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(move_x);
    xline(0);

    % Outcome
    nexttile; hold on; set(gca,'ColorOrder',prepost_colormap);
    curr_data_mean = ap.nestgroupfun({@nanmean,@nanmean}, ...
        striatum_mua(curr_trials,:,3),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    curr_data_sem = ap.nestgroupfun({@nanmean,@AP_sem}, ...
        striatum_mua(curr_trials,:,3),striatum_mua_grp.animal(curr_trials),plot_day_grp(curr_trials));
    ap.errorfill(psth_t,curr_data_mean',curr_data_sem');
    xlim(outcome_x);
    xline(0);

end
% (link all y, and x of same-alignment)
linkaxes(h.Children,'y');
for ax = 1:3
    linkaxes(h.Children(ax:3:end),'x');
end
% (set same data aspect to have same x-size)
[h.Children.DataAspectRatio] = deal(min(vertcat(h.Children.DataAspectRatio),[],1));
ap.prettyfig;


%% [Fig 2X] Cortex + striatum task trial heatmap (reaction-sorted)

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,Inf];
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);


% Plot heatmaps sorted by reaction times
heatmap_smooth = [20,1]; % ([trials,time] to smooth for graphics)

figure; 
h = tiledlayout(n_domains*2,max(cortex_plot_day_grp), ...
    'TileIndexing','ColumnMajor','TileSpacing','Tight');
for curr_day = unique(cortex_plot_day_grp)'
    for curr_domain = 1:n_domains

        % Cortex
        nexttile;
        curr_trials = find(cortex_plot_day_grp == curr_day);
       
        [sorted_rxn,sort_idx] = sort(wf_grp.rxn(curr_trials));
        imagesc(wf_t,[],movmean(wf_striatum_roi(curr_trials(sort_idx),:,curr_domain,1),heatmap_smooth));
        colormap(gca,ap.colormap('WG'));
        clim([0,1e-2]);
        axis off;

        hold on
        xline(0,'color','r');
        plot(sorted_rxn,1:length(curr_trials),'b');

        if curr_domain == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end

        % Striatum
        nexttile;
        curr_trials = find(striatum_plot_day_grp == curr_day & striatum_mua_grp.domain_idx == curr_domain);
      
        [sorted_rxn,sort_idx] = sort(striatum_mua_grp.rxn(curr_trials));      
        imagesc(psth_t,[],movmean(striatum_mua(curr_trials(sort_idx),:,1),heatmap_smooth));
        colormap(gca,ap.colormap('WK'));
        clim([0,3]);
        axis off;

        hold on
        xline(0,'color','r');
        plot(sorted_rxn,1:length(curr_trials),'b');

    end
end
linkaxes(h.Children,'x');
xlim(h.Children,[-0.5,1])
ap.prettyfig;


%% [Fig 2X] Cortex + striatum task average PSTH w/o movement

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-1,0,Inf];
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% NaN-out activity after movement onset (minus leeway time)
move_leeway = 0.1; % time pre-movement to exclude
striatum_mua_nomove = striatum_mua(:,:,1).*ap.nanout(psth_t > striatum_mua_grp.rxn-move_leeway);
wf_striatum_roi_nomove = wf_striatum_roi(:,:,:,1).*ap.nanout(wf_t > wf_grp.rxn-move_leeway);

% Get average and SEM no-movement activity
[striatum_mua_nomove_avg,striatum_mua_nomove_avg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},striatum_mua_nomove, ...
    striatum_mua_grp.animal,[striatum_plot_day_grp,striatum_mua_grp.domain_idx]);
striatum_mua_nomove_sem = ...
    ap.nestgroupfun({@nanmean,@AP_sem},striatum_mua_nomove, ...
    striatum_mua_grp.animal,[striatum_plot_day_grp,striatum_mua_grp.domain_idx]);

[wf_striatum_roi_nomove_avg,wf_striatum_roi_nomove_avg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},wf_striatum_roi_nomove, ...
    wf_grp.animal,cortex_plot_day_grp);
wf_striatum_roi_nomove_sem = ...
    ap.nestgroupfun({@nanmean,@AP_sem},wf_striatum_roi_nomove, ...
    wf_grp.animal,cortex_plot_day_grp);

str_day_color = ap.colormap('KW',length(plot_day_bins)-1);
ctx_day_color = ap.colormap('KG',length(plot_day_bins)-1);
figure; h = tiledlayout(n_domains*2,1,'TileSpacing','tight');
for curr_domain = 1:n_domains

    nexttile; hold on; set(gca,'ColorOrder',ctx_day_color);
    ap.errorfill(wf_t,wf_striatum_roi_nomove_avg(:,:,curr_domain)', ...
        wf_striatum_roi_nomove_sem(:,:,curr_domain)');

    nexttile; hold on; set(gca,'ColorOrder',str_day_color);
    plot_data_idx = striatum_mua_nomove_avg_grp(:,2) == curr_domain;
    ap.errorfill(psth_t,striatum_mua_nomove_avg(plot_data_idx,:)', ...
        striatum_mua_nomove_sem(plot_data_idx,:)');

end

linkaxes(h.Children(1:2:end));
linkaxes(h.Children(2:2:end));
xlim(h.Children,[-0.1,0.5]);
ap.prettyfig;

% ~stats~
fprintf('---- STATS ----\n')
compare_day_grps = [1,2];

stim_t = [0,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

[wf_striatum_roi_nomove_animal,wf_striatum_roi_nomove_animal_grp] = ...
    ap.groupfun(@nanmean,wf_striatum_roi_nomove, ...
    [wf_grp.animal,cortex_plot_day_grp]);
wf_striatum_roi_nomove_animal_tmax = max(wf_striatum_roi_nomove_animal(:,cortex_stim_t,:),[],2);

cortex_stat_usedata = ismember(wf_striatum_roi_nomove_animal_grp(:,2),compare_day_grps);
cortex_stat_meas = permute(diff(ap.groupfun(@mean,wf_striatum_roi_nomove_animal_tmax(cortex_stat_usedata,:,:), ...
    wf_striatum_roi_nomove_animal_grp(cortex_stat_usedata,2)),[],1),[3,2,1]);

[striatum_mua_nomove_animal,striatum_mua_nomove_animal_grp] = ...
    ap.groupfun(@nanmean,striatum_mua_nomove, ...
    [striatum_mua_grp.animal,striatum_plot_day_grp,striatum_mua_grp.domain_idx]);
striatum_mua_nomove_animal_tmax = max(striatum_mua_nomove_animal(:,striatum_stim_t),[],2);

striatum_stat_usedata = ismember(striatum_mua_nomove_animal_grp(:,2),compare_day_grps);
striatum_stat_meas = ap.nestgroupfun({@mean,@diff},striatum_mua_nomove_animal_tmax(striatum_stat_usedata), ...
    striatum_mua_nomove_animal_grp(striatum_stat_usedata,2),striatum_mua_nomove_animal_grp(striatum_stat_usedata,3));

n_shuff = 10000;
cortex_stat_null = nan(n_domains,n_shuff);
striatum_stat_null = nan(n_domains,n_shuff);
for curr_shuff = 1:n_shuff

    curr_ctx_shuff = ...
        ap.shake(wf_striatum_roi_nomove_animal_tmax(cortex_stat_usedata,:,:),1);
    cortex_stat_null(:,curr_shuff) = ...
        permute(diff(ap.groupfun(@mean,curr_ctx_shuff, ...
        wf_striatum_roi_nomove_animal_grp(cortex_stat_usedata,2)),[],1),[3,2,1]);

    curr_str_shuff = ...
        ap.shake(striatum_mua_nomove_animal_tmax(striatum_stat_usedata),1, ...
        striatum_mua_nomove_animal_grp(striatum_stat_usedata,3));
    striatum_stat_null(:,curr_shuff) = ...
        ap.nestgroupfun({@mean,@diff},curr_str_shuff, ...
        striatum_mua_nomove_animal_grp(striatum_stat_usedata,2),striatum_mua_nomove_animal_grp(striatum_stat_usedata,3));

end

fprintf('Day grps %d,%d:\n',compare_day_grps);

cortex_stat_rank = tiedrank([cortex_stat_meas,cortex_stat_null]')';
cortex_stat_p = 1-cortex_stat_rank(:,1)/(n_shuff+1);
cortex_stat_sig = discretize(cortex_stat_p < 0.05,[0,1,Inf],{'','*'});
for curr_domain = 1:n_domains
    fprintf('CTX%d p = %.2g%s\n',curr_domain,cortex_stat_p(curr_domain),cortex_stat_sig{curr_domain});
end

striatum_stat_rank = tiedrank([striatum_stat_meas,striatum_stat_null]')';
striatum_stat_p = 1-striatum_stat_rank(:,1)/(n_shuff+1);
striatum_stat_sig = discretize(striatum_stat_p < 0.05,[0,1,Inf],{'','*'});
for curr_domain = 1:n_domains
    fprintf('STR%d p = %.2g%s\n',curr_domain,striatum_stat_p(curr_domain),striatum_stat_sig{curr_domain});
end


%% [Fig 2X] Cortex + striatum task max, split within session

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

rxn_cutoff = 0.3; % only plot trials with slow reaction times

plot_day_bins = [-Inf,-2:1];
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% Split by trial percentile within day
n_split = 4;

striatum_split_idx = cell2mat(cellfun(@(n_trials,n_mua) ...
    reshape(repmat(ap.quantile_bin(n_trials,n_split),1,n_mua),[],1), ...
    num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

cortex_split_idx = cell2mat(cellfun(@(n_trials) ...
    ap.quantile_bin(n_trials,n_split), ...
    num2cell(cortex_trials_rec_n),'uni',false));

% For bins that contain multiple days: don't split within day
% (doesn't makes sense to average splits across days)
striatum_split_idx(ismember(striatum_plot_day_grp,find(diff(plot_day_bins)~=1))) = 1;
cortex_split_idx(ismember(cortex_plot_day_grp,find(diff(plot_day_bins)~=1))) = 1;

% Get activity: trial mean > time max > animal mean 
stim_t = [0,0.2];

% (striatum)
striatum_use_trials = striatum_mua_grp.rxn > rxn_cutoff;
[striatum_activity_mean,striatum_activity_mean_grp] = ap.groupfun(@mean,striatum_mua(striatum_use_trials,:,1), ...
    [striatum_mua_grp.animal(striatum_use_trials),striatum_plot_day_grp(striatum_use_trials), ...
    striatum_split_idx(striatum_use_trials),striatum_mua_grp.domain_idx(striatum_use_trials)]);

striatum_activity_max = max(striatum_activity_mean(:,isbetween(psth_t,stim_t(1),stim_t(2))),[],2);
[striatum_activity_max_avg,striatum_activity_max_grp] = ap.nestgroupfun({@mean,@mean},striatum_activity_max,striatum_activity_mean_grp(:,1),striatum_activity_mean_grp(:,2:end));
striatum_activity_max_sem = ap.nestgroupfun({@mean,@AP_sem},striatum_activity_max,striatum_activity_mean_grp(:,1),striatum_activity_mean_grp(:,2:end));

% (cortex)
cortex_use_trials = wf_grp.rxn > rxn_cutoff;
[cortex_activity_mean,cortex_activity_mean_grp] = ap.groupfun(@mean,wf_striatum_roi(cortex_use_trials,:,:,1), ...
    [wf_grp.animal(cortex_use_trials),cortex_plot_day_grp(cortex_use_trials),cortex_split_idx(cortex_use_trials)]);

cortex_activity_max = permute(max(cortex_activity_mean(:,isbetween(wf_t,stim_t(1),stim_t(2)),:),[],2),[1,3,2]);

[cortex_activity_max_avg,cortex_activity_max_grp] = ap.nestgroupfun({@mean,@mean},cortex_activity_max,cortex_activity_mean_grp(:,1),cortex_activity_mean_grp(:,2:end));
cortex_activity_max_sem = ap.nestgroupfun({@mean,@AP_sem},cortex_activity_max,cortex_activity_mean_grp(:,1),cortex_activity_mean_grp(:,2:end));

% Plot cortex and striatum overlaid
figure;
h = tiledlayout(n_domains,length(plot_day_bins)-1,'TileSpacing','compact');
for curr_domain = 1:n_domains
    for curr_day = 1:length(plot_day_bins)-1
        nexttile; hold on;

        yyaxis left;
        curr_cortex_data_idx = cortex_activity_max_grp(:,1) == curr_day;
        if sum(curr_cortex_data_idx) > 1
            ap.errorfill(1:sum(curr_cortex_data_idx),cortex_activity_max_avg(curr_cortex_data_idx,curr_domain), ...
                cortex_activity_max_sem(curr_cortex_data_idx),[0,0.6,0],0.5,false,2);
        else
            errorbar(cortex_activity_max_avg(curr_cortex_data_idx,curr_domain), ...
                cortex_activity_max_sem(curr_cortex_data_idx),'color',[0,0.6,0],'linewidth',4, ...
                'CapSize',0);
        end

        yyaxis right;
        curr_striatum_data_idx = striatum_activity_max_grp(:,1) == curr_day & ...
            striatum_activity_max_grp(:,3) == curr_domain;
        errorbar(striatum_activity_max_avg(curr_striatum_data_idx), ...
            striatum_activity_max_sem(curr_striatum_data_idx),'k','linewidth',2, ...
            'CapSize',0);

        if curr_domain == 1
            title(sprintf('%d:%d',plot_day_bins(curr_day),plot_day_bins(curr_day+1)));
        end
    end
end

% Set axis limits
for curr_ax = 1:length(h.Children)
    yyaxis(h.Children(curr_ax),'left');
end
linkaxes(h.Children);
ylim(h(1).Children,[0,8e-3]);
[h.Children.YColor] = deal([0,0.5,0]);

for curr_ax = 1:length(h.Children)
    yyaxis(h.Children(curr_ax),'right');
end
linkaxes(h.Children);
ylim(h(1).Children,[-0.2,8]);
[h.Children.YColor] = deal([0,0,0]);

xlim(h(1).Children,[0.8,n_split+0.2])

ap.prettyfig;

% ~stats~
fprintf('---- STATS ----\n')
fprintf('1-way ANOVA:');
for curr_domain = 1:n_domains
    for curr_day = 1:length(plot_day_bins)-1
        stat_data_idx = ismember(striatum_activity_mean_grp(:,[2,4]),[curr_day,curr_domain],'rows');
        p = anovan(striatum_activity_max(stat_data_idx),striatum_activity_mean_grp(stat_data_idx,3),'display','off');
        stat_sig = discretize(p < 0.05,[0,1,Inf],["","*"]);
        fprintf('STR %d, day grp %d, p = %.2g%s\n',curr_domain,curr_day,p,stat_sig)
    end
end
for curr_domain = 1:n_domains
    for curr_day = 1:length(plot_day_bins)-1
        stat_data_idx = cortex_activity_mean_grp(:,2) == curr_day;
        p = anovan(cortex_activity_max(stat_data_idx,curr_domain),cortex_activity_mean_grp(stat_data_idx,3),'display','off');
        stat_sig = discretize(p < 0.05,[0,1,Inf],["","*"]);
        fprintf('CTX %d, day grp %d, p = %.2g%s\n',curr_domain,curr_day,p,stat_sig)
    end
end


%% [Fig 3X] Striatum PSTH passive

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,2,Inf];
plot_days_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
day_colormap = unique(vertcat(flipud(ap.colormap('KB',sum(plot_day_bins(1:end-1)<=0))), ...
     ap.colormap('KR',sum(plot_day_bins(1:end-1)>=0))),'stable','rows');

[striatum_mua_avg,striatum_mua_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);
striatum_mua_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);

unique_stim = unique(striatum_mua_grp.stim);
figure;
h = tiledlayout(n_domains,max(plot_days_grp)*length(unique_stim),'TileSpacing','tight');
for curr_domain = unique(striatum_mua_avg_grp(:,3))'
    for curr_day = unique(striatum_mua_avg_grp(:,1))'
        for plot_stim = unique_stim'
            nexttile; axis off;
            plot_data = striatum_mua_avg_grp(:,1) == curr_day & ...
                striatum_mua_avg_grp(:,2) == plot_stim & ...
                striatum_mua_avg_grp(:,3) == curr_domain;
            if ~any(plot_data)
                continue
            end
            ap.errorfill(psth_t,striatum_mua_avg(plot_data,:),striatum_mua_sem(plot_data,:), ...
                day_colormap(curr_day,:));
        end
    end
end
linkaxes(h.Children,'xy');
xlim([-0.2,0.8]);
title(h,'Striatum');
ap.prettyfig;





%% [Fig 3X] Striatum fraction stim-responsive units

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

% Set days to group
plot_day_bins = [-Inf,-2:2,Inf];
plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

% Grab "responsive" striatal units
unit_responsive_p_thresh = 0.95;
unit_responsive = cell2mat(horzcat(ephys.unit_resp_p_value{:})') > unit_responsive_p_thresh;
striatum_units_responsive = unit_responsive(cell2mat(striatum_units),:);

% Get fraction of responsive unit (restrict to MSNs)
use_units = striatum_sua_grp.msn;

[unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
    +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
    [plot_day_grp(use_units),striatum_sua_grp.domain_idx(use_units)]);

unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
    [plot_day_grp(use_units),striatum_sua_grp.domain_idx(use_units)]);

figure;
h = tiledlayout(n_domains,1);
stim_colormap = ap.colormap('BKR',3);
binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');
for curr_domain = 1:n_domains
    nexttile; hold on;
    set(gca,'ColorOrder',stim_colormap);

    curr_data_idx = unit_responsive_mean_group(:,2) == curr_domain;
    errorbar(binned_days_x(unit_responsive_mean_group(curr_data_idx,1)), ...
        unit_responsive_mean(curr_data_idx,:), ...
        unit_responsive_sem(curr_data_idx,:),'linewidth',2);
    ylabel('Frac. responsive units');
    xline(0);
end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Fig 3X] Cortex map passive max

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,2,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

stim_t = wf_t > 0 & wf_t < 0.2;
[wf_avg,wf_avg_grp] = ap.groupfun(@mean,cell2mat(wf.V_event_align), ...
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


%% [Fig 3X] Cortex passive ROI average

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,2,Inf];
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
h = tiledlayout(n_domains,max(plot_day_grp)*length(unique_stim),'TileSpacing','tight');
for curr_domain = 1:size(striatum_wf_roi,3)
    for curr_ld = unique(plot_day_grp)'
        for curr_stim = unique_stim'
            nexttile; axis off;
            curr_data_idx = wf_striatum_roi_avg_grp(:,1) == curr_ld & ...
                wf_striatum_roi_avg_grp(:,2) == curr_stim;

            ap.errorfill(wf_t,wf_striatum_roi_avg(curr_data_idx,:,curr_domain), ...
                wf_striatum_roi_sem(curr_data_idx,:,curr_domain), ...
                day_colormap(curr_ld,:));
        end
    end
end
linkaxes(h.Children,'xy');
title(h,'Cortex');
ap.prettyfig;


%% [Fig 3X] Cortex + striatum passive ROI max

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-Inf,-2,0,2,Inf];
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
striatum_plot_days_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

% Average day, get max in time window
stim_t = [0,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

% (cortex)
[wf_striatum_roi_dayavg,wf_striatum_roi_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi, ...
    (1:size(wf_striatum_roi,1))', ...
    [wf_grp.animal,cortex_plot_day_grp,cell2mat(wf.trial_stim_values)]);

wf_striatum_roi_max = max(wf_striatum_roi_dayavg(:,cortex_stim_t,:),[],2);

[wf_striatum_roi_max_avg,wf_striatum_roi_max_avg_grp] = ...
    ap.nestgroupfun({@mean,@mean},wf_striatum_roi_max, ...
    wf_striatum_roi_grp(:,1),wf_striatum_roi_grp(:,2:3));

wf_striatum_roi_max_sem = ...
    ap.nestgroupfun({@mean,@AP_sem},wf_striatum_roi_max, ...
    wf_striatum_roi_grp(:,1),wf_striatum_roi_grp(:,2:3));

% (striatum)
[striatum_mua_dayavg,striatum_mua_dayavg_grp] = ...
    ap.groupfun(@mean,striatum_mua, ...
    [striatum_mua_grp.animal,striatum_plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);

striatum_mua_dayavg_tmax = max(striatum_mua_dayavg(:,striatum_stim_t),[],2);

[striatum_mua_max_avg,striatum_mua_max_grp] = ap.nestgroupfun({@mean,@mean}, ...
    striatum_mua_dayavg_tmax,striatum_mua_dayavg_grp(:,1), ...
    striatum_mua_dayavg_grp(:,2:end));
striatum_mua_max_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    striatum_mua_dayavg_tmax,striatum_mua_dayavg_grp(:,1), ...
    striatum_mua_dayavg_grp(:,2:end));

% Plot activity by day
figure;
h = tiledlayout(n_domains,2,'TileSpacing','tight');
stim_colormap = ap.colormap('BKR',3);
binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');
for curr_domain = 1:size(striatum_wf_roi,3)
    % (cortex)
    nexttile; hold on;
    set(gca,'ColorOrder',stim_colormap);
    errorbar(binned_days_x, ...
        reshape(wf_striatum_roi_max_avg(:,:,curr_domain),[],length(plot_day_bins)-1)', ...
        reshape(wf_striatum_roi_max_sem(:,:,curr_domain),[],length(plot_day_bins)-1)', ...
        'linewidth',2);
    ylabel('Max \DeltaF/F_0');
    axis padded
    xline(0);

    % (striatum)
    nexttile; hold on; 
    set(gca,'ColorOrder',stim_colormap);
    curr_str_idx = striatum_mua_max_grp(:,3) == curr_domain;
    errorbar(binned_days_x, ...
        reshape(striatum_mua_max_avg(curr_str_idx),[],length(plot_day_bins)-1)', ...
        reshape(striatum_mua_max_sem(curr_str_idx),[],length(plot_day_bins)-1)', ...
        'linewidth',2);
    ylabel('Max \DeltaR/R_0');
    axis padded
    xline(0);
end
linkaxes(h.Children(1:2:end),'xy');
linkaxes(h.Children(2:2:end),'xy');

ap.prettyfig;

% ~stats~
fprintf('---- STATS ----\n')
for curr_compare_day = 1:length(plot_day_bins)-2

    compare_day_grps = curr_compare_day+[0,1];

    cortex_stat_usedata = ismember(wf_striatum_roi_grp(:,2),compare_day_grps);
    [cortex_stat_meas,cortex_stat_grp] = ap.nestgroupfun({@mean,@diff},wf_striatum_roi_max(cortex_stat_usedata,:,:), ...
        wf_striatum_roi_grp(cortex_stat_usedata,2),wf_striatum_roi_grp(cortex_stat_usedata,3));

    striatum_stat_usedata = ismember(striatum_mua_dayavg_grp(:,2),compare_day_grps);
    [striatum_stat_meas,striatum_stat_grp] = ap.nestgroupfun({@mean,@diff},striatum_mua_dayavg_tmax(striatum_stat_usedata), ...
        striatum_mua_dayavg_grp(striatum_stat_usedata,2),striatum_mua_dayavg_grp(striatum_stat_usedata,3:end));

    [~,~,cortex_shuff_grp] = unique(wf_striatum_roi_grp(cortex_stat_usedata,[1,3]),'rows');
    [~,~,striatum_shuff_grp] = unique(striatum_mua_dayavg_grp(striatum_stat_usedata,[1,3,4]),'rows');

    n_shuff = 10000;
    cortex_stat_null = nan(size(cortex_stat_meas,1),n_shuff,n_domains);
    striatum_stat_null = nan(length(striatum_stat_meas),n_shuff);
    for curr_shuff = 1:n_shuff
        cortex_data_shuff = ap.shake(wf_striatum_roi_max(cortex_stat_usedata,:,:),1,cortex_shuff_grp);
        cortex_stat_null(:,curr_shuff,:) = ap.nestgroupfun({@mean,@diff},cortex_data_shuff, ...
            wf_striatum_roi_grp(cortex_stat_usedata,2),wf_striatum_roi_grp(cortex_stat_usedata,3));

        striatum_data_shuff = ap.shake(striatum_mua_dayavg_tmax(striatum_stat_usedata),1,striatum_shuff_grp);
        striatum_stat_null(:,curr_shuff) = ap.nestgroupfun({@mean,@diff},striatum_data_shuff, ...
            striatum_mua_dayavg_grp(striatum_stat_usedata,2),striatum_mua_dayavg_grp(striatum_stat_usedata,3:end));
    end

    cortex_stat_rank = permute(tiedrank(permute([cortex_stat_meas,cortex_stat_null],[2,1,3])),[2,1,3]);
    cortex_stat_p = 1-cortex_stat_rank(:,1,:)/(n_shuff+1);

    striatum_stat_rank = tiedrank([striatum_stat_meas,striatum_stat_null]')';
    striatum_stat_p = 1-striatum_stat_rank(:,1)/(n_shuff+1);

    fprintf('Cortex: day grps %d vs %d\n',compare_day_grps);
    stat_sig = discretize(cortex_stat_p < 0.05,[0,1,Inf],["","*"]);
    for curr_domain = 1:n_domains
        for curr_stim = unique(cortex_stat_grp(:,1))'
            curr_stat_idx = ismember(cortex_stat_grp,curr_stim,'rows');
            fprintf('ROI%d, Stim %3.f, p = %.2g%s\n', ...
                curr_domain,cortex_stat_grp(curr_stat_idx),cortex_stat_p(curr_stat_idx,1,curr_domain),stat_sig(curr_stat_idx,1,curr_domain));
        end
    end

    fprintf('Striatum: day grps %d vs %d\n',compare_day_grps);
    stat_sig = discretize(striatum_stat_p < 0.05,[0,1,Inf],["","*"]);
    for curr_domain = 1:n_domains
        for curr_stim = unique(striatum_stat_grp(:,1))'
            curr_stat_idx = ismember(striatum_stat_grp,[curr_stim,curr_domain],'rows');
            fprintf('D%d, Stim %3.f, p = %.2g%s\n', ...
                curr_domain,striatum_stat_grp(curr_stat_idx,1),striatum_stat_p(curr_stat_idx),stat_sig(curr_stat_idx));
        end
    end

end


%% [Fig 4X] Striatum passive unit responses

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

% Set striatal domains to plot (combines if multiple)
plot_domains = 1:2;

% Set days to group
plot_day_bins = [-inf,-2,0,Inf];
plot_day_grp = discretize(max(-inf,striatum_sua_grp.ld),plot_day_bins);

% Get max activity to stim onset
stim_t = [0.05,0.15];
striatum_sua_tavg = squeeze(mean(striatum_sua(:,isbetween(psth_t,stim_t(1),stim_t(2)),:),2));

n_stim = size(striatum_sua,3);

% Plot heatmap
[~,max_stim] = max(striatum_sua_tavg,[],2);
for curr_celltype = celltypes
    figure;
    colormap(ap.colormap('BWR',[],2));
    h = tiledlayout(max(plot_day_grp)*n_stim,n_stim,'TileSpacing','none');
    title(h,curr_celltype);
    for curr_day_grp = 1:length(plot_day_bins)-1
        for curr_max_stim = 1:n_stim

            % Get units and sorting
            curr_units = find( ...
                ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
                plot_day_grp == curr_day_grp & ...
                max_stim == curr_max_stim & ...
                striatum_sua_grp.(curr_celltype) & ...
                any(striatum_sua,[2,3]));
     
            [~,sort_idx] = sort(max(striatum_sua_tavg(curr_units,curr_max_stim),[],2),'descend');
            plot_units = curr_units(sort_idx);

            % (get rec/IDs of cells for single-unit PSTH plotting)
            curr_sorted_unit_coordinate = ...
                [ephys.animal(striatum_sua_grp.rec(plot_units)), ...
                ephys.rec_day(striatum_sua_grp.rec(plot_units)), ...
                num2cell(striatum_sua_grp.unit_id(plot_units))];

            % (smooth relative to number of number of units);
            max_n_cells = max(ap.groupfun(@sum, ...
                +(ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
                striatum_sua_grp.(curr_celltype)),plot_day_grp));
            smooth_n = 5*max_n_cells/500;

            for curr_stim = 1:n_stim
                nexttile;
                imagesc(psth_t,[],movmean(striatum_sua(plot_units,:,curr_stim),smooth_n,1));
                clim([-1,1])
                xlim([-0.2,0.8])
                axis off;
            end
        end
    end
    linkaxes(h.Children,'xy');
    ap.prettyfig;
end

% Plot max by stim response
[~,max_stim] = max(striatum_sua_tavg,[],2);
for curr_celltype = celltypes
    figure;
    h = tiledlayout(n_stim*max(plot_day_grp),1,'TileSpacing','compact');
    title(h,curr_celltype);

    stim_colormap = ap.colormap('BKR',n_stim);
    for curr_day_grp = 1:length(plot_day_bins)-1
        for curr_max_stim = 1:3

            curr_units = ...
                ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
                plot_day_grp == curr_day_grp & striatum_sua_grp.(curr_celltype) & ...
                max_stim == curr_max_stim;

            curr_act_mean = ap.nestgroupfun({@nanmean,@nanmean},striatum_sua_tavg(curr_units,:), ...
                striatum_sua_grp.animal(curr_units),plot_day_grp(curr_units));

            curr_act_sem = ap.nestgroupfun({@nanmean,@AP_sem},striatum_sua_tavg(curr_units,:), ...
                striatum_sua_grp.animal(curr_units),plot_day_grp(curr_units));

            nexttile; hold on;
            b = bar(curr_act_mean,'FaceColor','flat','CData',stim_colormap);
            errorbar(curr_act_mean,curr_act_sem, ...
                'marker','none','linestyle','none','color','k','linewidth',1);

        end
    end
    set(h.Children,'YTickLabel','')
    set(h.Children(2:end),'XTickLabel','')
    linkaxes(h.Children,'xy');
    ap.prettyfig;
end

% Average stim response grouped by celltype
plot_celltypes = ["tan","fsi","msn"];
celltype_id = sum(cell2mat(arrayfun(@(x) ...
    striatum_sua_grp.(plot_celltypes(x)).*x, ...
    1:length(plot_celltypes),'uni',false)),2);

stim_colormap = ap.colormap('BWR',n_stim);

figure;
h = tiledlayout(1,length(celltypes),'TileSpacing','compact');
for curr_celltype = celltypes

    curr_units =  ...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        striatum_sua_grp.(curr_celltype);

    curr_act_mean = ap.nestgroupfun({@nanmean,@nanmean},striatum_sua_tavg(curr_units,:), ...
        striatum_sua_grp.animal(curr_units),plot_day_grp(curr_units));

    curr_act_sem = ap.nestgroupfun({@nanmean,@AP_sem},striatum_sua_tavg(curr_units,:), ...
        striatum_sua_grp.animal(curr_units),plot_day_grp(curr_units));

    nexttile; hold on;
    b = bar(curr_act_mean','FaceColor','flat');
    [b.CData] = deal(stim_colormap);
    errorbar(vertcat(b.XEndPoints)',curr_act_mean',curr_act_sem', ...
        'marker','none','linestyle','none','color','k','linewidth',1);
    title(curr_celltype);

    % ~stats~
    fprintf('---- STATS ----\n')
    compare_day_grps = [1,2];
    stat_meas = diff(curr_act_mean(compare_day_grps,:),[],1);

    curr_stat_units = curr_units & ismember(plot_day_grp,compare_day_grps);
    n_shuff = 1000;
    stat_null = nan(n_shuff,3);
    for curr_shuff = 1:n_shuff
        plot_day_grp_shuff = ap.shake(plot_day_grp(curr_stat_units),1,striatum_sua_grp.animal(curr_stat_units));
        curr_act_mean_shuff = ap.nestgroupfun({@nanmean,@nanmean},striatum_sua_tavg(curr_stat_units,:), ...
            striatum_sua_grp.animal(curr_stat_units),plot_day_grp_shuff);
        stat_null(curr_shuff,:) = diff(curr_act_mean_shuff,[],1);
    end
    stat_rank = tiedrank([stat_meas;stat_null]);
    stat_p = 1-stat_rank(1,:)/(n_shuff+1);
    fprintf('%s day grps %d,%d p = %.2g L, %.2g C, %.2g R\n',curr_celltype,compare_day_grps,stat_p);

end
linkaxes(h.Children,'xy');
ap.prettyfig;


%% [Supp. Fig 1x] Striatal domain clustering and classification

%%% Load data for figure
load_dataset = 'noact';
AP_longstriatum_load_data;
%%%

% Choose animal and day to plot
use_animal = 'AM026';
use_ld = 0;
use_rec = strcmp(bhv.animal,use_animal) & bhv.days_from_learning == use_ld;
use_cortex_kernel = ctx_str_maps.cortex_striatum_map{use_rec};

domain_color = {'R','G','B'};

% Plot all domains
% (grayscale and colored by domain)
figure;
h = tiledlayout(size(use_cortex_kernel,3),2,'TileSpacing','none');
for curr_depth=1:size(use_cortex_kernel, 3)
    nexttile;
    imagesc(use_cortex_kernel(:,:,curr_depth));
    axis image off;
    ap.wf_draw('cortex',[0.5,0.5,0.5]);
    colormap(gca,ap.colormap('WK',[],2));

    nexttile;
    imagesc(use_cortex_kernel(:,:,curr_depth));
    axis image off
    ap.wf_draw('cortex',[0.5,0.5,0.5]);
    colormap(gca,ap.colormap(['W' domain_color{domain_idx_rec{use_rec}(curr_depth)}],[],2));
end
clim(h.Children,[0,0.01]);
ap.prettyfig;

% Plot domain means
figure;
h = tiledlayout(n_domains,1,'tilespacing','none');
for curr_domain = 1:n_domains
    nexttile;
    imagesc(kmeans_cluster_mean(:,:,curr_domain));
    axis image off
    colormap(gca,ap.colormap(['W',domain_color{curr_domain}],[],2));
    ap.wf_draw('cortex',[0.5,0.5,0.5]);
end
clim(h.Children,[0,0.01]);
ap.prettyfig;

% Load data from example recording
animal = bhv.animal{use_rec};
rec_day = bhv.rec_day{use_rec};
recordings = plab.find_recordings(animal,rec_day,'stim_wheel*');
rec_time = recordings.recording{end};
load_parts.ephys = true;
ap.load_recording;

% Plot units and overlay clustering
figure; 
ax = axes; hold on;

domain_color_rgb = [1,0,0;0,1,0;0,0,1];
domain_im = permute(domain_color_rgb(domain_idx_rec{use_rec},:),[1,3,2]);
imagesc(ax,[],ctx_str_maps.depth_group_edges{use_rec},domain_im);
ax.YDir = 'reverse';

ap.plot_unit_depthrate(spike_times_timelite,spike_templates,template_depths,[],ax);
yline(ctx_str_maps.depth_group_edges{use_rec},'linewidth',2,'color',[0.5,0.5,0.5]);

ap.prettyfig;

% Plot clustered MUA
depth_group = discretize(spike_depths,ctx_str_maps.depth_group_edges{use_rec});

plot_t = [100,140];
bin_t = 0.1;
domain_mua = zeros(n_domains,diff(plot_t)/bin_t);
for curr_domain = 1:n_domains
    curr_spikes = spike_times_timelite(isbetween(spike_times_timelite,plot_t(1),plot_t(2)) & ...
        ismember(depth_group,find(domain_idx_rec{use_rec} == curr_domain)));
    domain_mua(curr_domain,:) = histcounts(curr_spikes,plot_t(1):bin_t:plot_t(2))./bin_t;
end

figure;
plot(conv(plot_t(1):bin_t:plot_t(2),[0.5,0.5],'valid'),domain_mua','linewidth',2);
set(gca,'ColorOrder',domain_color_rgb);
axis off;
ap.scalebar(10,500);
ap.prettyfig


%% (end timer)

toc







