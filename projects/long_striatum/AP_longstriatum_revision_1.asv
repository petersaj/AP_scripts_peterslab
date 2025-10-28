%% Activity by learning for individual animals
% (also fold in task vs passive?)

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

plot_day_bins = [-Inf,-3,-2,-1,0,1,Inf];
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
figure('Name','Fig 2 psth'); h = tiledlayout(n_domains*2,1,'TileSpacing','tight');
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


%%%%% Get t-max point
stim_t = [0,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

[wf_striatum_roi_nomove_animal,wf_striatum_roi_nomove_animal_grp] = ...
    ap.groupfun(@nanmean,wf_striatum_roi_nomove, ...
    [wf_grp.animal,cortex_plot_day_grp]);
wf_striatum_roi_nomove_animal_tmax = max(wf_striatum_roi_nomove_animal(:,cortex_stim_t,:),[],2);

[striatum_mua_nomove_animal,striatum_mua_nomove_animal_grp] = ...
    ap.groupfun(@nanmean,striatum_mua_nomove, ...
    [striatum_mua_grp.animal,striatum_plot_day_grp,striatum_mua_grp.domain_idx]);
[~,~,striatum_shuffgroup] = unique(striatum_mua_nomove_animal_grp(:,[1,3]),'rows');
striatum_mua_nomove_animal_tmax = max(striatum_mua_nomove_animal(:,striatum_stim_t),[],2);


%%%% repackage into grid - dirty code
wf_striatum_roi_nomove_animal_tmax_grid = nan(max(wf_grp.animal),length(plot_day_bins)-1,3);
for curr_roi = 1:n_domains
    wf_striatum_roi_nomove_animal_tmax_grid(sub2ind(size(wf_striatum_roi_nomove_animal_tmax_grid), ...
        wf_striatum_roi_nomove_animal_grp(:,1),wf_striatum_roi_nomove_animal_grp(:,2), ...
        repelem(curr_roi,size(wf_striatum_roi_nomove_animal_grp,1),1))) = ...
        wf_striatum_roi_nomove_animal_tmax(:,:,curr_roi);
end

striatum_mua_nomove_animal_tmax_grid = nan(max(striatum_mua_grp.animal),length(plot_day_bins)-1,3);
striatum_mua_nomove_animal_tmax_grid(sub2ind(size(striatum_mua_nomove_animal_tmax_grid), ...
    striatum_mua_nomove_animal_grp(:,1),striatum_mua_nomove_animal_grp(:,2),striatum_mua_nomove_animal_grp(:,3))) = ...
    striatum_mua_nomove_animal_tmax;

figure; hold on;
for curr_daygrp = 1:length(plot_day_bins)-1
    plot(reshape(striatum_mua_nomove_animal_tmax_grid(:,curr_daygrp,1),[],1), ...
        reshape(wf_striatum_roi_nomove_animal_tmax_grid(:,curr_daygrp,2),[],1),'.','MarkerSize',10)
    plot(nanmean(reshape(striatum_mua_nomove_animal_tmax_grid(:,curr_daygrp,1),[],1)), ...
        nanmean(reshape(wf_striatum_roi_nomove_animal_tmax_grid(:,curr_daygrp,2),[],1)),'.','MarkerSize',30)
end
xlabel('Vis-Str');
ylabel('mPFC');

% normalize day 0 = 1
wf_striatum_roi_nomove_animal_tmax_grid_norm = ...
    wf_striatum_roi_nomove_animal_tmax_grid./ ...
    wf_striatum_roi_nomove_animal_tmax_grid(:,plot_day_bins(1:end-1)==0,:);

striatum_mua_nomove_animal_tmax_grid_norm = ...
    striatum_mua_nomove_animal_tmax_grid./ ...
    striatum_mua_nomove_animal_tmax_grid(:,plot_day_bins(1:end-1)==0,:);

figure; hold on;
for curr_daygrp = 1:length(plot_day_bins)-1
    plot(reshape(striatum_mua_nomove_animal_tmax_grid(:,curr_daygrp,1),[],1), ...
        reshape(wf_striatum_roi_nomove_animal_tmax_grid_norm(:,curr_daygrp,2),[],1),'.','MarkerSize',10)
end
xlabel('Vis-Str');
ylabel('mPFC');

%% [Fig 3B-I] Passive Cortex + striatum PSTH

%%% Load data for figure
load_dataset = 'passive';
Marica_2025.figures.load_data;
%%%

plot_day_bins = [-Inf,-2,-1,0,1,Inf];
striatum_plot_days_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

[striatum_mua_avg,striatum_mua_avg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [striatum_plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);
striatum_mua_sem = ap.nestgroupfun({@nanmean,@AP_sem}, ...
    striatum_mua,striatum_mua_grp.animal, ...
    [striatum_plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);

[wf_striatum_roi_avg,wf_striatum_roi_avg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},wf_striatum_roi, ...
    wf_grp.animal,[cortex_plot_day_grp,cell2mat(wf.trial_stim_values)]);
wf_striatum_roi_sem = ...
    ap.nestgroupfun({@nanmean,@AP_sem},wf_striatum_roi, ...
    wf_grp.animal,[cortex_plot_day_grp,cell2mat(wf.trial_stim_values)]);

% unique_stim = unique(striatum_mua_grp.stim);
unique_stim = [90,0];
stim_color = {'KR','KW'};

figure('Name','Fig 3 psth');
h = tiledlayout(n_domains*2,max(striatum_plot_days_grp)*length(unique_stim), ...
    'TileIndexing','ColumnMajor','TileSpacing','tight');
for curr_stim = unique_stim

    [~,curr_stim_idx] = ismember(curr_stim,unique_stim);
    day_colormap = ap.colormap(stim_color{curr_stim_idx},length(plot_day_bins)-1);

    for curr_day = 1:length(plot_day_bins)-1
        for curr_domain = 1:n_domains

            % (cortex)
            nexttile; axis off;
            curr_data_idx = wf_striatum_roi_avg_grp(:,1) == curr_day & ...
                wf_striatum_roi_avg_grp(:,2) == curr_stim;
            ap.errorfill(wf_t,wf_striatum_roi_avg(curr_data_idx,:,curr_domain), ...
                wf_striatum_roi_sem(curr_data_idx,:,curr_domain), ...
                day_colormap(curr_day,:));

            % (striatum)
            nexttile; axis off;
            curr_data_idx = striatum_mua_avg_grp(:,1) == curr_day & ...
                striatum_mua_avg_grp(:,2) == curr_stim & ...
                striatum_mua_avg_grp(:,3) == curr_domain;
            ap.errorfill(psth_t,striatum_mua_avg(curr_data_idx,:),striatum_mua_sem(curr_data_idx,:), ...
                day_colormap(curr_day,:));

        end
    end
end
linkaxes(h.Children(1:2:end),'y');
linkaxes(h.Children(2:2:end),'y');
xlim(h.Children,[0,0.5]);
ap.prettyfig;


%% [Fig 3J-K] Passive cortex + striatum max

%%% Load data for figure
load_dataset = 'passive';
Marica_2025.figures.load_data;
%%%

plot_day_bins = [-Inf,-2,-1,0,1,Inf];
% plot_day_bins = [-Inf,-2,0,Inf];

cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
striatum_plot_days_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

% Average day, get max in time window
stim_t = [0,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

% (cortex)
[wf_striatum_roi_dayavg,wf_striatum_roi_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},wf_striatum_roi, ...
    (1:size(wf_striatum_roi,1))', ...
    [wf_grp.animal,cortex_plot_day_grp,cell2mat(wf.trial_stim_values)]);

wf_striatum_roi_max = max(wf_striatum_roi_dayavg(:,cortex_stim_t,:),[],2);

[wf_striatum_roi_max_avg,wf_striatum_roi_max_avg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},wf_striatum_roi_max, ...
    wf_striatum_roi_grp(:,1),wf_striatum_roi_grp(:,2:3));

wf_striatum_roi_max_sem = ...
    ap.nestgroupfun({@nanmean,@AP_sem},wf_striatum_roi_max, ...
    wf_striatum_roi_grp(:,1),wf_striatum_roi_grp(:,2:3));

% (striatum)
[striatum_mua_dayavg,striatum_mua_dayavg_grp] = ...
    ap.groupfun(@nanmean,striatum_mua, ...
    [striatum_mua_grp.animal,striatum_plot_days_grp,striatum_mua_grp.stim,striatum_mua_grp.domain_idx]);

striatum_mua_dayavg_tmax = max(striatum_mua_dayavg(:,striatum_stim_t),[],2);

[striatum_mua_max_avg,striatum_mua_max_grp] = ap.nestgroupfun({@nanmean,@nanmean}, ...
    striatum_mua_dayavg_tmax,striatum_mua_dayavg_grp(:,1), ...
    striatum_mua_dayavg_grp(:,2:end));
striatum_mua_max_sem = ap.nestgroupfun({@nanmean,@AP_sem}, ...
    striatum_mua_dayavg_tmax,striatum_mua_dayavg_grp(:,1), ...
    striatum_mua_dayavg_grp(:,2:end));

% Plot activity by day
figure('Name','Fig 3 tmax');
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


% trying dstr vs dmpfc

use_stim = 90;

max_dim = [length(unique(bhv.animal)),length(plot_day_bins)-1];

s = accumarray(striatum_mua_dayavg_grp(striatum_mua_dayavg_grp(:,3)==use_stim,[1,2,4]), ...
    striatum_mua_dayavg_tmax(striatum_mua_dayavg_grp(:,3)==use_stim),[max_dim,3],[],NaN);

c = accumarray(wf_striatum_roi_grp(wf_striatum_roi_grp(:,3)==use_stim,[1,2]), ...
    double(wf_striatum_roi_max(wf_striatum_roi_grp(:,3)==use_stim,:,2)),max_dim,[],NaN);

figure; tiledlayout;
nexttile; hold on;
plot(s(:,:,1)',c(:,:,1)');
plot(nanmean(s(:,:,1)),nanmean(c(:,:,1)),'k','linewidth',2);

sn = s./s(:,end,:);
cn = c./c(:,end,:);
nexttile; hold on;
plot(sn(:,:,1)',cn(:,:,1)')
plot(nanmean(sn(:,:,1)),nanmean(cn(:,:,1)),'k','linewidth',2)


%% Task cortex + striatum max

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

plot_day_bins = [-3,-2,-1,0,1];
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

rxn_cutoff = 0.3;

% Average day, get max in time window
stim_t = [0,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

% (cortex)
use_wf_trials = wf_grp.rxn > rxn_cutoff;
[wf_striatum_roi_dayavg,wf_striatum_roi_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},wf_striatum_roi(use_wf_trials,:,1), ...
    (1:sum(use_wf_trials,1))', ...
    [wf_grp.animal(use_wf_trials),cortex_plot_day_grp(use_wf_trials)]);

wf_striatum_roi_max = max(wf_striatum_roi_dayavg(:,cortex_stim_t,:),[],2);

% (striatum)
use_striatum_trials = striatum_mua_grp.rxn > rxn_cutoff;
[striatum_mua_dayavg,striatum_mua_dayavg_grp] = ...
    ap.groupfun(@nanmean,striatum_mua(use_striatum_trials,:,:), ...
    [striatum_mua_grp.animal(use_striatum_trials),striatum_plot_day_grp(use_striatum_trials),striatum_mua_grp.domain_idx(use_striatum_trials)]);

striatum_mua_dayavg_tmax = max(striatum_mua_dayavg(:,striatum_stim_t),[],2);

% reshape into grid
max_dim = [length(unique(bhv.animal)),length(plot_day_bins)-1];

c = accumarray(wf_striatum_roi_grp, ...
    double(wf_striatum_roi_max),max_dim,[],NaN);

s = accumarray(striatum_mua_dayavg_grp, ...
    striatum_mua_dayavg_tmax,[max_dim,3],[],NaN);

figure; tiledlayout;
nexttile; hold on;
plot(s(:,:,1)',c(:,:,1)');
plot(nanmean(s(:,:,1)),nanmean(c(:,:,1)),'k','linewidth',2);

sn = s./s(:,end,:);
cn = c./c(:,end,:);
nexttile; hold on;
plot(sn(:,:,1)',cn(:,:,1)')
plot(nanmean(sn(:,:,1)),nanmean(cn(:,:,1)),'k','linewidth',2)


% try averaging trial-trial trends normalized to training time
figure; hold on;
n_points = 50;
data_interp = nan(2,n_points,14);
for use_animal = 1:14

    s = mean(striatum_mua(striatum_mua_grp.animal==use_animal & striatum_mua_grp.domain_idx==1,isbetween(psth_t,0.05,0.15),1),2);

    s_days = unique(striatum_mua_grp.ld(striatum_mua_grp.animal==use_animal & striatum_mua_grp.domain_idx==1));
    c = mean(wf_striatum_roi(wf_grp.animal==use_animal & ismember(wf_grp.ld,s_days),isbetween(wf_t,0.05,0.15),1),2);

    if isempty(c)
        continue
    end

    s_smooth = smoothdata(s,'movmean',100);
    c_smooth = smoothdata(c,'movmean',100);
    plot(smoothdata(s,'movmean',100),smoothdata(c,'movmean',100));

    data_interp(1,:,use_animal) = s_smooth(round(linspace(1,length(s_smooth),n_points)));
    data_interp(2,:,use_animal) = c_smooth(round(linspace(1,length(c_smooth),n_points)));
end

figure; hold on;
plot(permute(data_interp(1,:,:),[2,3,1]),permute(data_interp(2,:,:),[2,3,1]));
plot(nanmean(data_interp(1,:,:),3),nanmean(data_interp(2,:,:),3),'k','linewidth',2)


%% Cortex movie passive

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-2,-1,0,0];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

wf_V_event_align_cat = cell2mat(wf.V_event_align);
wf_trial_stim_values_cat = cell2mat(wf.trial_stim_values);


% use_trials = wf_grp.animal == 2;
% use_trials = true(size(wf_grp.animal));
% use_trials = wf_trial_stim_values_cat==90;

n_learned_day = cellfun(@(x) max([0, ...
    find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'));
use_trials = wf_trial_stim_values_cat==90 & n_learned_day(wf_grp.animal) > 3;


[wf_avg,wf_avg_grp] = ap.nestgroupfun({@nanmean,@nanmean},wf_V_event_align_cat(use_trials,:,:), ...
    wf_grp.animal(use_trials),[plot_day_grp(use_trials),wf_trial_stim_values_cat(use_trials)]);

% [wf_avg,wf_avg_grp] = ap.groupfun(@nanmean,wf_V_event_align_cat(use_trials,:,:), ...
%     [wf_grp.animal(use_trials),plot_day_grp(use_trials),wf_trial_stim_values_cat(use_trials)]);

baseline_t = [-0.2,-0.1];
baseline_t_idx = isbetween(wf_t,baseline_t(1),baseline_t(2));
wf_px = plab.wf.svd2px(U_master,permute(wf_avg-nanmean(wf_avg(:,baseline_t_idx,:),2),[3,2,1]));

ap.imscroll(wf_px,wf_t);
colormap(ap.colormap('PWG',[],1.5));
clim([-1,1]*3e-3);
axis image off;
ap.wf_draw('ccf',[0.5,0.5,0.5]);


%% Cortex movie task

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-2,-1,0,0];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

move_leeway = 0.1; % time pre-movement to exclude
wf_V_event_align_cat = cell2mat(cellfun(@(x) x(:,:,:,1), ...
    wf.V_event_align,'uni',false)).*ap.nanout(wf_t > wf_grp.rxn-move_leeway);

use_trials = true(size(wf_grp.animal));

% n_learned_day = cellfun(@(x) max([0, ...
%     find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
%     unique(bhv.animal,'stable'));
% use_trials = n_learned_day(wf_grp.animal)>3;


[wf_avg,wf_avg_grp] = ap.groupfun(@nanmean, ...
    wf_V_event_align_cat(use_trials,:,:), ...
    [wf_grp.animal(use_trials),plot_day_grp(use_trials)]);

% [wf_avg,wf_avg_grp] = ap.nestgroupfun({@nanmean,@nanmean}, ...
%     wf_V_event_align_cat(use_trials,:,:), ...
%     wf_grp.animal(use_trials),[plot_day_grp(use_trials)]);

baseline_t = [-0.2,-0.1];
baseline_t_idx = isbetween(wf_t,baseline_t(1),baseline_t(2));
wf_px = plab.wf.svd2px(U_master,permute(wf_avg-nanmean(wf_avg(:,baseline_t_idx,:),2),[3,2,1]));

ap.imscroll(wf_px,wf_t);
colormap(ap.colormap('PWG',[],1.5));
clim([-1,1]*3e-3);
axis image off;
ap.wf_draw('ccf',[0.5,0.5,0.5]);


%% mPFC TASK: get day it's different from first

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

rxn_cutoff = 0.3;

% Average day, get max in time window
stim_t = [0.1,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

mpfc_avg = max(wf_striatum_roi(:,cortex_stim_t,2,1),[],2);

td_cat = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));

wf_grp.td = cell2mat(cellfun(@(animal,data) ...
    repmat(animal,size(data,1),1), ...
    num2cell(td_cat),wf.V_event_align,'uni',false));

mpfc_diffday = cell(length(unique(bhv.animal)),1);
for curr_animal = 1:length(unique(bhv.animal,'stable'))

    curr_animal_tds = unique(wf_grp.td(wf_grp.animal == curr_animal));

    curr_baseline_trials = wf_grp.animal == curr_animal & ...
        wf_grp.td == min(curr_animal_tds) & ...
        wf_grp.rxn >= rxn_cutoff;

    n_shuff = 1000;
    mpfc_p = nan(size(curr_animal_tds));
    for curr_day = curr_animal_tds(2:end)'

        curr_day_trials = wf_grp.animal == curr_animal & ...
            wf_grp.td == curr_day & ...
            wf_grp.rxn >= rxn_cutoff;

        if sum(curr_day_trials) < 10
            continue
        end

        curr_compare_idx = curr_baseline_trials + curr_day_trials.*2;
        curr_mpfc_diff = diff(ap.groupfun(@mean,mpfc_avg,curr_compare_idx.*ap.nanout(curr_compare_idx==0)));

        curr_compare_idx_shuff = ap.shake(repmat(curr_compare_idx,1,n_shuff),1,curr_compare_idx>0);
        curr_mpfc_diff_shuff = ...
            (mpfc_avg'*(curr_compare_idx_shuff==2))./sum(curr_compare_idx_shuff==2) - ...
            (mpfc_avg'*(curr_compare_idx_shuff==1))./sum(curr_compare_idx_shuff==1);

        stat_rank = tiedrank([curr_mpfc_diff,curr_mpfc_diff_shuff]');
        stat_p = 1-stat_rank(1)/(n_shuff+1);

        mpfc_p(curr_animal_tds == curr_day) = stat_p;
    end

    mpfc_diffday{curr_animal} = (1:length(curr_animal_tds))' - min([find(mpfc_p < 0.05,1),nan]);
    disp(curr_animal);

end

n_learned_day = cellfun(@(x) max([0, ...
    find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'));

disp([n_learned_day,cellfun(@(x) min([find(x==0),NaN]),mpfc_diffday)]);

% Sub LD for mPFC day
wf_grp.ld = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
        num2cell(cell2mat(mpfc_diffday)),wf.V_event_align,'uni',false));
striatum_mua_grp.ld = cell2mat(cellfun(@(grp,n_trials,n_mua) repmat(grp,n_trials*n_mua,1), ...
        num2cell(cell2mat(mpfc_diffday)),num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

% % (to put back)
% wf_grp.ld = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
%         num2cell(bhv.days_from_learning),wf.V_event_align,'uni',false));
% striatum_mua_grp.ld = cell2mat(cellfun(@(grp,n_trials,n_mua) repmat(grp,n_trials*n_mua,1), ...
%         num2cell(bhv.days_from_learning),num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));


%% mPFC PASSIVE: get day it's different from first

%%% Load data for figure
load_dataset = 'passive';
Marica_2025.figures.load_data;
%%%

% Average day, get max in time window
stim_t = [0,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

mpfc_avg = max(wf_striatum_roi(:,cortex_stim_t,2,1),[],2);

td_cat = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));

wf_grp.td = cell2mat(cellfun(@(animal,data) ...
    repmat(animal,size(data,1),1), ...
    num2cell(td_cat),wf.V_event_align,'uni',false));

mpfc_diffday = cell(length(unique(bhv.animal)),1);
for curr_animal = 1:length(unique(bhv.animal,'stable'))

    curr_animal_tds = unique(wf_grp.td(wf_grp.animal == curr_animal));

    curr_baseline_trials = ...
        wf_grp.animal == curr_animal & ...
        ismember(wf_grp.td, [1]) & ...
        wf_grp.stim == 90;

    n_shuff = 1000;
    mpfc_p = nan(size(curr_animal_tds));
    for curr_day = curr_animal_tds(2:end)'

        curr_day_trials = wf_grp.animal == curr_animal & ...
            wf_grp.td == curr_day & ...
            wf_grp.stim == 90;

        if sum(curr_day_trials) < 10
            continue
        end

        curr_compare_idx = curr_baseline_trials + curr_day_trials.*2;
        curr_mpfc_diff = diff(ap.groupfun(@mean,mpfc_avg,curr_compare_idx.*ap.nanout(curr_compare_idx==0)));

        curr_compare_idx_shuff = ap.shake(repmat(curr_compare_idx,1,n_shuff),1,curr_compare_idx>0);
        curr_mpfc_diff_shuff = ...
            (mpfc_avg'*(curr_compare_idx_shuff==2))./sum(curr_compare_idx_shuff==2) - ...
            (mpfc_avg'*(curr_compare_idx_shuff==1))./sum(curr_compare_idx_shuff==1);

        stat_rank = tiedrank([curr_mpfc_diff,curr_mpfc_diff_shuff]');
        stat_p = 1-stat_rank(1)/(n_shuff+1);

        mpfc_p(curr_animal_tds == curr_day) = stat_p;
    end

    mpfc_diffday{curr_animal} = (1:length(curr_animal_tds))' - min([find(mpfc_p < 0.05,1),nan]);
    disp(curr_animal);

end

n_learned_day = cellfun(@(x) max([0, ...
    find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'));

disp([n_learned_day,cellfun(@(x) min([find(x==0),NaN]),mpfc_diffday)]);

% Sub LD for mPFC day
wf_grp.ld = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
        num2cell(cell2mat(mpfc_diffday)),wf.V_event_align,'uni',false));
striatum_mua_grp.ld = cell2mat(cellfun(@(grp,n_trials,n_mua) repmat(grp,n_trials*n_mua,1), ...
        num2cell(cell2mat(mpfc_diffday)),num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

% % (to put back)
% wf_grp.ld = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
%         num2cell(bhv.days_from_learning),wf.V_event_align,'uni',false));
% striatum_mua_grp.ld = cell2mat(cellfun(@(grp,n_trials,n_mua) repmat(grp,n_trials*n_mua,1), ...
%         num2cell(bhv.days_from_learning),num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));


%% Load passive naive

load('\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Marica_2025\data\wf_passive_naive.mat');

% Grab aligned data time
wf_t = wf.wf_stim_time{1};

% Create group indicies
wf_grp = struct;

wf_grp.animal = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
    num2cell(grp2idx(wf.animal)),wf.V_event_align,'uni',false));

wf_grp.stim = cell2mat(wf.trial_stim_values);

cortex_trials_rec_n = cellfun(@(x) size(x,1),wf.V_event_align);















