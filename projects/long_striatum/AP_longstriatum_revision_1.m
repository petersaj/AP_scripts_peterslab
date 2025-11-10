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

%% Passive cortex + striatum PSTH

%%% Load data for figure
load_dataset = 'passive';
Marica_2025.figures.load_data;
%%%

plot_day_bins = [-8,-2,0,Inf];
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


%% Passive cortex + striatum max

%%% Load data for figure
load_dataset = 'passive';
Marica_2025.figures.load_data;
%%%

% plot_day_bins = [-Inf,-2,-1,0,1,Inf];
% plot_day_bins = [-8,-2,0,Inf];
plot_day_bins = [-2,0,Inf];

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
    [wf_grp.animal,cortex_plot_day_grp,wf_grp.stim]);

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
for curr_domain = 1:3
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


%%%% stats - shuffling across, not just within animals
% ~~~ STATS ~~~
% (compare day i to i+1)
print_stat('\n--FIG 3--\n');
print_stat('PSTH max\n');
for curr_compare_day = 1:length(plot_day_bins)-2

    compare_day_grps = curr_compare_day+[0,1];

    cortex_stat_usedata = ismember(wf_striatum_roi_grp(:,2),compare_day_grps);
    [cortex_stat_meas,cortex_stat_grp] = ap.nestgroupfun({@mean,@diff},wf_striatum_roi_max(cortex_stat_usedata,:,:), ...
        wf_striatum_roi_grp(cortex_stat_usedata,2),wf_striatum_roi_grp(cortex_stat_usedata,3));

    striatum_stat_usedata = ismember(striatum_mua_dayavg_grp(:,2),compare_day_grps);
    [striatum_stat_meas,striatum_stat_grp] = ap.nestgroupfun({@mean,@diff},striatum_mua_dayavg_tmax(striatum_stat_usedata), ...
        striatum_mua_dayavg_grp(striatum_stat_usedata,2),striatum_mua_dayavg_grp(striatum_stat_usedata,3:end));

    [~,~,cortex_shuff_grp] = unique(wf_striatum_roi_grp(cortex_stat_usedata,[3]),'rows');
    [~,~,striatum_shuff_grp] = unique(striatum_mua_dayavg_grp(striatum_stat_usedata,[3,4]),'rows');

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

    print_stat('Cortex: day grps %d vs %d\n',compare_day_grps);
    stat_sig = discretize(cortex_stat_p < 0.05,[0,1,Inf],["","*"]);
    for curr_domain = 1:n_domains
        for curr_stim = unique(cortex_stat_grp(:,1))'
            curr_stat_idx = ismember(cortex_stat_grp,curr_stim,'rows');
            print_stat('ROI%d, Stim %3.f, p = %.2g%s\n', ...
                curr_domain,cortex_stat_grp(curr_stat_idx),cortex_stat_p(curr_stat_idx,1,curr_domain),stat_sig(curr_stat_idx,1,curr_domain));
        end
    end

    print_stat('Striatum: day grps %d vs %d\n',compare_day_grps);
    stat_sig = discretize(striatum_stat_p < 0.05,[0,1,Inf],["","*"]);
    for curr_domain = 1:n_domains
        for curr_stim = unique(striatum_stat_grp(:,1))'
            curr_stat_idx = ismember(striatum_stat_grp,[curr_stim,curr_domain],'rows');
            print_stat('D%d, Stim %3.f, p = %.2g%s\n', ...
                curr_domain,striatum_stat_grp(curr_stat_idx,1),striatum_stat_p(curr_stat_idx),stat_sig(curr_stat_idx));
        end
    end
end

%% Passive stats (shuffle days across animals)


curr_compare_day = 1;

compare_day_grps = curr_compare_day+[0,1];

cortex_stat_usedata = ismember(wf_striatum_roi_grp(:,2),compare_day_grps);

dat = ap.nestgroupfun({@mean,@diff},wf_striatum_roi_max(cortex_stat_usedata,:,:), ...
    wf_striatum_roi_grp(cortex_stat_usedata,2),wf_striatum_roi_grp(cortex_stat_usedata,3));

n_shuff = 10000;
dat_shuff = nan(3,n_shuff,3);
for curr_shuff = 1:n_shuff
    dat_shuff(:,curr_shuff,:) = ...
        ap.nestgroupfun({@mean,@diff},wf_striatum_roi_max(cortex_stat_usedata,:,:), ...
        ap.shake(wf_striatum_roi_grp(cortex_stat_usedata,2),1, ...
        wf_striatum_roi_grp(cortex_stat_usedata,3)), ...
        wf_striatum_roi_grp(cortex_stat_usedata,3));
end

stat_rank = tiedrank(permute([dat,dat_shuff],[2,1,3]));
stat_p = permute(1-stat_rank(1,:,:)/(n_shuff+1),[2,3,1]);
% (roi x stim [-90,0,90])
disp(stat_p<0.05);



striatum_stat_usedata = ismember(striatum_mua_dayavg_grp(:,2),compare_day_grps);

[dat,x] = ap.nestgroupfun({@mean,@diff},striatum_mua_dayavg_tmax(striatum_stat_usedata,:,:), ...
    striatum_mua_dayavg_grp(striatum_stat_usedata,2),striatum_mua_dayavg_grp(striatum_stat_usedata,3:4));

n_shuff = 10000;
dat_shuff = nan(9,n_shuff);
for curr_shuff = 1:n_shuff
    [~,~,a] = unique(striatum_mua_dayavg_grp(striatum_stat_usedata,3:4),'rows');

    dat_shuff(:,curr_shuff) = ...
        ap.nestgroupfun({@mean,@diff},striatum_mua_dayavg_tmax(striatum_stat_usedata,:,:), ...
        ap.shake(striatum_mua_dayavg_grp(striatum_stat_usedata,2),1, ...
        a), ...
        a);
end

stat_rank = tiedrank([dat,dat_shuff]');
stat_p = 1-stat_rank(1,:,:)/(n_shuff+1);
% (stim [-90,0,90] x roi)
disp(reshape(stat_p,3,3)<0.05);


a = striatum_mua_dayavg_tmax(ismember(striatum_mua_dayavg_grp(:,2:4),[1,90,1],'rows'));
b = striatum_mua_dayavg_tmax(ismember(striatum_mua_dayavg_grp(:,2:4),[2,90,1],'rows'));


use_roi = 2;

a_idx = ismember(wf_striatum_roi_grp(:,2:3),[compare_day_grps(1),90],'rows');
a = wf_striatum_roi_max(a_idx,:,use_roi);

b_idx = ismember(wf_striatum_roi_grp(:,2:3),[compare_day_grps(2),90],'rows');
b = wf_striatum_roi_max(b_idx,:,use_roi);


data = [a;b];
idx = [ones(size(a));2.*ones(size(b))];

d = diff(ap.groupfun(@mean,data,idx));

n_shuff = 10000;
d_shuff = nan(n_shuff,1);
for curr_shuff = 1:n_shuff
    d_shuff(curr_shuff) = diff(ap.groupfun(@mean,data,ap.shake(idx)));
end

stat_rank = tiedrank([d;d_shuff]);
stat_p = 1-stat_rank(:,1)/(n_shuff+1);
disp(stat_p(1));




%% Task cortex + striatum trial heatmap

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

plot_day_bins = [-inf,-4,-2,1,Inf];
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

% plot_day_bins = [1:6];
% cortex_plot_day_grp = discretize(max(wf_grp.td,-inf),plot_day_bins);
% striatum_plot_day_grp = discretize(max(striatum_mua_grp.td,-inf),plot_day_bins);

% Plot heatmaps sorted by reaction times
heatmap_smooth = [5,1]; % ([trials,time] to smooth for graphics)

figure('Name','Fig 2 heatmaps'); 
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


%% Task cortex + striatum average PSTH w/o movement

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

plot_day_bins = [-2,-1,0,1,Inf];
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% NaN-out activity after movement onset (minus leeway time)
move_leeway = 0.1; % time pre-movement to exclude
striatum_mua_nomove = striatum_mua(:,:,1).*ap.nanout(psth_t > striatum_mua_grp.rxn-move_leeway);
wf_striatum_roi_nomove = wf_striatum_roi(:,:,:,1).*ap.nanout(wf_t > wf_grp.rxn-move_leeway);

% % (to exclude whole trial if short reaction time)
% striatum_mua_nomove = striatum_mua(:,:,1).*ap.nanout(striatum_mua_grp.rxn < 0.3);
% wf_striatum_roi_nomove = wf_striatum_roi(:,:,:,1).*ap.nanout(wf_grp.rxn < 0.3);

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


%% Task cortex + striatum max

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

% plot_day_bins = [-8,-1,0,Inf];
plot_day_bins = [-2,0,Inf];
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
plot(s(:,:,1)',c(:,:,1)','.');
plot(nanmean(s(:,:,1)),nanmean(c(:,:,1)),'k','linewidth',2);

sn = s./s(:,end,:);
cn = c./c(:,end,:);
nexttile; hold on;
plot(sn(:,:,1)',cn(:,:,1)','.')
plot(nanmean(sn(:,:,1)),nanmean(cn(:,:,1)),'k','linewidth',2)


% try averaging trial-trial trends normalized to training time
figure; hold on;
n_points = 50;
data_interp = nan(2,n_points,14);
for use_animal = 1:14

    use_s = striatum_mua_grp.animal==use_animal & ...
        striatum_mua_grp.domain_idx==1 & ...
        striatum_mua_grp.ld > -inf;

    s = mean(striatum_mua(use_s,isbetween(psth_t,0.05,0.15),1),2);

    s_days = unique(striatum_mua_grp.ld(use_s));
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


% ~~~ STATS ~~~

% % NaN-out activity after movement onset (minus leeway time)
% move_leeway = 0.1; % time pre-movement to exclude
% striatum_mua_nomove = striatum_mua(:,:,1).*ap.nanout(psth_t > striatum_mua_grp.rxn-move_leeway);
% wf_striatum_roi_nomove = wf_striatum_roi(:,:,:,1).*ap.nanout(wf_t > wf_grp.rxn-move_leeway);

% (to exclude whole trial if short reaction time)
striatum_mua_nomove = striatum_mua(:,:,1).*ap.nanout(striatum_mua_grp.rxn < rxn_cutoff);
wf_striatum_roi_nomove = wf_striatum_roi(:,:,:,1).*ap.nanout(wf_grp.rxn < rxn_cutoff);

% (ACROSS OR WITHIN ANIMALS)
print_stat('\n--FIG 2--\n');
for curr_compare_day = 1:length(plot_day_bins)-2

    compare_day_grps = curr_compare_day+[0,1];

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
    % [~,~,striatum_shuffgroup] = unique(striatum_mua_nomove_animal_grp(:,[1,3]),'rows'); % within animals
    [~,~,striatum_shuffgroup] = unique(striatum_mua_nomove_animal_grp(:,[3]),'rows'); % across animals
    striatum_mua_nomove_animal_tmax = max(striatum_mua_nomove_animal(:,striatum_stim_t),[],2);

    striatum_stat_usedata = ismember(striatum_mua_nomove_animal_grp(:,2),compare_day_grps);
    striatum_stat_meas = ap.nestgroupfun({@mean,@diff},striatum_mua_nomove_animal_tmax(striatum_stat_usedata), ...
        striatum_mua_nomove_animal_grp(striatum_stat_usedata,2),striatum_mua_nomove_animal_grp(striatum_stat_usedata,3));

    n_shuff = 10000;
    cortex_stat_null = nan(n_domains,n_shuff);
    striatum_stat_null = nan(n_domains,n_shuff);
    for curr_shuff = 1:n_shuff

        % curr_ctx_shuff = ... % within animals
        %     ap.shake(wf_striatum_roi_nomove_animal_tmax(cortex_stat_usedata,:,:),1, ...
        %     wf_striatum_roi_nomove_animal_grp(cortex_stat_usedata,1));
        curr_ctx_shuff = ... % across animals
            ap.shake(wf_striatum_roi_nomove_animal_tmax(cortex_stat_usedata,:,:),1);

        cortex_stat_null(:,curr_shuff) = ...
            permute(diff(ap.groupfun(@mean,curr_ctx_shuff, ...
            wf_striatum_roi_nomove_animal_grp(cortex_stat_usedata,2)),[],1),[3,2,1]);     

        curr_str_shuff = ...
            ap.shake(striatum_mua_nomove_animal_tmax(striatum_stat_usedata),1, ...
            striatum_shuffgroup(striatum_stat_usedata));
        striatum_stat_null(:,curr_shuff) = ...
            ap.nestgroupfun({@mean,@diff},curr_str_shuff, ...
            striatum_mua_nomove_animal_grp(striatum_stat_usedata,2),striatum_mua_nomove_animal_grp(striatum_stat_usedata,3));

    end

    print_stat('Average PSTH tmax, day grps %d,%d:\n',compare_day_grps);

    cortex_stat_rank = tiedrank([cortex_stat_meas,cortex_stat_null]')';
    cortex_stat_p = 1-cortex_stat_rank(:,1)/(n_shuff+1);
    cortex_stat_sig = discretize(cortex_stat_p < 0.05,[0,1,Inf],{'','*'});
    for curr_domain = 1:n_domains
        print_stat('CTX%d p = %.2g%s\n',curr_domain,cortex_stat_p(curr_domain),cortex_stat_sig{curr_domain});
    end

    striatum_stat_rank = tiedrank([striatum_stat_meas,striatum_stat_null]')';
    striatum_stat_p = 1-striatum_stat_rank(:,1)/(n_shuff+1);
    striatum_stat_sig = discretize(striatum_stat_p < 0.05,[0,1,Inf],{'','*'});
    for curr_domain = 1:n_domains
        print_stat('STR%d p = %.2g%s\n',curr_domain,striatum_stat_p(curr_domain),striatum_stat_sig{curr_domain});
    end

end


%% Task stats: shuffle days across animals (not just within)

plot_day_bins = [-8,-1,0,Inf];
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);


% NaN-out activity after movement onset (minus leeway time)
move_leeway = 0.1; % time pre-movement to exclude
striatum_mua_nomove = striatum_mua(:,:,1).*ap.nanout(psth_t > striatum_mua_grp.rxn-move_leeway);
wf_striatum_roi_nomove = wf_striatum_roi(:,:,:,1).*ap.nanout(wf_t > wf_grp.rxn-move_leeway);

% % (to exclude whole trial if short reaction time)
% striatum_mua_nomove = striatum_mua(:,:,1).*ap.nanout(striatum_mua_grp.rxn < 0.3);
% wf_striatum_roi_nomove = wf_striatum_roi(:,:,:,1).*ap.nanout(wf_grp.rxn < 0.3);


curr_compare_day = 1;

compare_day_grps = curr_compare_day+[0,1];

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

print_stat('Average PSTH tmax, day grps %d,%d:\n',compare_day_grps);

cortex_stat_rank = tiedrank([cortex_stat_meas,cortex_stat_null]')';
cortex_stat_p = 1-cortex_stat_rank(:,1)/(n_shuff+1);
cortex_stat_sig = discretize(cortex_stat_p < 0.05,[0,1,Inf],{'','*'});
for curr_domain = 1:n_domains
    print_stat('CTX%d p = %.2g%s\n',curr_domain,cortex_stat_p(curr_domain),cortex_stat_sig{curr_domain});
end

striatum_stat_rank = tiedrank([striatum_stat_meas,striatum_stat_null]')';
striatum_stat_p = 1-striatum_stat_rank(:,1)/(n_shuff+1);
striatum_stat_sig = discretize(striatum_stat_p < 0.05,[0,1,Inf],{'','*'});
for curr_domain = 1:n_domains
    print_stat('STR%d p = %.2g%s\n',curr_domain,striatum_stat_p(curr_domain),striatum_stat_sig{curr_domain});
end




% Testing out shuffle day stats (this works - incorporated above)

data_idx = striatum_mua_nomove_animal_grp(:,2) == 1 & striatum_mua_nomove_animal_grp(:,3) == 2;
a = max(striatum_mua_nomove_animal(data_idx,striatum_stim_t),[],2);

data_idx = striatum_mua_nomove_animal_grp(:,2) == 2 & striatum_mua_nomove_animal_grp(:,3) == 2;
b = max(striatum_mua_nomove_animal(data_idx,striatum_stim_t),[],2);

% data_idx = wf_striatum_roi_nomove_animal_grp(:,2) == 1;
% a = max(wf_striatum_roi_nomove_animal(data_idx,cortex_stim_t,2),[],2);
% 
% data_idx = wf_striatum_roi_nomove_animal_grp(:,2) == 2;
% b = max(wf_striatum_roi_nomove_animal(data_idx,cortex_stim_t,2),[],2);

data = [a;b];
idx = [ones(size(a));2.*ones(size(b))];

d = diff(ap.groupfun(@mean,data,idx));

n_shuff = 10000;
d_shuff = nan(n_shuff,1);
for curr_shuff = 1:n_shuff
    d_shuff(curr_shuff) = diff(ap.groupfun(@mean,data,ap.shake(idx)));
end

stat_rank = tiedrank([d;d_shuff]);
stat_p = 1-stat_rank(:,1)/(n_shuff+1);
disp(stat_p(1));




%% Cortex movie passive

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-2,-1,0,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% plot_day_bins = [1:9,inf];
% plot_day_grp = discretize(max(wf_grp.td,-inf),plot_day_bins);

wf_V_event_align_cat = cell2mat(wf.V_event_align);
wf_trial_stim_values_cat = cell2mat(wf.trial_stim_values);


use_trials = wf_grp.animal == 2 & wf_trial_stim_values_cat==90;
% use_trials = true(size(wf_grp.animal));
% use_trials = wf_trial_stim_values_cat == 90;

% n_learned_day = cellfun(@(x) max([0, ...
%     find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
%     unique(bhv.animal,'stable'));
% use_trials = wf_trial_stim_values_cat==90 & n_learned_day(wf_grp.animal) > 3;


[wf_avg,wf_avg_grp] = ap.nestgroupfun({@nanmean,@nanmean},wf_V_event_align_cat(use_trials,:,:), ...
    wf_grp.animal(use_trials),[plot_day_grp(use_trials),wf_trial_stim_values_cat(use_trials)]);

% [wf_avg,wf_avg_grp] = ap.groupfun(@nanmean,wf_V_event_align_cat(use_trials,:,:), ...
%     [wf_grp.animal(use_trials),plot_day_grp(use_trials),wf_trial_stim_values_cat(use_trials)]);

baseline_t = [-0.05,-0];
baseline_t_idx = isbetween(wf_t,baseline_t(1),baseline_t(2));
wf_px = plab.wf.svd2px(U_master,permute(wf_avg-nanmean(wf_avg(:,baseline_t_idx,:),2),[3,2,1]));

ap.imscroll(wf_px,wf_t);
colormap(ap.colormap('PWG',[],1.5));
clim([-1,1]*3e-3);
axis image off;
ap.wf_draw('ccf',[0.5,0.5,0.5]);

figure;
imagesc(reshape(squeeze(max(wf_px(:,:,isbetween(wf_t,0,0.5),:),[],3)),size(wf_px,1),[]))
colormap(ap.colormap('PWG',[],1.5));
clim([-1,1]*3e-3);
axis image off;


%% Cortex movie task

%%% Load data for figure
load_dataset = 'task';
AP_longstriatum_load_data;
%%%

plot_day_bins = [-2,-1,0,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

move_leeway = 0.1; % time pre-movement to exclude
wf_V_event_align_cat = cell2mat(cellfun(@(x) x(:,:,:,1), ...
    wf.V_event_align,'uni',false)).*ap.nanout(wf_t > (wf_grp.rxn-move_leeway));

use_trials = true(size(wf_grp.animal)) & wf_grp.animal == 2;

% n_learned_day = cellfun(@(x) max([0, ...
%     find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
%     unique(bhv.animal,'stable'));
% use_trials = n_learned_day(wf_grp.animal)>3;

% [wf_avg,wf_avg_grp] = ap.groupfun(@nanmean, ...
%     wf_V_event_align_cat(use_trials,:,:), ...
%     [wf_grp.animal(use_trials),plot_day_grp(use_trials)]);

[wf_avg,wf_avg_grp] = ap.nestgroupfun({@nanmean,@nanmean}, ...
    wf_V_event_align_cat(use_trials,:,:), ...
    wf_grp.animal(use_trials),[plot_day_grp(use_trials)]);

baseline_t = [-0.2,-0.1];
baseline_t_idx = isbetween(wf_t,baseline_t(1),baseline_t(2));
wf_px = plab.wf.svd2px(U_master,permute(wf_avg-nanmean(wf_avg(:,baseline_t_idx,:),2),[3,2,1]));

ap.imscroll(wf_px,wf_t);
colormap(ap.colormap('PWG',[],1.5));
clim([-1,1]*3e-3);
axis image off;
ap.wf_draw('ccf',[0.5,0.5,0.5]);

figure;
imagesc(reshape(squeeze(max(wf_px(:,:,isbetween(wf_t,0,0.5),:),[],3)),size(wf_px,1),[]))
colormap(ap.colormap('PWG',[],1.5));
clim([-1,1]*3e-3);
axis image off;


%% mPFC TASK: get day it's different from first

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

rxn_cutoff = 0.3;

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
        curr_compare_idx_shuff = ap.shake(repmat(curr_compare_idx,1,n_shuff),1,curr_compare_idx>0);

        % % trial max, mean
        % curr_mpfc_diff = diff(ap.groupfun(@mean,mpfc_avg,curr_compare_idx.*ap.nanout(curr_compare_idx==0)));    
        % curr_mpfc_diff_shuff = ...
        %     (mpfc_avg'*(curr_compare_idx_shuff==2))./sum(curr_compare_idx_shuff==2) - ...
        %     (mpfc_avg'*(curr_compare_idx_shuff==1))./sum(curr_compare_idx_shuff==1);

        % trial mean, max
        curr_mpfc_diff = diff(max(ap.groupfun(@mean,wf_striatum_roi(:,cortex_stim_t,2), ...
            curr_compare_idx.*ap.nanout(curr_compare_idx==0)),[],2));
        curr_mpfc_diff_shuff = nan(1,n_shuff);
        for curr_shuff = 1:n_shuff
            curr_mpfc_diff_shuff(curr_shuff) = diff(max(ap.groupfun(@mean,wf_striatum_roi(:,cortex_stim_t,2), ...
                curr_compare_idx_shuff(:,curr_shuff).*ap.nanout(curr_compare_idx_shuff(:,curr_shuff)==0)),[],2));
        end

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
n_mpfc_diff_day = cellfun(@(x) min([find(x==0),NaN]),mpfc_diffday);
disp([n_learned_day,n_mpfc_diff_day,n_learned_day-n_mpfc_diff_day]);

% %%%% (to replace no diff days with learning day)
% for curr_animal = find(isnan(n_mpfc_diff_day))'
%     mpfc_diffday{curr_animal} = (1:length(mpfc_diffday{curr_animal}))'-n_learned_day(curr_animal);
% end
% %%%%

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

%% mPFC TASK: get day it's different from baseline time

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

rxn_cutoff = 0.3;

% Average day, get max in time window
stim_t = [0.05,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
cortex_baseline_t = isbetween(wf_t,-stim_t(2),-stim_t(1));

td_cat = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));

wf_grp.td = cell2mat(cellfun(@(animal,data) ...
    repmat(animal,size(data,1),1), ...
    num2cell(td_cat),wf.V_event_align,'uni',false));

mpfc_diffday = cell(length(unique(bhv.animal)),1);
for curr_animal = 1:length(unique(bhv.animal,'stable'))

    curr_animal_tds = unique(wf_grp.td(wf_grp.animal == curr_animal));

    n_shuff = 1000;
    mpfc_p = nan(size(curr_animal_tds));
    for curr_day = curr_animal_tds(2:end)'

        curr_day_trials = wf_grp.animal == curr_animal & ...
            wf_grp.td == curr_day & ...
            wf_grp.rxn >= rxn_cutoff;

        if sum(curr_day_trials) < 10
            continue
        end

        % trial mean, max
        curr_act = [...
            mat2cell(wf_striatum_roi(curr_day_trials,cortex_baseline_t,2),ones(sum(curr_day_trials),1),sum(cortex_baseline_t)), ...
            mat2cell(wf_striatum_roi(curr_day_trials,cortex_stim_t,2),ones(sum(curr_day_trials),1),sum(cortex_stim_t))];

        curr_mpfc_diff = diff(arrayfun(@(x) max(nanmean(vertcat(curr_act{:,x}),1)),1:2));

        curr_mpfc_diff_shuff = nan(1,n_shuff);
        for curr_shuff = 1:n_shuff
            curr_act_shuff = ap.shake(curr_act,2);
            curr_mpfc_diff_shuff(curr_shuff) = ...
                diff(arrayfun(@(x) max(nanmean(vertcat(curr_act_shuff{:,x}),1)),1:2));
        end

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
n_mpfc_diff_day = cellfun(@(x) min([find(x==0),NaN]),mpfc_diffday);
disp([n_learned_day,n_mpfc_diff_day,n_learned_day-n_mpfc_diff_day]);

% %%%% (to replace no diff days with learning day)
% for curr_animal = find(isnan(n_mpfc_diff_day))'
%     mpfc_diffday{curr_animal} = (1:length(mpfc_diffday{curr_animal}))'-n_learned_day(curr_animal);
% end
% %%%%

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


%% mPFC TASK: get day it's different from R hemi

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

rxn_cutoff = 0.3;

% Average day, get max in time window
stim_t = [0.05,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
cortex_baseline_t = isbetween(wf_t,-stim_t(2),-stim_t(1));

td_cat = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));

wf_grp.td = cell2mat(cellfun(@(animal,data) ...
    repmat(animal,size(data,1),1), ...
    num2cell(td_cat),wf.V_event_align,'uni',false));

mpfc_diffday = cell(length(unique(bhv.animal)),1);
for curr_animal = 1:length(unique(bhv.animal,'stable'))

    curr_animal_tds = unique(wf_grp.td(wf_grp.animal == curr_animal));

    n_shuff = 1000;
    mpfc_p = nan(size(curr_animal_tds));
    for curr_day = curr_animal_tds(2:end)'

        curr_day_trials = wf_grp.animal == curr_animal & ...
            wf_grp.td == curr_day & ...
            wf_grp.rxn >= rxn_cutoff;

        if sum(curr_day_trials) < 10
            continue
        end

        % trial mean, max
        curr_act = [...
            mat2cell(wf_striatum_roi(curr_day_trials,cortex_baseline_t,5,1),ones(sum(curr_day_trials),1),sum(cortex_stim_t)), ...
            mat2cell(wf_striatum_roi(curr_day_trials,cortex_stim_t,2,1),ones(sum(curr_day_trials),1),sum(cortex_stim_t))];

        curr_mpfc_diff = diff(arrayfun(@(x) max(nanmean(vertcat(curr_act{:,x}),1)),1:2));

        curr_mpfc_diff_shuff = nan(1,n_shuff);
        for curr_shuff = 1:n_shuff
            curr_act_shuff = ap.shake(curr_act,2);
            curr_mpfc_diff_shuff(curr_shuff) = ...
                diff(arrayfun(@(x) max(nanmean(vertcat(curr_act_shuff{:,x}),1)),1:2));
        end

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
n_mpfc_diff_day = cellfun(@(x) min([find(x==0),NaN]),mpfc_diffday);
disp([n_learned_day,n_mpfc_diff_day,n_learned_day-n_mpfc_diff_day]);

% %%%% (to replace no diff days with learning day)
% for curr_animal = find(isnan(n_mpfc_diff_day))'
%     mpfc_diffday{curr_animal} = (1:length(mpfc_diffday{curr_animal}))'-n_learned_day(curr_animal);
% end
% %%%%

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

% trial max,mean or trial mean,max
% vs naive or TD1

%%% Load data for figure
load_dataset = 'passive';
Marica_2025.figures.load_data;
%%%

%%% Load passive naive
wf_passive_naive = load("C:\Users\petersa\Desktop\wf_passive_naive.mat");
%%%

use_stim = 90; 

stim_t = [0,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

td_cat = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));

wf_grp.td = cell2mat(cellfun(@(animal,data) ...
    repmat(animal,size(data,1),1), ...
    num2cell(td_cat),wf.V_event_align,'uni',false));

mpfc_diffday = cell(length(unique(bhv.animal)),1);
for curr_animal = 1:length(unique(bhv.animal,'stable'))

    curr_animal_tds = unique(wf_grp.td(wf_grp.animal == curr_animal));

    % % (baseline = naive)
    % curr_baseline_data = wf_passive_naive.wf_striatum_roi(...
    %     wf_passive_naive.wf_grp.animal == curr_animal & ...
    %     wf_passive_naive.wf_grp.stim == use_stim,:,2);

    % % (baseline = naive hemidiff)
    % curr_baseline_data = wf_passive_naive.wf_striatum_roi(...
    %     wf_passive_naive.wf_grp.animal == curr_animal & ...
    %     wf_passive_naive.wf_grp.stim == use_stim,:,2) - ...
    %     wf_passive_naive.wf_striatum_roi(...
    %     wf_passive_naive.wf_grp.animal == curr_animal & ...
    %     wf_passive_naive.wf_grp.stim == use_stim,:,5);

    % % (baseline = day 1)
    % curr_baseline_data = wf_striatum_roi(...
    %     wf_grp.td == 1 & ...
    %     wf_grp.animal == curr_animal & ...
    %     wf_grp.stim == use_stim,:,2);

     % (baseline = day 1 hemidiff)
    curr_baseline_data = wf_striatum_roi(...
        wf_grp.td == 1 & ...
        wf_grp.animal == curr_animal & ...
        wf_grp.stim == use_stim,:,2) - ...
        wf_striatum_roi(...
        wf_grp.td == 1 & ...
        wf_grp.animal == curr_animal & ...
        wf_grp.stim == use_stim,:,5);


    n_shuff = 1000;
    mpfc_p = nan(size(curr_animal_tds));
    for curr_day = curr_animal_tds'

        % curr_day_data = wf_striatum_roi( ...
        %     wf_grp.animal == curr_animal & ...
        %     wf_grp.td == curr_day & ...
        %     wf_grp.stim == use_stim,:,2);

        % (hemidiff)
        curr_day_data = wf_striatum_roi( ...
            wf_grp.animal == curr_animal & ...
            wf_grp.td == curr_day & ...
            wf_grp.stim == use_stim,:,2) - ...
            wf_striatum_roi( ...
            wf_grp.animal == curr_animal & ...
            wf_grp.td == curr_day & ...
            wf_grp.stim == use_stim,:,5);

        if size(curr_day_data,1) < 10
            continue
        end

        curr_data_cat = [curr_baseline_data;curr_day_data];
        curr_data_idx = [ones(size(curr_baseline_data,1),1); ...
            2.*ones(size(curr_day_data,1),1)];

        % (trial mean, max)
        mpfc_diff = diff(max(ap.groupfun(@mean, ...
            curr_data_cat(:,cortex_stim_t),curr_data_idx),[],2));
        mpfc_diff_shuff = nan(n_shuff,1);
        for curr_shuff = 1:n_shuff
            mpfc_diff_shuff(curr_shuff) = ...
                diff(max(ap.groupfun(@mean, ...
                curr_data_cat(:,cortex_stim_t),ap.shake(curr_data_idx)),[],2));
        end

        % % (trial max, median)
        % mpfc_diff = diff(ap.groupfun(@median, ...
        %     max(curr_data_cat(:,cortex_stim_t),[],2),curr_data_idx));
        % mpfc_diff_shuff = nan(n_shuff,1);
        % for curr_shuff = 1:n_shuff
        %     mpfc_diff_shuff(curr_shuff) = ...
        %         diff(ap.groupfun(@median, ...
        %         max(curr_data_cat(:,cortex_stim_t),[],2),ap.shake(curr_data_idx)));
        % end

        stat_rank = tiedrank([mpfc_diff;mpfc_diff_shuff]);
        stat_p = 1-stat_rank(1)/(n_shuff+1);

        mpfc_p(curr_animal_tds == curr_day) = stat_p;
    end

    mpfc_diffday{curr_animal} = (1:length(curr_animal_tds))' - min([find(mpfc_p < 0.05,1),nan]);
    disp(curr_animal);

end

n_learned_day = cellfun(@(x) max([0, ...
    find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'));
n_mpfc_diff_day = cellfun(@(x) min([find(x==0),NaN]),mpfc_diffday);
disp([n_learned_day,n_mpfc_diff_day,n_learned_day-n_mpfc_diff_day]);

% %%%% (to replace no diff days with learning day)
% for curr_animal = find(isnan(n_mpfc_diff_day))'
%     mpfc_diffday{curr_animal} = (1:length(mpfc_diffday{curr_animal}))'-n_learned_day(curr_animal);
% end
% %%%%

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


%% mPFC PASSIVE: get day R different from L

% trial max,mean or trial mean,max
% vs naive or TD1

%%% Load data for figure
load_dataset = 'passive';
Marica_2025.figures.load_data;
%%%

%%% Load passive naive
wf_passive_naive = load("C:\Users\petersa\Desktop\wf_passive_naive.mat");
%%%

stim_t = [0,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

td_cat = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));

wf_grp.td = cell2mat(cellfun(@(animal,data) ...
    repmat(animal,size(data,1),1), ...
    num2cell(td_cat),wf.V_event_align,'uni',false));

mpfc_diffday = cell(length(unique(bhv.animal)),1);
for curr_animal = 1:length(unique(bhv.animal,'stable'))

    curr_animal_tds = unique(wf_grp.td(wf_grp.animal == curr_animal));

    n_shuff = 1000;
    mpfc_p = nan(size(curr_animal_tds));
    for curr_day = curr_animal_tds'

        curr_baseline_data = wf_striatum_roi( ...
            wf_grp.animal == curr_animal & ...
            wf_grp.td == curr_day & ...
            wf_grp.stim == -90,:,2);

        curr_day_data = wf_striatum_roi( ...
            wf_grp.animal == curr_animal & ...
            wf_grp.td == curr_day & ...
            wf_grp.stim == 90,:,2);


        if size(curr_day_data,1) < 10
            continue
        end

        curr_data_cat = [curr_baseline_data;curr_day_data];
        curr_data_idx = [ones(size(curr_baseline_data,1),1); ...
            2.*ones(size(curr_day_data,1),1)];

        % (trial mean, max)
        mpfc_diff = diff(max(ap.groupfun(@mean, ...
            curr_data_cat(:,cortex_stim_t),curr_data_idx),[],2));
        mpfc_diff_shuff = nan(n_shuff,1);
        for curr_shuff = 1:n_shuff
            mpfc_diff_shuff(curr_shuff) = ...
                diff(max(ap.groupfun(@mean, ...
                curr_data_cat(:,cortex_stim_t),ap.shake(curr_data_idx)),[],2));
        end

        % % (trial max, median)
        % mpfc_diff = diff(ap.groupfun(@median, ...
        %     max(curr_data_cat(:,cortex_stim_t),[],2),curr_data_idx));
        % mpfc_diff_shuff = nan(n_shuff,1);
        % for curr_shuff = 1:n_shuff
        %     mpfc_diff_shuff(curr_shuff) = ...
        %         diff(ap.groupfun(@median, ...
        %         max(curr_data_cat(:,cortex_stim_t),[],2),ap.shake(curr_data_idx)));
        % end

        stat_rank = tiedrank([mpfc_diff;mpfc_diff_shuff]);
        stat_p = 1-stat_rank(1)/(n_shuff+1);

        mpfc_p(curr_animal_tds == curr_day) = stat_p;
    end

    mpfc_diffday{curr_animal} = (1:length(curr_animal_tds))' - min([find(mpfc_p < 0.05,1),nan]);
    disp(curr_animal);

end

n_learned_day = cellfun(@(x) max([0, ...
    find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'));
n_mpfc_diff_day = cellfun(@(x) min([find(x==0),NaN]),mpfc_diffday);
disp([n_learned_day,n_mpfc_diff_day,n_learned_day-n_mpfc_diff_day]);

% %%%% (to replace no diff days with learning day)
% for curr_animal = find(isnan(n_mpfc_diff_day))'
%     mpfc_diffday{curr_animal} = (1:length(mpfc_diffday{curr_animal}))'-n_learned_day(curr_animal);
% end
% %%%%

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



%% Widefield striatum ROIs (L/R hemi)

% Create ROIs by striatum cluster maps
kmeans_centroid_blur = imgaussfilt(kmeans_cluster_mean,10);
striatum_wf_roi = kmeans_centroid_blur > prctile(kmeans_centroid_blur,100,[1,2])*0.75;

%%%% second set of L/R flip rois
striatum_wf_roi = cat(3,striatum_wf_roi,fliplr(striatum_wf_roi));

% Get ROI activity
wf_striatum_roi_lr = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master,permute(x,[3,2,1,4]),[],[],striatum_wf_roi), ...
    [3,2,1,4]),wf.V_event_align(~cellfun(@isempty,wf.V_event_align)),'uni',false));

% hemidiff
wf_striatum_roi = wf_striatum_roi_lr(:,:,1:3,:) - wf_striatum_roi_lr(:,:,4:6,:);

% (baseline-subtract)
baseline_t = wf_t < 0;
wf_striatum_roi = wf_striatum_roi - nanmean(wf_striatum_roi(:,baseline_t,:,:),2);


%% Make indicies for training day

td_cat = cell2mat(cellfun(@(x) (1:sum(strcmp(bhv.animal,x)))', ...
    unique(bhv.animal,'stable'),'uni',false));

wf_grp.td = cell2mat(cellfun(@(animal,data) ...
    repmat(animal,size(data,1),1), ...
    num2cell(td_cat),wf.V_event_align,'uni',false));

striatum_mua_grp.td = cell2mat(cellfun(@(grp,n_trials,n_mua) ...
    repmat(grp,n_trials*n_mua,1),num2cell(td_cat), ...
    num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));








