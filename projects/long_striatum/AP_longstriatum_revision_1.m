%% Exploratory/testing for revision

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

plot_day_bins = [-10,-1,0,Inf];
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

plot_day_bins = [-Inf,-2,0,2];
% plot_day_bins = [-10,-1,0,Inf];
% plot_day_bins = -10:10;

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

%% Passive trial-smoothed

% try averaging trial-trial trends normalized to training time
figure; hold on;
n_points = 50;
data_interp = nan(2,n_points,14);
for use_animal = 1:14

    use_s = striatum_mua_grp.animal==use_animal & ...
        striatum_mua_grp.domain_idx==1 & ...
        striatum_mua_grp.ld > -inf & ...
        striatum_mua_grp.stim == 90;

    s = mean(striatum_mua(use_s,isbetween(psth_t,0.05,0.15),1),2);

    s_days = unique(striatum_mua_grp.ld(use_s));
    c = mean(wf_striatum_roi(wf_grp.stim == 90 & wf_grp.animal==use_animal & ismember(wf_grp.ld,s_days),isbetween(wf_t,0.05,0.15),2),2);

    if isempty(c)
        continue
    end

    s_smooth = smoothdata(s,'movmean',20);
    c_smooth = smoothdata(c,'movmean',20);
    plot(smoothdata(s,'movmean',100),smoothdata(c,'movmean',100));

    data_interp(1,:,use_animal) = s_smooth(round(linspace(1,length(s_smooth),n_points)));
    data_interp(2,:,use_animal) = c_smooth(round(linspace(1,length(c_smooth),n_points)));
end

figure; hold on;
plot(permute(data_interp(1,:,:),[2,3,1]),permute(data_interp(2,:,:),[2,3,1]));
plot(nanmean(data_interp(1,:,:),3),nanmean(data_interp(2,:,:),3),'k','linewidth',2)



%% Task cortex + striatum trial heatmap

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

% plot_day_bins = [-Inf,-1,0,Inf];
plot_day_bins = [-10,-1,0,Inf];
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
for curr_day = 1:length(plot_day_bins)-1
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

plot_day_bins = [-10,-1,0,Inf];
% plot_day_bins = [-Inf,-2,0,Inf];
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

plot_day_bins = [-10,-1,0,Inf];
% plot_day_bins = [-Inf,-1,0,Inf];
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
    [~,~,striatum_shuffgroup] = unique(striatum_mua_nomove_animal_grp(:,[1,3]),'rows'); % within animals
    % [~,~,striatum_shuffgroup] = unique(striatum_mua_nomove_animal_grp(:,[3]),'rows'); % across animals
    striatum_mua_nomove_animal_tmax = max(striatum_mua_nomove_animal(:,striatum_stim_t),[],2);

    striatum_stat_usedata = ismember(striatum_mua_nomove_animal_grp(:,2),compare_day_grps);
    striatum_stat_meas = ap.nestgroupfun({@mean,@diff},striatum_mua_nomove_animal_tmax(striatum_stat_usedata), ...
        striatum_mua_nomove_animal_grp(striatum_stat_usedata,2),striatum_mua_nomove_animal_grp(striatum_stat_usedata,3));

    n_shuff = 10000;
    cortex_stat_null = nan(n_domains,n_shuff);
    striatum_stat_null = nan(n_domains,n_shuff);
    for curr_shuff = 1:n_shuff

        curr_ctx_shuff = ... % within animals
            ap.shake(wf_striatum_roi_nomove_animal_tmax(cortex_stat_usedata,:,:),1, ...
            wf_striatum_roi_nomove_animal_grp(cortex_stat_usedata,1));
        % curr_ctx_shuff = ... % across animals
        %     ap.shake(wf_striatum_roi_nomove_animal_tmax(cortex_stat_usedata,:,:),1);

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

% plot_day_bins = [-8,-1,0,Inf];
plot_day_bins = [-8,-2,0,2];
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

plot_day_bins = [-Inf,0,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% plot_day_bins = [1:9,inf];
% plot_day_grp = discretize(max(wf_grp.td,-inf),plot_day_bins);

wf_V_event_align_cat = cell2mat(wf.V_event_align);
wf_trial_stim_values_cat = cell2mat(wf.trial_stim_values);


% use_trials = wf_grp.animal == 2 & wf_trial_stim_values_cat==90;
% use_trials = true(size(wf_grp.animal));
use_trials = wf_trial_stim_values_cat == 90;

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

plot_day_bins = [-Inf,0,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

move_leeway = 0.1; % time pre-movement to exclude
wf_V_event_align_cat = cell2mat(cellfun(@(x) x(:,:,:,1), ...
    wf.V_event_align,'uni',false)).*ap.nanout(wf_t > (wf_grp.rxn-move_leeway));

% use_trials = true(size(wf_grp.animal));
use_trials = wf_grp.rxn > 0.3;

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


%% striatum TASK: get day it's different from previous
% (not working at the moment)

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

rxn_cutoff = 0.3;

% Average day, get max in time window
stim_t = [0,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));


str_diffday = cell(length(unique(bhv.animal)),1);
for curr_animal = 1:length(unique(bhv.animal,'stable'))

    curr_animal_tds = unique(wf_grp.td(wf_grp.animal == curr_animal));

    n_shuff = 500;
    mpfc_p = nan(size(curr_animal_tds));
    for curr_day = curr_animal_tds(2:end)'

        curr_baseline_trials = striatum_mua_grp.animal == curr_animal & ...
            striatum_mua_grp.td < curr_day & ...
            striatum_mua_grp.rxn >= rxn_cutoff & ...
            striatum_mua_grp.domain_idx == 1;

        curr_day_trials = striatum_mua_grp.animal == curr_animal & ...
            striatum_mua_grp.td == curr_day & ...
            striatum_mua_grp.rxn >= rxn_cutoff  &...
            striatum_mua_grp.domain_idx == 1;

        if sum(curr_baseline_trials) < 10 || sum(curr_day_trials) < 10
            continue
        end

        curr_compare_idx = curr_baseline_trials + curr_day_trials.*2;
        curr_compare_idx_shuff = ap.shake(repmat(curr_compare_idx,1,n_shuff),1,curr_compare_idx>0);

        % trial mean, max
        curr_mpfc_diff = diff(max(ap.groupfun(@mean,striatum_mua(:,striatum_stim_t,1), ...
            curr_compare_idx.*ap.nanout(curr_compare_idx==0)),[],2));
        curr_mpfc_diff_shuff = nan(1,n_shuff);
        for curr_shuff = 1:n_shuff
            curr_mpfc_diff_shuff(curr_shuff) = diff(max(ap.groupfun(@mean,striatum_mua(:,striatum_stim_t,1), ...
                curr_compare_idx_shuff(:,curr_shuff).*ap.nanout(curr_compare_idx_shuff(:,curr_shuff)==0)),[],2));
        end

        stat_rank = tiedrank([curr_mpfc_diff,curr_mpfc_diff_shuff]');
        stat_p = 1-stat_rank(1)/(n_shuff+1);

        mpfc_p(curr_animal_tds == curr_day) = stat_p;
    end

    str_diffday{curr_animal} = (1:length(curr_animal_tds))' - min([find(mpfc_p < 0.01,1),nan]);
    disp(curr_animal);

end

n_learned_day = cellfun(@(x) max([0, ...
    find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'));
n_mpfc_diff_day = cellfun(@(x) min([find(x==0),NaN]),str_diffday);
disp([n_learned_day,n_mpfc_diff_day,n_learned_day-n_mpfc_diff_day]);

% %%%% (to replace no diff days with learning day)
% for curr_animal = find(isnan(n_mpfc_diff_day))'
%     mpfc_diffday{curr_animal} = (1:length(mpfc_diffday{curr_animal}))'-n_learned_day(curr_animal);
% end
% %%%%

% Sub LD for mPFC day
wf_grp.ld = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
    num2cell(cell2mat(str_diffday)),wf.V_event_align,'uni',false));
striatum_mua_grp.ld = cell2mat(cellfun(@(grp,n_trials,n_mua) repmat(grp,n_trials*n_mua,1), ...
    num2cell(cell2mat(str_diffday)),num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

% % (to put back)
% wf_grp.ld = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
%         num2cell(bhv.days_from_learning),wf.V_event_align,'uni',false));
% striatum_mua_grp.ld = cell2mat(cellfun(@(grp,n_trials,n_mua) repmat(grp,n_trials*n_mua,1), ...
%         num2cell(bhv.days_from_learning),num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));


%% mPFC TASK: get day it's different from first N

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

rxn_cutoff = 0.3;

% Average day, get max in time window
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

    % curr_baseline_trials = wf_grp.animal == curr_animal & ...
    %     ismember(wf_grp.td,[1,2]) & ...
    %     wf_grp.rxn >= rxn_cutoff;

    n_shuff = 2000;
    mpfc_p = nan(size(curr_animal_tds));
    for curr_day = curr_animal_tds(2:end)'

        curr_baseline_trials = wf_grp.animal == curr_animal & ...
            wf_grp.td < curr_day & ...
            wf_grp.rxn >= rxn_cutoff;

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
        curr_mpfc_diff = diff(max(ap.groupfun(@mean,wf_striatum_roi(:,cortex_stim_t,2,1), ...
            curr_compare_idx.*ap.nanout(curr_compare_idx==0)),[],2));
        curr_mpfc_diff_shuff = nan(1,n_shuff);
        for curr_shuff = 1:n_shuff
            curr_mpfc_diff_shuff(curr_shuff) = diff(max(ap.groupfun(@mean,wf_striatum_roi(:,cortex_stim_t,2,1), ...
                curr_compare_idx_shuff(:,curr_shuff).*ap.nanout(curr_compare_idx_shuff(:,curr_shuff)==0)),[],2));
        end

        stat_rank = tiedrank([curr_mpfc_diff,curr_mpfc_diff_shuff]');
        stat_p = 1-stat_rank(1)/(n_shuff+1);

        mpfc_p(curr_animal_tds == curr_day) = stat_p;
    end

    mpfc_diffday{curr_animal} = (1:length(curr_animal_tds))' - min([find(mpfc_p < 0.01,1),nan]);
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

% (to put back)
wf_grp.ld = cell2mat(cellfun(@(animal,data) repmat(animal,size(data,1),1), ...
    num2cell(bhv.days_from_learning),wf.V_event_align,'uni',false));
striatum_mua_grp.ld = cell2mat(cellfun(@(grp,n_trials,n_mua) repmat(grp,n_trials*n_mua,1), ...
    num2cell(bhv.days_from_learning),num2cell(striatum_trials_rec_n),num2cell(striatum_mua_rec_n),'uni',false));

%% *mPFC TASK: get day it's different from baseline time (do hemidiff first)

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

    n_shuff = 5000;
    mpfc_p = nan(size(curr_animal_tds));
    for curr_day = curr_animal_tds'

        curr_day_trials = wf_grp.animal == curr_animal & ...
            wf_grp.td == curr_day & ...
            wf_grp.rxn >= rxn_cutoff;

        if sum(curr_day_trials) < 10
            continue
        end

        % trial mean, max
        curr_act = [...
            mat2cell(wf_striatum_roi(curr_day_trials,cortex_baseline_t,2,1),ones(sum(curr_day_trials),1),sum(cortex_baseline_t)), ...
            mat2cell(wf_striatum_roi(curr_day_trials,cortex_stim_t,2,1),ones(sum(curr_day_trials),1),sum(cortex_stim_t))];

        curr_mpfc_diff = diff(arrayfun(@(x) max(nanmean(vertcat(curr_act{:,x}),1)),1:2)); % max
        % curr_mpfc_diff = diff(arrayfun(@(x) nanmean(nanmean(vertcat(curr_act{:,x}),1)),1:2)); % mean

        curr_mpfc_diff_shuff = nan(1,n_shuff);
        for curr_shuff = 1:n_shuff
            curr_act_shuff = ap.shake(curr_act,2);
            curr_mpfc_diff_shuff(curr_shuff) = ...
                diff(arrayfun(@(x) max(nanmean(vertcat(curr_act_shuff{:,x}),1)),1:2)); % max
            % curr_mpfc_diff_shuff(curr_shuff) = ...
            %     diff(arrayfun(@(x) max(nanmean(vertcat(curr_act_shuff{:,x}),1)),1:2)); % mean
        end

        stat_rank = tiedrank([curr_mpfc_diff,curr_mpfc_diff_shuff]');
        stat_p = 1-stat_rank(1)/(n_shuff+1);

        mpfc_p(curr_animal_tds == curr_day) = stat_p;
    end

    mpfc_diffday{curr_animal} = (1:length(curr_animal_tds))' - min([find(mpfc_p < 0.01,1),nan]); % single day
    % mpfc_diffday{curr_animal} = (1:length(curr_animal_tds))' - min([find(movsum(max(mpfc_p,0)<0.01,[0,1])==2,1),nan]); % 2 days in a row
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


%% mPFC TASK: get day L hemi different from R hemi

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

    mpfc_diffday{curr_animal} = (1:length(curr_animal_tds))' - min([find(mpfc_p < 0.01,1),nan]);
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


%% mPFC PASSIVE: get day it's different from baseline time

%%% Load data for figure
load_dataset = 'passive';
Marica_2025.figures.load_data;
%%%

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

    n_shuff = 5000;
    mpfc_p = nan(size(curr_animal_tds));
    for curr_day = curr_animal_tds(2:end)'

        curr_day_trials = wf_grp.animal == curr_animal & ...
            wf_grp.td == curr_day & ...
            wf_grp.stim == 90;

        if sum(curr_day_trials) < 10
            continue
        end

        % trial mean, max
        curr_act = [...
            mat2cell(wf_striatum_roi(curr_day_trials,cortex_baseline_t,2,1),ones(sum(curr_day_trials),1),sum(cortex_baseline_t)), ...
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


%% mPFC PASSIVE: get day it's different from naive or td1

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


%% mPFC PASSIVE: get day R stim different from L stim

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

%% (testing: plot psth for individual animal)

spikes_norm_smooth_fcn = @(spikes,baseline) ...
    smoothdata((spikes-baseline)./(baseline+softnorm),2, ...
    'gaussian',[100,0]);

striatum_mua_rec = cellfun(spikes_norm_smooth_fcn, ...
    striatum_mua_sum,mua_baseline,'uni',false);


striatum_mua_max_rec = cellfun(@(data,domain,trials) ...
    nanmean(data(trials,:,domain==1,1),1), ...
    striatum_mua_rec,striatum_domain_idx,use_trials, ...
    'ErrorHandler',@(varargin) NaN(1,1501),'uni',false);

striatum_mua_max_rec{25} = nan(1,1501);
x = cell2mat(striatum_mua_max_rec);

figure;plot(psth_t,x(strcmp(bhv.animal,'AM022'),:)');


%% Activity v performance (task or passive)

% Get performance as reaction index
rxn_meas = cellfun(@nanmean,bhv.stim_to_move);
rxn_null = cellfun(@nanmean,bhv.stim_to_move_nullmean);

rxn_idx = (rxn_null-rxn_meas)./(rxn_null+rxn_meas);

% rxn_idx = cellfun(@nanmedian,bhv.stim_to_move);
% rxn_idx = bhv.days_from_learning;

n_learned_day = cellfun(@(x) max([0, ...
    find(bhv.learned_days(strcmp(bhv.animal,x)),1)]), ...
    unique(bhv.animal,'stable'));

figure; h = tiledlayout('flow');
for curr_animal = unique(bhv.animal)'
    nexttile;
    plot(rxn_idx(strcmp(bhv.animal,curr_animal)),'linewidth',2);
    xline(n_learned_day(find(strcmp(curr_animal,unique(bhv.animal)),1)),'r');
end
linkaxes(get(h,'Children'),'y');

figure; hold on
for curr_animal = unique(bhv.animal)'
    plot(bhv.days_from_learning(strcmp(bhv.animal,curr_animal)), ...
        rxn_idx(strcmp(bhv.animal,curr_animal)));
end


% Get activity

% Average day, get max in time window
stim_t = [0,0.2];
cortex_stim_t = isbetween(wf_t,stim_t(1),stim_t(2));
striatum_stim_t = isbetween(psth_t,stim_t(1),stim_t(2));

switch load_dataset
    case 'task'
        rxn_cutoff = 0.3;
        use_trials = cellfun(@(x) x > rxn_cutoff,bhv.stim_to_move,'uni',false);
    case 'passive'
        use_trials = cellfun(@(x) x == 90,wf.trial_stim_values,'uni',false);
end

wf_striatum_roi_max_rec = cell2mat(cellfun(@(data,trials) ...
    squeeze(max(nanmean(data(trials,cortex_stim_t,:),1),[],2))', ...
    mat2cell(wf_striatum_roi(:,:,:,1),cortex_trials_rec_n,length(wf_t),3), ...
    use_trials,'uni',false));

spikes_norm_smooth_fcn = @(spikes,sub_baseline,div_baseline) ...
    smoothdata((spikes-sub_baseline)./(div_baseline+softnorm),2, ...
    'gaussian',[100,0]);
striatum_mua_rec = cellfun(spikes_norm_smooth_fcn, ...
    striatum_mua_sum,mua_sub_baseline,mua_div_baseline,'uni',false);

striatum_mua_max_rec = cellfun(@(data,domain,trials) ...
    max([nanmean(data(trials,striatum_stim_t,domain==1,1),1),NaN],[],2), ...
    striatum_mua_rec,striatum_domain_idx,use_trials, ...
    'ErrorHandler',@(varargin) NaN);


figure; hold on
yyaxis left; plot(rxn_idx,striatum_mua_max_rec,'.','MarkerSize',15);
ylabel('Str1');
yyaxis right; plot(rxn_idx,wf_striatum_roi_max_rec(:,2),'.','MarkerSize',15)
ylabel('mPFC');
xlabel('Performance')


n_rec = cellfun(@(x) sum(strcmp(x,bhv.animal)),unique(bhv.animal,'stable'));
figure;
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'))
cellfun(@(x,y) plot(x(~isnan(y)),y(~isnan(y)),'.-','MarkerSize',15),mat2cell(rxn_idx,n_rec),mat2cell(striatum_mua_max_rec,n_rec));
plot(ap.groupfun(@nanmean,rxn_idx,bhv.days_from_learning),ap.groupfun(@nanmean,striatum_mua_max_rec,bhv.days_from_learning),'k','linewidth',2)
ylabel('Str1');
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'))
cellfun(@(x,y) plot(x(~isnan(y)),y(~isnan(y)),'.-','MarkerSize',15),mat2cell(rxn_idx,n_rec),mat2cell(wf_striatum_roi_max_rec(:,2),n_rec));
plot(ap.groupfun(@nanmean,rxn_idx,bhv.days_from_learning),ap.groupfun(@nanmean,wf_striatum_roi_max_rec(:,2),bhv.days_from_learning),'k','linewidth',2)
ylabel('mPFC');




figure; hold on; set(gca,'ColorOrder',ap.colormap('tube'));
cellfun(@(x,y) plot(x(~isnan(x)),y(~isnan(x)),'.','MarkerSize',15), ...
    mat2cell(striatum_mua_max_rec,n_rec),mat2cell(wf_striatum_roi_max_rec(:,2),n_rec));
xlabel('Striatum');
ylabel('Cortex');

figure; hold on; set(gca,'ColorOrder',ap.colormap('tube'));
cellfun(@(x,y,z) plot3(x(~isnan(y)),y(~isnan(y)),z(~isnan(y)),'.','MarkerSize',15), ...
    mat2cell(rxn_idx,n_rec),mat2cell(striatum_mua_max_rec,n_rec),mat2cell(wf_striatum_roi_max_rec(:,2),n_rec));
xlabel('Performance');
ylabel('Striatum');
zlabel('Cortex');
view([45,45]); axis vis3d; grid on;


n_rec = cellfun(@(x) sum(strcmp(x,bhv.animal)),unique(bhv.animal,'stable'));
figure;
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'))
cellfun(@(x,y) plot(y,'.-','MarkerSize',15),mat2cell(rxn_idx,n_rec),mat2cell(striatum_mua_max_rec,n_rec));
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'))
cellfun(@(x,y) plot(y,'.-','MarkerSize',15),mat2cell(rxn_idx,n_rec),mat2cell(wf_striatum_roi_max_rec(:,2),n_rec));


str_padcat = AP_padcatcell(mat2cell(striatum_mua_max_rec,n_rec));
mpfc_padcat = AP_padcatcell(mat2cell(wf_striatum_roi_max_rec(:,2),n_rec));
figure;
yyaxis left; hold on;
plot(str_padcat)
plot(nanmean(str_padcat,2),'linewidth',3);
yyaxis right; hold on;
plot(mpfc_padcat)
plot(nanmean(mpfc_padcat,2),'linewidth',3);


%% Widefield kernels

%%% Load non-activity data
load_dataset = 'noact';
Marica_2025.figures.load_data;
%%%

% Load kernels
U_master = plab.wf.load_master_U;
load(fullfile(data_path,'wf_kernels'));

% Plot grand average
n_vs = 200;
task_kernel_avg = plab.wf.svd2px(U_master(:,:,1:n_vs),nanmean(cat(4,wf_kernels.task_kernels{:}),4));
passive_kernel_avg = plab.wf.svd2px(U_master(:,:,1:n_vs),nanmean(cat(4,wf_kernels.passive_kernels{:}),4));

ap.imscroll(task_kernel_avg);axis image
colormap(gca,ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);

ap.imscroll(passive_kernel_avg);axis image
colormap(gca,ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);

% CHANGE THIS: average days within animals > average animals
% Plot pre/post learning
wf_kernel_task_prepost = plab.wf.svd2px(U_master(:,:,1:n_vs),ap.groupfun(@nanmean,cat(4,wf_kernels.task_kernels{:}),[],[],[],bhv.days_from_learning>=0));
ap.imscroll(squeeze(wf_kernel_task_prepost));
axis image
colormap(gca,ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);

wf_kernel_passive_prepost = plab.wf.svd2px(U_master(:,:,1:n_vs),ap.groupfun(@nanmean,cat(4,wf_kernels.passive_kernels{:}),[],[],[],bhv.days_from_learning>=0));
ap.imscroll(squeeze(wf_kernel_passive_prepost(:,:,:,3,:)));
axis image
colormap(gca,ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);

% Plot max kernels in task and passive
wf_kernel_max = cat(3,squeeze(max(wf_kernel_task_prepost,[],3)), ...
    squeeze(max(wf_kernel_passive_prepost(:,:,:,3,:),[],3)));
figure; tiledlayout(2,2,'TileSpacing','none');
c = [0,3e-4];
for curr_im = 1:size(wf_kernel_max,3)
    nexttile;
    imagesc(wf_kernel_max(:,:,curr_im));
    ap.wf_draw('ccf',[0.5,0.5,0.5]);
    colormap(gca,ap.colormap('WG'));
    axis image off;
end


% Plot by LD
ld_x = unique(bhv.days_from_learning(~isnan(bhv.days_from_learning)));
% r = squeeze(max(plab.wf.svd2px(U_master(:,:,1:n_vs),ap.groupfun(@nanmean,cat(4,wf_kernels.task_kernels{:}),[],[],[],bhv.days_from_learning)),[],3));
% r = squeeze(max(plab.wf.svd2px(U_master(:,:,1:n_vs),ap.groupfun(@nanmean,cat(4,wf_kernels.passive_kernels{:}),[],[],[],bhv.days_from_learning)),[],3));
ap.imscroll(r);
axis image;
colormap(gca,ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);


n_rec = cellfun(@(x) sum(strcmp(x,bhv.animal)),unique(bhv.animal,'stable'));
m_t = cell2mat(cellfun(@(x) ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi(:,:,2)),wf_kernels.task_kernels,'uni',false));
m_p = cell2mat(cellfun(@(x) ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi(:,:,2)),cellfun(@(x) x(:,:,3),wf_kernels.passive_kernels,'uni',false),'uni',false));

m_t_max = mat2cell(max(m_t,[],2),n_rec);
m_p_max = mat2cell(max(m_p,[],2),n_rec);

ld_rec = mat2cell(bhv.days_from_learning,n_rec);

figure;
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'));
cellfun(@(x,y) plot(x,y,'.-','MarkerSize',15),ld_rec,m_t_max);
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'));
cellfun(@(x,y) plot(x,y,'.-','MarkerSize',15),ld_rec,m_p_max);

figure;
hold on; set(gca,'ColorOrder',ap.colormap('tube'));
cellfun(@(x,y) plot(x,y,'.','MarkerSize',15),m_p_max,m_t_max);

rxn_meas = cellfun(@nanmean,bhv.stim_to_move);
rxn_null = cellfun(@nanmean,bhv.stim_to_move_nullmean);
rxn_idx = (rxn_null-rxn_meas)./(rxn_null+rxn_meas);
rxn_idx_rec = mat2cell(rxn_idx,n_rec);
figure;
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'));
cellfun(@(x,y) plot(x,y,'.','MarkerSize',15),rxn_idx_rec,m_t_max);
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'));
cellfun(@(x,y) plot(x,y,'.','MarkerSize',15),rxn_idx_rec,m_p_max);

figure;
hold on; set(gca,'ColorOrder',ap.colormap('tube'));
cellfun(@(p,x,y) plot3(p(~isnan(y)),x(~isnan(y)),y(~isnan(y)),'.-','MarkerSize',15),rxn_idx_rec,m_p_max,mat2cell(striatum_mua_max_rec,n_rec));
view([45,45]); axis vis3d; grid on;
xlabel('Performance');ylabel('mPFC');zlabel('Striatum');


switch load_dataset
    case 'passive'
        plot_k = m_p_max;
    case 'task'
        plot_k = m_t_max;
end

figure;
nexttile;hold on; set(gca,'ColorOrder',ap.colormap('tube'));
cellfun(@(x,y) plot(x,y,'.','MarkerSize',15),plot_k,mat2cell(striatum_mua_max_rec,n_rec));
plot(ap.groupfun(@nanmean,cell2mat(plot_k),bhv.days_from_learning), ...
    ap.groupfun(@nanmean,striatum_mua_max_rec,bhv.days_from_learning),'k','linewidth',2);
xlabel('mPFC - kernel')

nexttile;hold on; set(gca,'ColorOrder',ap.colormap('tube'));
cellfun(@(x,y) plot(x,y,'.','MarkerSize',15), ...
    mat2cell(wf_striatum_roi_max_rec(:,2),n_rec),mat2cell(striatum_mua_max_rec,n_rec));
plot(ap.groupfun(@nanmean,wf_striatum_roi_max_rec(:,2),bhv.days_from_learning), ...
    ap.groupfun(@nanmean,striatum_mua_max_rec,bhv.days_from_learning),'k','linewidth',2);
xlabel('mPFC - dff')


% Normalize task/passive by LD0 (save striatum_mua_max_rec passive in as matlab.mat on desktop)
x = load('C:\Users\petersa\Desktop\matlab.mat');
norm_val = cellfun(@(x,y) max([(x(y==0)),nan]),mat2cell(striatum_mua_max_rec,n_rec),mat2cell(bhv.days_from_learning,n_rec),'uni',false,'ErrorHandler',@(varargin) NaN);
str_task_norm = cellfun(@(x,y) x./y,mat2cell(striatum_mua_max_rec,n_rec),norm_val,'uni',false,'ErrorHandler',@(varargin) []);
str_passive_norm = cellfun(@(x,y) x./y,mat2cell(x.striatum_mua_max_rec,n_rec),norm_val,'uni',false,'ErrorHandler',@(varargin) []);

norm_val = cellfun(@(x,y) max([(x(y==0)),nan]),m_t_max,mat2cell(bhv.days_from_learning,n_rec),'uni',false,'ErrorHandler',@(varargin) NaN);
m_t_max_norm = cellfun(@(x,y) x./y,m_t_max,norm_val,'uni',false,'ErrorHandler',@(varargin) []);
m_p_max_norm = cellfun(@(x,y) x./y,m_p_max,norm_val,'uni',false,'ErrorHandler',@(varargin) []);

% % (dff)
% y = load('C:\Users\petersa\Desktop\matlab2.mat');
% norm_val = cellfun(@(x,y) max([(x(y==0)),nan]),mat2cell(wf_striatum_roi_max_rec(:,2),n_rec),mat2cell(bhv.days_from_learning,n_rec),'uni',false,'ErrorHandler',@(varargin) NaN);
% m_t_max_norm = cellfun(@(x,y) x./y,mat2cell(wf_striatum_roi_max_rec(:,2),n_rec),norm_val,'uni',false,'ErrorHandler',@(varargin) []);
% m_p_max_norm = cellfun(@(x,y) x./y,mat2cell(y.wf_striatum_roi_max_rec(:,2),n_rec),norm_val,'uni',false,'ErrorHandler',@(varargin) []);

figure;
nexttile;hold on; set(gca,'ColorOrder',ap.colormap('tube'));
plot(ap.groupfun(@nanmean,cell2mat(m_t_max_norm),bhv.days_from_learning), ...
    ap.groupfun(@nanmean,cell2mat(str_task_norm),bhv.days_from_learning),'k','linewidth',2);
plot(ap.groupfun(@nanmean,cell2mat(m_p_max_norm),bhv.days_from_learning), ...
    ap.groupfun(@nanmean,cell2mat(str_passive_norm),bhv.days_from_learning),'r','linewidth',2);
xlabel('mPFC - kernel')



figure; tiledlayout(2,2);
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'))
cellfun(@(x,y) plot(x(~isnan(y)),y(~isnan(y)),'.-','MarkerSize',15),mat2cell(bhv.days_from_learning,n_rec),str_task_norm);
plot(ap.groupfun(@nanmean,rxn_idx,bhv.days_from_learning),ap.groupfun(@nanmean,cell2mat(str_task_norm),bhv.days_from_learning),'k','linewidth',2)
title('task');
ylabel('Str1');
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'))
cellfun(@(x,y) plot(x(~isnan(y)),y(~isnan(y)),'.-','MarkerSize',15),mat2cell(bhv.days_from_learning,n_rec),str_passive_norm);
plot(ap.groupfun(@nanmean,rxn_idx,bhv.days_from_learning),ap.groupfun(@nanmean,cell2mat(str_passive_norm),bhv.days_from_learning),'k','linewidth',2)
ylabel('Str1');
title('passive');
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'))
cellfun(@(x,y) plot(x(~isnan(y)),y(~isnan(y)),'.-','MarkerSize',15),mat2cell(bhv.days_from_learning,n_rec),m_t_max_norm);
plot(ap.groupfun(@nanmean,rxn_idx,bhv.days_from_learning),ap.groupfun(@nanmean,cell2mat(m_t_max_norm),bhv.days_from_learning),'k','linewidth',2)
ylabel('mPFC');
title('task');
nexttile; hold on; set(gca,'ColorOrder',ap.colormap('tube'))
cellfun(@(x,y) plot(x(~isnan(y)),y(~isnan(y)),'.-','MarkerSize',15),mat2cell(bhv.days_from_learning,n_rec),m_p_max_norm);
plot(ap.groupfun(@nanmean,rxn_idx,bhv.days_from_learning),ap.groupfun(@nanmean,cell2mat(m_p_max_norm),bhv.days_from_learning),'k','linewidth',2)
ylabel('mPFC');
title('passive');

%% Plot striatum vs mPFC kernel in task/passive

load_dataset_retain = true;

% ~~ Set up data structure params
load_dataset = 'passive';
Marica_2025.figures.load_data;

% Set up parameters for activity grids [animal x day x domain x stim]
data_grid_params = struct;

data_grid_params.stim_t = [0,0.2];
data_grid_params.cortex_stim_t = isbetween(wf_t,data_grid_params.stim_t(1),data_grid_params.stim_t(2));
data_grid_params.striatum_stim_t = isbetween(psth_t,data_grid_params.stim_t(1),data_grid_params.stim_t(2));
data_grid_params.rxn_cutoff = 0.3;

data_grid_params.ld_unique = unique((bhv.days_from_learning(~isnan(bhv.days_from_learning))));
data_grid_params.grid_size = [length(unique(bhv.animal)),length(data_grid_params.ld_unique),n_domains];

% Set up data structure
data_grids = struct;

clearvars -except load_dataset_retain data_grid_params data_grids

% ~~ Get task activity
load_dataset = 'task';
Marica_2025.figures.load_data;

% (striatum)
striatum_use_trials = striatum_mua_grp.rxn > data_grid_params.rxn_cutoff;
[~,striatum_ld_idx] = ismember(striatum_mua_grp.ld,data_grid_params.ld_unique);

[striatum_rec,striatum_rec_grp] = ap.groupfun(@nanmean,striatum_mua(:,:,1), ...
    [striatum_mua_grp.animal,striatum_ld_idx,striatum_mua_grp.domain_idx].* ...
    ap.nanout(~(striatum_use_trials & striatum_ld_idx)));
striatum_rec_tmax = max(striatum_rec(:,data_grid_params.striatum_stim_t),[],2);

data_grids.striatum_task = accumarray(striatum_rec_grp,striatum_rec_tmax,data_grid_params.grid_size,[],NaN);

% (widefield)
wf_use_trials = wf_grp.rxn > data_grid_params.rxn_cutoff;
[~,wf_ld_idx] = ismember(wf_grp.ld,data_grid_params.ld_unique);

[wf_roi_rec,wf_roi_rec_grp] = ap.groupfun(@nanmean,wf_striatum_roi(:,:,:,1), ...
    [wf_grp.animal,wf_ld_idx].* ...
    ap.nanout(~(wf_use_trials & wf_ld_idx)));
wf_roi_rec_tmax = permute(max(wf_roi_rec(:,data_grid_params.cortex_stim_t,:),[],2),[1,3,2]);

data_grids.wf_roi_task = cell2mat(permute(arrayfun(@(domain) accumarray(wf_roi_rec_grp, ...
    wf_roi_rec_tmax(:,domain),data_grid_params.grid_size(1:2),[],NaN('single')),1:n_domains,'uni',false),[1,3,2]));

clearvars -except load_dataset_retain data_grid_params data_grids

% Get passive activity
load_dataset = 'passive';
Marica_2025.figures.load_data;

% (striatum)
[~,striatum_ld_idx] = ismember(striatum_mua_grp.ld,data_grid_params.ld_unique);
[~,striatum_stim_idx] = ismember(striatum_mua_grp.stim,unique(striatum_mua_grp.stim));

[striatum_rec,striatum_rec_grp] = ap.groupfun(@nanmean,striatum_mua(:,:,1), ...
    [striatum_mua_grp.animal,striatum_ld_idx,striatum_mua_grp.domain_idx,striatum_stim_idx].* ...
    ap.nanout(~(striatum_ld_idx)));
striatum_rec_tmax = max(striatum_rec(:,data_grid_params.striatum_stim_t),[],2);

data_grids.striatum_passive = accumarray(striatum_rec_grp,striatum_rec_tmax,[data_grid_params.grid_size,max(striatum_stim_idx)],[],NaN);

% (widefield)
[~,wf_ld_idx] = ismember(wf_grp.ld,data_grid_params.ld_unique);
[~,wf_stim_idx] = ismember(wf_grp.stim,unique(wf_grp.stim));

[wf_roi_rec,wf_roi_rec_grp] = ap.groupfun(@nanmean,wf_striatum_roi(:,:,:,1), ...
    [wf_grp.animal,wf_ld_idx,wf_stim_idx].* ...
    ap.nanout(~(wf_ld_idx)));
wf_roi_rec_tmax = permute(max(wf_roi_rec(:,data_grid_params.cortex_stim_t,:),[],2),[1,3,2]);

data_grids.wf_roi_passive = permute(cell2mat(permute(arrayfun(@(domain) accumarray(wf_roi_rec_grp, ...
    wf_roi_rec_tmax(:,domain),[data_grid_params.grid_size(1:2),length(unique(wf_grp.stim))], ...
    [],NaN('single')),1:n_domains,'uni',false),[1,3,4,2])),[1,2,4,3]);

clearvars -except load_dataset_retain data_grid_params data_grids

% ~~ Get widefield stim kernels
load_dataset = 'noact';
Marica_2025.figures.load_data;

U_master = plab.wf.load_master_U;
load(fullfile(data_path,'wf_kernels'));

n_vs = size(wf_kernels.task_kernels{1},1);

wf_grid_idx = [grp2idx(bhv.animal),grp2idx(bhv.days_from_learning)];
wf_grid_idx_use = ~any(isnan(wf_grid_idx),2);

% (task)
wf_kernel_roi_task = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi),[3,2,1]), ...
    wf_kernels.task_kernels,'uni',false));

wf_kernel_roi_task_tmax = permute(max(wf_kernel_roi_task,[],2),[1,3,2]);

data_grids.wf_kernel_roi_task = ...
    cell2mat(permute(arrayfun(@(domain) accumarray(wf_grid_idx(wf_grid_idx_use,:), ...
    wf_kernel_roi_task_tmax(wf_grid_idx_use,domain),data_grid_params.grid_size(1:2),[],NaN('single')), ...
    1:n_domains,'uni',false),[1,3,2]));

% (passive)
wf_kernel_roi_passive = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi),[4,2,1,3]), ...
    wf_kernels.passive_kernels,'uni',false));
wf_kernel_roi_passive_tmax = permute(max(wf_kernel_roi_passive,[],2),[1,3,4,2]);

data_grids.wf_kernel_roi_passive = cell2mat(permute(arrayfun(@(stim) ...
    cell2mat(permute(arrayfun(@(domain) accumarray(wf_grid_idx(wf_grid_idx_use,:), ...
    wf_kernel_roi_passive_tmax(wf_grid_idx_use,domain,stim),data_grid_params.grid_size(1:2),[],NaN('single')), ...
    1:n_domains,'uni',false),[1,3,2])),1:size(wf_kernel_roi_passive_tmax,3),'uni',false), [1,3,4,2]));

clearvars -except load_dataset_retain data_grid_params data_grids


% ~~ Plot

% Normalize to task LD0
str_normval = data_grids.striatum_task(:,data_grid_params.ld_unique==0,:);
wf_normval = data_grids.wf_roi_task(:,data_grid_params.ld_unique==0,:);
wf_kernel_normval = data_grids.wf_kernel_roi_task(:,data_grid_params.ld_unique==0,:);

% Set days to plot (n>3)
plot_ld_idx = sum(~isnan(data_grids.striatum_task(:,:,1))) > 3;

% Plot task vs passive for mPFC and striatum
max_ld = max(abs(data_grid_params.ld_unique(plot_ld_idx)));
ld_colors = ap.colormap('BKR',max_ld*2+1);
plot_ld_colors = ld_colors(ismember(-max_ld:max_ld, ...
    data_grid_params.ld_unique(plot_ld_idx)),:);

plot_str = 1;
plot_wf_roi = 2;
plot_stim = 3;

figure; tiledlayout(1,2);
nexttile; hold on;
scatter(reshape(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,plot_wf_roi,plot_stim),[],1), ...
    reshape(data_grids.wf_kernel_roi_task(:,plot_ld_idx,plot_wf_roi),[],1), ...
    20,repelem(plot_ld_colors,size(data_grids.striatum_task,1),1),'filled');
plot(nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,plot_wf_roi,plot_stim),1), ...
    nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,plot_wf_roi),1),'color',[0.5,0.5,0.5],'linewidth',2);
scatter(nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,plot_wf_roi,plot_stim),1), ...
    nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,plot_wf_roi),1),80,plot_ld_colors,'linewidth',4);
xlim(prctile([xlim,ylim],[0,100]));ylim(xlim)
axis square
line(ylim,ylim,'linestyle','--','color',[0.5,0.5,0.5]);
xlabel('Passive');ylabel('Task');
title(sprintf('Cortex kernel %d',plot_wf_roi));

nexttile; hold on;
scatter(reshape(data_grids.striatum_passive(:,plot_ld_idx,plot_str,plot_stim),[],1), ...
    reshape(data_grids.striatum_task(:,plot_ld_idx,plot_str),[],1), ...
    20,repelem(plot_ld_colors,size(data_grids.striatum_task,1),1),'filled');
plot(nanmean(data_grids.striatum_passive(:,plot_ld_idx,plot_str,plot_stim),1), ...
    nanmean(data_grids.striatum_task(:,plot_ld_idx,plot_str),1),'color',[0.5,0.5,0.5],'linewidth',2);
scatter(nanmean(data_grids.striatum_passive(:,plot_ld_idx,plot_str,plot_stim),1), ...
    nanmean(data_grids.striatum_task(:,plot_ld_idx,plot_str),1),80,plot_ld_colors,'linewidth',4);
xlim(prctile([xlim,ylim],[0,100]));ylim(xlim)
axis square
line(ylim,ylim,'linestyle','--','color',[0.5,0.5,0.5]);
xlabel('Passive');ylabel('Task');
title(sprintf('Striatum %d',plot_str));
ap.prettyfig;

% Plot striatum vs mPFC for task and passive
outline_ld_cols = [0.5,0.5,0.5;0,0,0];
outline_ld = @(h) scatter(h.XData(ismember(data_grid_params.ld_unique(plot_ld_idx),[-1,0])), ...
    h.YData(ismember(data_grid_params.ld_unique(plot_ld_idx),[-1,0])), ...
    100,outline_ld_cols,'linewidth',3);

plot_task = @(roi,str,col) ...
    plot(nanmean(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
    nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,roi)./wf_kernel_normval(:,:,roi),1), ...
    '.-','color',col,'linewidth',2,'MarkerSize',30);

plot_passive = @(roi,str,plot_stim,col) ...
    plot(nanmean(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
    nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_kernel_normval(:,:,roi),1), ...
    '.-','color',col,'linewidth',2,'MarkerSize',30);

plot_str = 1;
plot_wf_roi = 2;

figure; hold on;
h = plot_task(plot_wf_roi,plot_str,[0.5,0,0]); outline_ld(h);
stim_col = [0.8,0.8,0.3;0.3,0.3,0.8;0.8,0.3,0.3];
for curr_stim = 1:3
    h = plot_passive(plot_wf_roi,plot_str,curr_stim,stim_col(curr_stim,:)); outline_ld(h);
end
xlabel(sprintf('Striatum %d (LD0-norm)',plot_str));
ylabel(sprintf('Cortex kernel %d (LD0-norm)',plot_wf_roi));
ap.prettyfig;

% ~~~ STATS ~~~
n_shuff = 10000;
sig_flag = @(p) discretize(p < 0.05,[0,1,Inf],["","*"]);

stat_day = -1; % day to compare vs mean(<day)
plot_str = 1;
plot_wf_roi = 2;
plot_stim = 3;

pre_days =  data_grid_params.ld_unique < -1;
post_days = ismember(data_grid_params.ld_unique,-1);

stat_label = {...
    sprintf('striatum %d task',plot_str), ...
    sprintf('striatum %d passive',plot_str), ...
    sprintf('cortex kernel %d task',plot_wf_roi), ...
    sprintf('cortex_kernel %d passive',plot_wf_roi)};
pre_data = [...
    nanmean(data_grids.striatum_task(:,pre_days,plot_str)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.striatum_passive(:,pre_days,plot_str,plot_stim)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.wf_kernel_roi_task(:,pre_days,plot_wf_roi)./wf_kernel_normval(:,:,plot_wf_roi),2), ...
    nanmean(data_grids.wf_kernel_roi_passive(:,pre_days,plot_wf_roi,plot_stim)./wf_kernel_normval(:,:,plot_wf_roi),2)];
post_data = [...
    nanmean(data_grids.striatum_task(:,post_days,plot_str)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.striatum_passive(:,post_days,plot_str,plot_stim)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.wf_kernel_roi_task(:,post_days,plot_wf_roi)./wf_kernel_normval(:,:,plot_wf_roi),2), ...
    nanmean(data_grids.wf_kernel_roi_passive(:,post_days,plot_wf_roi,plot_stim)./wf_kernel_normval(:,:,plot_wf_roi),2)];

data_meas = nanmean(diff(cat(3,pre_data,post_data),[],3),1);
data_shuff = nan(n_shuff,size(data_meas,2));
for curr_shuff = 1:n_shuff
    data_shuff(curr_shuff,:) = nanmean(diff(ap.shake(cat(3,pre_data,post_data),3),[],3),1);
end
stat_rank = tiedrank([data_meas;data_shuff]);
stat_p = 1-stat_rank(1,:)/(n_shuff+1);
for curr_stat = 1:length(stat_p)
    fprintf('Shuffle day %d vs mean(<): %s p = %.2g%s\n', ...
        stat_day,stat_label{curr_stat},stat_p(curr_stat),sig_flag(stat_p(curr_stat)));
end



%%%%%%% TESTING

% plot_task = @(roi,str,col) ...
%     plot(nanmean(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
%     nanmean(data_grids.wf_roi_task(:,plot_ld_idx,roi)./wf_normval(:,:,roi),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30);

plot_task = @(roi,str,col) ...
    plot(nanmean(data_grids.striatum_task(:,plot_ld_idx,str),1), ...
    nanmean(data_grids.wf_roi_task(:,plot_ld_idx,roi),1), ...
    '.-','color',col,'linewidth',2,'MarkerSize',30);

figure; hold on;
h = plot_task(2,1,'k'); outline_ld(h);
h = plot_task(2,2,'r'); outline_ld(h);


%% --> ALT VERSION: LD bins

load_dataset_retain = true;

% ~~ Set up data structure params
load_dataset = 'passive';
Marica_2025.figures.load_data;

% Set up parameters for activity grids [animal x day x domain x stim]
data_grid_params = struct;

data_grid_params.stim_t = [0,0.2];
data_grid_params.cortex_stim_t = isbetween(wf_t,data_grid_params.stim_t(1),data_grid_params.stim_t(2));
data_grid_params.striatum_stim_t = isbetween(psth_t,data_grid_params.stim_t(1),data_grid_params.stim_t(2));
data_grid_params.rxn_cutoff = 0.3;

data_grid_params.ld_bins = [-Inf,-2:0,1,Inf];
data_grid_params.ld_unique = 1:(length(data_grid_params.ld_bins)-1);
data_grid_params.grid_size = [length(unique(bhv.animal)),length(data_grid_params.ld_unique),n_domains];

% Set up data structure
data_grids = struct;

clearvars -except load_dataset_retain data_grid_params data_grids

% ~~ Get task activity
load_dataset = 'task';
Marica_2025.figures.load_data;

% (striatum)
striatum_use_trials = striatum_mua_grp.rxn > data_grid_params.rxn_cutoff;
striatum_ld_idx = discretize(striatum_mua_grp.ld,data_grid_params.ld_bins);

[striatum_rec,striatum_rec_grp] = ap.groupfun(@nanmean,striatum_mua(:,:,1), ...
    [striatum_mua_grp.animal,striatum_ld_idx,striatum_mua_grp.domain_idx].* ...
    ap.nanout(~(striatum_use_trials & ~isnan(striatum_ld_idx))));
striatum_rec_tmax = max(striatum_rec(:,data_grid_params.striatum_stim_t),[],2);

data_grids.striatum_task = accumarray(striatum_rec_grp,striatum_rec_tmax,data_grid_params.grid_size,[],NaN);

% (widefield)
wf_use_trials = wf_grp.rxn > data_grid_params.rxn_cutoff;
wf_ld_idx = discretize(wf_grp.ld,data_grid_params.ld_bins);

[wf_roi_rec,wf_roi_rec_grp] = ap.groupfun(@nanmean,wf_striatum_roi(:,:,:,1), ...
    [wf_grp.animal,wf_ld_idx].* ...
    ap.nanout(~(wf_use_trials & ~isnan(wf_ld_idx))));
wf_roi_rec_tmax = permute(max(wf_roi_rec(:,data_grid_params.cortex_stim_t,:),[],2),[1,3,2]);

data_grids.wf_roi_task = cell2mat(permute(arrayfun(@(domain) accumarray(wf_roi_rec_grp, ...
    wf_roi_rec_tmax(:,domain),data_grid_params.grid_size(1:2),[],NaN('single')),1:n_domains,'uni',false),[1,3,2]));

clearvars -except load_dataset_retain data_grid_params data_grids

% Get passive activity
load_dataset = 'passive';
Marica_2025.figures.load_data;

% (striatum)
striatum_ld_idx = discretize(striatum_mua_grp.ld,data_grid_params.ld_bins);
[~,striatum_stim_idx] = ismember(striatum_mua_grp.stim,unique(striatum_mua_grp.stim));

[striatum_rec,striatum_rec_grp] = ap.groupfun(@nanmean,striatum_mua(:,:,1), ...
    [striatum_mua_grp.animal,striatum_ld_idx,striatum_mua_grp.domain_idx,striatum_stim_idx].* ...
    ap.nanout(isnan(striatum_ld_idx)));
striatum_rec_tmax = max(striatum_rec(:,data_grid_params.striatum_stim_t),[],2);

data_grids.striatum_passive = accumarray(striatum_rec_grp,striatum_rec_tmax,[data_grid_params.grid_size,max(striatum_stim_idx)],[],NaN);

% (widefield)
wf_ld_idx = discretize(wf_grp.ld,data_grid_params.ld_bins);
[~,wf_stim_idx] = ismember(wf_grp.stim,unique(wf_grp.stim));

[wf_roi_rec,wf_roi_rec_grp] = ap.groupfun(@nanmean,wf_striatum_roi(:,:,:,1), ...
    [wf_grp.animal,wf_ld_idx,wf_stim_idx].* ...
    ap.nanout(isnan(wf_ld_idx)));
wf_roi_rec_tmax = permute(max(wf_roi_rec(:,data_grid_params.cortex_stim_t,:),[],2),[1,3,2]);

data_grids.wf_roi_passive = permute(cell2mat(permute(arrayfun(@(domain) accumarray(wf_roi_rec_grp, ...
    wf_roi_rec_tmax(:,domain),[data_grid_params.grid_size(1:2),length(unique(wf_grp.stim))], ...
    [],NaN('single')),1:n_domains,'uni',false),[1,3,4,2])),[1,2,4,3]);

clearvars -except load_dataset_retain data_grid_params data_grids

% ~~ Get widefield stim kernels
load_dataset = 'noact';
Marica_2025.figures.load_data;

U_master = plab.wf.load_master_U;
load(fullfile(data_path,'wf_kernels'));

n_vs = size(wf_kernels.task_kernels{1},1);

wf_grid_idx = [grp2idx(bhv.animal),discretize(bhv.days_from_learning,data_grid_params.ld_bins)];
wf_grid_idx_use = ~any(isnan(wf_grid_idx),2);

% (task)
wf_kernel_roi_task = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi),[3,2,1]), ...
    wf_kernels.task_kernels,'uni',false));

wf_kernel_roi_task_tmax = permute(max(wf_kernel_roi_task,[],2),[1,3,2]);

data_grids.wf_kernel_roi_task = ...
    cell2mat(permute(arrayfun(@(domain) accumarray(wf_grid_idx(wf_grid_idx_use,:), ...
    wf_kernel_roi_task_tmax(wf_grid_idx_use,domain),data_grid_params.grid_size(1:2),@nanmean,NaN('single')), ...
    1:n_domains,'uni',false),[1,3,2]));

% (passive)
wf_kernel_roi_passive = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master(:,:,1:n_vs),x,[],[],striatum_wf_roi),[4,2,1,3]), ...
    wf_kernels.passive_kernels,'uni',false));
wf_kernel_roi_passive_tmax = permute(max(wf_kernel_roi_passive,[],2),[1,3,4,2]);

data_grids.wf_kernel_roi_passive = cell2mat(permute(arrayfun(@(stim) ...
    cell2mat(permute(arrayfun(@(domain) accumarray(wf_grid_idx(wf_grid_idx_use,:), ...
    wf_kernel_roi_passive_tmax(wf_grid_idx_use,domain,stim),data_grid_params.grid_size(1:2),@nanmean,NaN('single')), ...
    1:n_domains,'uni',false),[1,3,2])),1:size(wf_kernel_roi_passive_tmax,3),'uni',false), [1,3,4,2]));

clearvars -except load_dataset_retain data_grid_params data_grids


% ~~ Plot

% Normalize to task LD0
str_normval = data_grids.striatum_task(:,find(data_grid_params.ld_bins>=0,1),:);
wf_normval = data_grids.wf_roi_task(:,find(data_grid_params.ld_bins>=0,1),:);
wf_kernel_normval = data_grids.wf_kernel_roi_task(:,find(data_grid_params.ld_bins>=0,1),:);

% Set days to plot (n>3)
% plot_ld_idx = sum(~isnan(data_grids.striatum_task(:,:,1))) > 3;
plot_ld_idx = 1:size(data_grids.striatum_task,2);

% Plot task vs passive for mPFC and striatum
max_ld = max(abs(data_grid_params.ld_unique(plot_ld_idx)));
ld_colors = ap.colormap('BKR',max_ld+(1-mod(max_ld,2)));
% plot_ld_colors = ld_colors(ismember(-max_ld:max_ld, ...
%     data_grid_params.ld_unique(plot_ld_idx)),:);
plot_ld_colors = ld_colors(1:max_ld,:);

plot_str = 2;
plot_wf_roi = 2;
plot_stim = 3;

figure; tiledlayout(1,2);
nexttile; hold on;
scatter(reshape(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,plot_wf_roi,plot_stim),[],1), ...
    reshape(data_grids.wf_kernel_roi_task(:,plot_ld_idx,plot_wf_roi),[],1), ...
    20,repelem(plot_ld_colors,size(data_grids.striatum_task,1),1),'filled');
plot(nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,plot_wf_roi,plot_stim),1), ...
    nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,plot_wf_roi),1),'color',[0.5,0.5,0.5],'linewidth',2);
scatter(nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,plot_wf_roi,plot_stim),1), ...
    nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,plot_wf_roi),1),80,plot_ld_colors,'linewidth',4);
xlim(prctile([xlim,ylim],[0,100]));ylim(xlim)
axis square
line(ylim,ylim,'linestyle','--','color',[0.5,0.5,0.5]);
xlabel('Passive');ylabel('Task');
title(sprintf('Cortex kernel %d',plot_wf_roi));

nexttile; hold on;
scatter(reshape(data_grids.striatum_passive(:,plot_ld_idx,plot_str,plot_stim),[],1), ...
    reshape(data_grids.striatum_task(:,plot_ld_idx,plot_str),[],1), ...
    20,repelem(plot_ld_colors,size(data_grids.striatum_task,1),1),'filled');
plot(nanmean(data_grids.striatum_passive(:,plot_ld_idx,plot_str,plot_stim),1), ...
    nanmean(data_grids.striatum_task(:,plot_ld_idx,plot_str),1),'color',[0.5,0.5,0.5],'linewidth',2);
scatter(nanmean(data_grids.striatum_passive(:,plot_ld_idx,plot_str,plot_stim),1), ...
    nanmean(data_grids.striatum_task(:,plot_ld_idx,plot_str),1),80,plot_ld_colors,'linewidth',4);
xlim(prctile([xlim,ylim],[0,100]));ylim(xlim)
axis square
line(ylim,ylim,'linestyle','--','color',[0.5,0.5,0.5]);
xlabel('Passive');ylabel('Task');
title(sprintf('Striatum %d',plot_str));
ap.prettyfig;

% Plot striatum vs mPFC for task and passive
outline_ld_cols = [0.5,0.5,0.5;0,0,0];
outline_ld = @(h) scatter(h.XData(find(data_grid_params.ld_bins==0)+[-1,0]), ...
    h.YData(find(data_grid_params.ld_bins==0)+[-1,0]), ...
    100,outline_ld_cols,'linewidth',3);

% % (wf kernel, LD-0 norm)
% plot_task = @(roi,str,col) ...
%     plot(nanmean(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
%     nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,roi)./wf_kernel_normval(:,:,roi),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30);
% plot_passive = @(roi,str,plot_stim,col) ...
%     plot(nanmean(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
%     nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_kernel_normval(:,:,roi),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30);

% (wf kernel, LD-0 norm,errorbar)
plot_task = @(roi,str,col) ...
    errorbar(nanmean(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
    nanmean(data_grids.wf_kernel_roi_task(:,plot_ld_idx,roi)./wf_kernel_normval(:,:,roi),1), ...
    ap.sem(data_grids.wf_kernel_roi_task(:,plot_ld_idx,roi)./wf_kernel_normval(:,:,roi),1), ...
    ap.sem(data_grids.wf_kernel_roi_task(:,plot_ld_idx,roi)./wf_kernel_normval(:,:,roi),1), ...
    ap.sem(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
    ap.sem(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
    '.-','color',col,'linewidth',2,'MarkerSize',30,'capsize',0);
plot_passive = @(roi,str,plot_stim,col) ...
    errorbar(nanmean(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
    nanmean(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_kernel_normval(:,:,roi),1), ...
    ap.sem(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_kernel_normval(:,:,roi),1), ...
    ap.sem(data_grids.wf_kernel_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_kernel_normval(:,:,roi),1), ...
    ap.sem(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
    ap.sem(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
    '.-','color',col,'linewidth',2,'MarkerSize',30,'capsize',0);

% % (wf dff, non-norm)
% plot_task = @(roi,str,col) ...
%     plot(nanmean(data_grids.striatum_task(:,plot_ld_idx,str),1), ...
%     nanmean(data_grids.wf_roi_task(:,plot_ld_idx,roi),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30);
% plot_passive = @(roi,str,plot_stim,col) ...
%     plot(nanmean(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim),1), ...
%     nanmean(data_grids.wf_roi_passive(:,plot_ld_idx,roi,plot_stim),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30);

% % (wf dff, norm)
% plot_task = @(roi,str,col) ...
%     plot(nanmean(data_grids.striatum_task(:,plot_ld_idx,str)./str_normval(:,:,str),1), ...
%     nanmean(data_grids.wf_roi_task(:,plot_ld_idx,roi)./wf_normval(:,:,roi),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30);
% plot_passive = @(roi,str,plot_stim,col) ...
%     plot(nanmean(data_grids.striatum_passive(:,plot_ld_idx,str,plot_stim)./str_normval(:,:,str),1), ...
%     nanmean(data_grids.wf_roi_passive(:,plot_ld_idx,roi,plot_stim)./wf_normval(:,:,roi),1), ...
%     '.-','color',col,'linewidth',2,'MarkerSize',30);

plot_str = 1;
plot_wf_roi = 2;

figure; hold on;
h = plot_task(plot_wf_roi,plot_str,[0.5,0,0]); outline_ld(h);
stim_col = [0.8,0.8,0.3;0.3,0.3,0.8;0.8,0.3,0.3];
for curr_stim = 1:3
    h = plot_passive(plot_wf_roi,plot_str,curr_stim,stim_col(curr_stim,:)); outline_ld(h);
end
xlabel(sprintf('Striatum %d (LD0-norm)',plot_str));
ylabel(sprintf('Cortex kernel %d (LD0-norm)',plot_wf_roi));
ap.prettyfig;


figure; hold on;
plot_str = [1,2];
plot_wf_roi = 2;
h1 = plot_task(plot_wf_roi,plot_str(1),[0.5,0,0]);
h2 = plot_task(plot_wf_roi,plot_str(2),'k');
outline_ld(h1); outline_ld(h2);
axis equal square
line(xlim,xlim,'color',[0.5,0.5,0.5]);
legend(arrayfun(@(x) sprintf('Str %d',x),plot_str,'uni',false));
xlabel(sprintf('Cortex kernel %d (LD0-norm)',plot_wf_roi));
ylabel('Striatum (LD0-norm)')
ap.prettyfig;


% ~~~ STATS ~~~
n_shuff = 10000;
sig_flag = @(p) discretize(p < 0.05,[0,1,Inf],["","*"]);

plot_str = 1;
plot_wf_roi = 2;
plot_stim = 3;

pre_days =  data_grid_params.ld_unique == 1;
post_days = ismember(data_grid_params.ld_unique,3);

stat_label = {...
    sprintf('striatum %d task',plot_str), ...
    sprintf('striatum %d passive',plot_str), ...
    sprintf('cortex kernel %d task',plot_wf_roi), ...
    sprintf('cortex_kernel %d passive',plot_wf_roi)};
pre_data = [...
    nanmean(data_grids.striatum_task(:,pre_days,plot_str)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.striatum_passive(:,pre_days,plot_str,plot_stim)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.wf_kernel_roi_task(:,pre_days,plot_wf_roi)./wf_kernel_normval(:,:,plot_wf_roi),2), ...
    nanmean(data_grids.wf_kernel_roi_passive(:,pre_days,plot_wf_roi,plot_stim)./wf_kernel_normval(:,:,plot_wf_roi),2)];
post_data = [...
    nanmean(data_grids.striatum_task(:,post_days,plot_str)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.striatum_passive(:,post_days,plot_str,plot_stim)./str_normval(:,:,plot_str),2), ...
    nanmean(data_grids.wf_kernel_roi_task(:,post_days,plot_wf_roi)./wf_kernel_normval(:,:,plot_wf_roi),2), ...
    nanmean(data_grids.wf_kernel_roi_passive(:,post_days,plot_wf_roi,plot_stim)./wf_kernel_normval(:,:,plot_wf_roi),2)];

data_meas = nanmean(diff(cat(3,pre_data,post_data),[],3),1);
data_shuff = nan(n_shuff,size(data_meas,2));
for curr_shuff = 1:n_shuff
    data_shuff(curr_shuff,:) = nanmean(diff(ap.shake(cat(3,pre_data,post_data),3),[],3),1);
end
stat_rank = tiedrank([data_meas;data_shuff]);
stat_p = 1-stat_rank(1,:)/(n_shuff+1);
for curr_stat = 1:length(stat_p)
    fprintf('Shuffle: %s p = %.2g%s\n', ...
        stat_label{curr_stat},stat_p(curr_stat),sig_flag(stat_p(curr_stat)));
end

for curr_stat = 1:size(pre_data,2)
    stat_signrank = ranksum(pre_data(:,curr_stat),post_data(:,curr_stat));
    fprintf('Ranksum: %s p = %.2g%s\n', ...
        stat_label{curr_stat},stat_signrank,sig_flag(stat_signrank));
end


%% R1: Example animal

%%% Load data for figure
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

% %%% Load data for figure
% load_dataset = 'passive';
% Marica_2025.figures.load_data;
% %%%


% Get PSTH by recording
rxn_cutoff = 0.3;
% ld_unique = unique((bhv.days_from_learning(~isnan(bhv.days_from_learning))));

% plot_day_bins = [-Inf,-2,-1,0,1,Inf];
plot_day_bins = [-Inf,-2,0,2];
cortex_plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);
striatum_plot_day_grp = discretize(max(striatum_mua_grp.ld,-inf),plot_day_bins);

% (task)
striatum_use_trials = striatum_mua_grp.rxn > rxn_cutoff;
wf_use_trials = wf_grp.rxn > rxn_cutoff;

% (passive)
% striatum_use_trials = striatum_mua_grp.stim == 90;
% wf_use_trials = wf_grp.stim == 90;

[striatum_rec,striatum_rec_grp] = ap.groupfun(@nanmean,striatum_mua(:,:,1), ...
    [striatum_mua_grp.animal,striatum_plot_day_grp,striatum_mua_grp.domain_idx].* ...
    ap.nanout(~(striatum_use_trials)));

[wf_roi_rec,wf_roi_rec_grp] = ap.groupfun(@nanmean,wf_striatum_roi(:,:,:,1), ...
    [wf_grp.animal,cortex_plot_day_grp].*ap.nanout(~(wf_use_trials)));

% Get grid of animal/day/domain recording presence
grid_size = [length(unique(bhv.animal)),length(plot_day_bins)-1,n_domains];
str_present_grid = accumarray(striatum_rec_grp,1,grid_size,@sum,0);
ap.imscroll(str_present_grid);
axis on;

figure;tiledlayout(4,14,'TileSpacing','none','TileIndexing','ColumnMajor');
for plot_animal = 1:14

    plot_t = [0,0.5];
    cortex_plot_t = isbetween(wf_t,plot_t(1),plot_t(2));
    striatum_plot_t = isbetween(psth_t,plot_t(1),plot_t(2));

    plot_domain = 1;
    nexttile;
    plot(reshape(padarray(wf_roi_rec(wf_roi_rec_grp(:,1) == plot_animal, ...
        cortex_plot_t,plot_domain),[0,5],NaN,'post')',[],1));
    title(sprintf('Cortex %d',plot_domain)); axis off;
    nexttile;
    plot(reshape(padarray(striatum_rec(striatum_rec_grp(:,1) == plot_animal & ...
        striatum_rec_grp(:,3) == plot_domain,striatum_plot_t),[0,5],NaN,'post')',[],1));
    title(sprintf('Striatum %d',plot_domain)); axis off;

    plot_domain = 2;
    nexttile;
    plot(reshape(padarray(wf_roi_rec(wf_roi_rec_grp(:,1) == plot_animal, ...
        cortex_plot_t,plot_domain),[0,5],NaN,'post')',[],1));
    title(sprintf('Cortex %d',plot_domain)); axis off;
    nexttile; 
    plot(reshape(padarray(striatum_rec(striatum_rec_grp(:,1) == plot_animal & ...
        striatum_rec_grp(:,3) == plot_domain,striatum_plot_t),[0,5],NaN,'post')',[],1));
    title(sprintf('Striatum %d',plot_domain)); axis off;

end

figure; tiledlayout(4,1);
plot_animal = 12;

plot_t = [-0.5,0.5];
cortex_plot_t = isbetween(wf_t,plot_t(1),plot_t(2));
striatum_plot_t = isbetween(psth_t,plot_t(1),plot_t(2));

plot_domain = 1;
nexttile;
plot(reshape(padarray(wf_roi_rec(wf_roi_rec_grp(:,1) == plot_animal, ...
    cortex_plot_t,plot_domain),[0,5],NaN,'post')',[],1));
title(sprintf('Cortex %d',plot_domain)); axis off;
nexttile;
plot(reshape(padarray(striatum_rec(striatum_rec_grp(:,1) == plot_animal & ...
    striatum_rec_grp(:,3) == plot_domain,striatum_plot_t),[0,5],NaN,'post')',[],1));
title(sprintf('Striatum %d',plot_domain)); axis off;

plot_domain = 2;
nexttile;
plot(reshape(padarray(wf_roi_rec(wf_roi_rec_grp(:,1) == plot_animal, ...
    cortex_plot_t,plot_domain),[0,5],NaN,'post')',[],1));
title(sprintf('Cortex %d',plot_domain)); axis off;
nexttile;
plot(reshape(padarray(striatum_rec(striatum_rec_grp(:,1) == plot_animal & ...
    striatum_rec_grp(:,3) == plot_domain,striatum_plot_t),[0,5],NaN,'post')',[],1));
title(sprintf('Striatum %d',plot_domain)); axis off;


%% (testing) plot domain map for animal/day

x = cat(3,ctx_str_maps.cortex_striatum_map{:});

idx = strcmp(bhv.animal,'AM029') & bhv.days_from_learning == 0;

idx2 = cellfun(@(x) false(size(x)),domain_idx_rec,'uni',false);
idx2{idx}(domain_idx_rec{idx} == 2) = true;

k2 = x(:,:,cell2mat(idx2));

disp(sum(cell2mat(idx2)));

ap.imscroll(k2);
axis image off
clim(max(abs(clim)).*[-1,1]);
colormap(gca,ap.colormap('PWG'));


%% Rate of non-stim movements

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

use_stat = 'mean';

% Grab learning day for each mouse
bhv = struct;

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};

    use_workflow = 'stim_wheel_right_stage\d';
    recordings = plab.find_recordings(animal,[],use_workflow);

    n_bins = 3;

    move_rate = nan(length(recordings),n_bins);
    move_rate_stim = nan(length(recordings),n_bins);
    move_rate_nonstim = nan(length(recordings),n_bins);
    p_cue_given_move = nan(length(recordings),n_bins);
    rxn_stat_p = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get P(cue|wheel_move)
        move_onsets = timelite.timestamps(diff(wheel_move) == 1);
        move_onset_prev_stim = interp1(stimOn_times,stimOn_times,move_onsets,'previous','extrap');
        move_prevstim_t = move_onsets - move_onset_prev_stim;

        move_poststim = move_prevstim_t < 0.4;

        % Bin by time within session
        session_bin_edges = linspace(timelite.timestamps(1),timelite.timestamps(end),n_bins+1);
        move_session_bins = discretize(move_onsets,session_bin_edges);

        move_rate(curr_recording,:) = ap.groupfun(@length,+move_poststim,move_session_bins)./diff(session_bin_edges)';
        move_rate_stim(curr_recording,:) = ap.groupfun(@sum,+move_poststim,move_session_bins)./diff(session_bin_edges)';
        move_rate_nonstim(curr_recording,:) = ap.groupfun(@sum,+~move_poststim,move_session_bins)./diff(session_bin_edges)';
        p_cue_given_move(curr_recording,:) = ap.groupfun(@mean,+move_poststim,move_session_bins);

        % Get association stat
        [rxn_stat_p(curr_recording), ...
            rxn_stat(curr_recording),rxn_null_stat(curr_recording)] = ...
            AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    bhv(curr_animal_idx).move_rate = move_rate;
    bhv(curr_animal_idx).move_rate_stim = move_rate_stim;
    bhv(curr_animal_idx).move_rate_nonstim = move_rate_nonstim;
    bhv(curr_animal_idx).p_cue_given_move = p_cue_given_move;
    bhv(curr_animal_idx).rxn_stat_p = rxn_stat_p;

    % Clear vars except pre-load for next loop
    clearvars('-except',preload_vars{:});
    AP_print_progress_fraction(curr_recording,length(recordings));

end

ld = cellfun(@(x) ((1:size(x,1)) - find(x<0.05,1))',{bhv.rxn_stat_p}','uni',false);
use_animals = ~cellfun(@isempty,ld);

move_rate_cat = cell2mat({bhv(use_animals).move_rate}');
move_rate_stim_cat = cell2mat({bhv(use_animals).move_rate_stim}');
move_rate_nonstim_cat = cell2mat({bhv(use_animals).move_rate_nonstim}');
p_c2m_cat = cell2mat({bhv(use_animals).p_cue_given_move}');

ld_cat = cell2mat(ld(use_animals));
ld_split = ld_cat + linspace(0,(n_bins-1)/n_bins,n_bins);

[move_rate_avg,move_rate_grp] = ap.groupfun(@mean,move_rate_cat(:),ld_split(:));
move_rate_sem = ap.groupfun(@AP_sem,move_rate_cat(:),ld_split(:));

[move_rate_stim_avg,move_rate_stim_grp] = ap.groupfun(@mean,move_rate_stim_cat(:),ld_split(:));
move_rate_stim_sem = ap.groupfun(@AP_sem,move_rate_stim_cat(:),ld_split(:));

[move_rate_nonstim_avg,move_rate_nonstim_grp] = ap.groupfun(@mean,move_rate_nonstim_cat(:),ld_split(:));
move_rate_nonstim_sem = ap.groupfun(@AP_sem,move_rate_nonstim_cat(:),ld_split(:));

[p_c2m_avg,p_c2m_avg_grp] = ap.groupfun(@mean,p_c2m_cat(:),ld_split(:));
p_c2m_sem = ap.groupfun(@AP_sem,p_c2m_cat(:),ld_split(:));

figure('Name','Fig S1 p stim move'); tiledlayout(2,1);

ax1 = nexttile; hold on;
h1 = errorbar( ...
    padarray(reshape(move_rate_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_sem,n_bins,[]),[1,0],NaN,'post'),'color',[0.5,0.5,0.5],'linewidth',2);
h2 = errorbar( ...
    padarray(reshape(move_rate_stim_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_stim_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_stim_sem,n_bins,[]),[1,0],NaN,'post'),'color',[0,0.5,0],'linewidth',2);
h3 = errorbar( ...
    padarray(reshape(move_rate_nonstim_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_nonstim_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_nonstim_sem,n_bins,[]),[1,0],NaN,'post'),'color',[0.5,0,0],'linewidth',2);
xlabel('Learned day')
ylabel(' Moves/s');
xline(0,'r');
legend([h1(1),h2(1),h3(1)],{'All','Non-stim','Stim'});

ax2 = nexttile;
errorbar( ...
    padarray(reshape(p_c2m_avg_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(p_c2m_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(p_c2m_sem,n_bins,[]),[1,0],NaN,'post'),'k','linewidth',2);
xlabel('Learned day')
ylabel('P(stim|move)');
xline(0,'r');

linkaxes([ax1,ax2],'x');
ap.prettyfig;


%% Non-stim move activity
% (in progress)

%%% Load non-activity data
load_dataset = 'task';
Marica_2025.figures.load_data;
%%%

% Load nonstim move activity
U_master = plab.wf.load_master_U;
load(fullfile(data_path,'nonstim_move'));

% (quick test: pre/post wf nonstim move)
use_nostim_move_recordings = ~cellfun(@isempty,nonstim_move.V_move_nostim_align);

n_vs = 2000;
nostim_move_prepost = plab.wf.svd2px(U_master(:,:,1:n_vs), ...
    ap.groupfun(@nanmean,cat(3,nonstim_move.V_move_nostim_align{:}),[],[], ...
    bhv.days_from_learning(use_nostim_move_recordings)>=0));

stim_move_prepost = plab.wf.svd2px(U_master(:,:,1:n_vs), ...
    ap.groupfun(@nanmean,cell2mat(permute(cellfun(@(x) ...
    permute(nanmean(x(:,:,:,2),1),[3,2,1,4]),wf.V_event_align(use_nostim_move_recordings),'uni',false),[2,3,1])),[],[], ...
    bhv.days_from_learning(use_nostim_move_recordings)>=0));

ap.imscroll(cat(4,stim_move_prepost,nostim_move_prepost));
axis image
colormap(gca,ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);

% Get widefield ROIs for no stim moves
wf_nostim_move_striatum_roi = cell2mat(cellfun(@(x) ...
    permute(ap.wf_roi(U_master,x,[],[],striatum_wf_roi), ...
    [3,2,1]),nonstim_move.V_move_nostim_align(use_nostim_move_recordings), ...
    'uni',false));

% Get nonstim move ephys
% (sum into domain multiunit)
striatum_nostim_move_mua_sum = cellfun(@(mua,domain_idx) ...
    permute(ap.groupfun(@sum,mua,domain_idx,[]),[3,2,1]), ...
    nonstim_move.binned_msn_spikes_move_nostim_align,domain_idx_rec,'uni',false);

striatum_nostim_move_mua = cellfun(mua_norm_smooth_reshape_fcn, ...
    striatum_nostim_move_mua_sum,cellfun(@(x) nanmean(x,1),mua_sub_baseline,'uni',false),mua_div_baseline,'uni',false, ...
    'ErrorHandler',@(varargin) []); %#ok<FUNFUN>

% Plot average stim/non-stim move for domain across learning day
rxn_bin = [-Inf,Inf];
striatum_mua_rxngrp = discretize(striatum_mua_grp.rxn,rxn_bin);

[striatum_mua_moveavg,striatum_mua_moveavg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},striatum_mua(:,:,2), ...
    striatum_mua_grp.animal, ...
    [striatum_mua_grp.domain_idx, striatum_mua_grp.ld,striatum_mua_rxngrp]);

figure;h = tiledlayout(n_domains,5);
for curr_domain = 1:n_domains
    for curr_ld = -2:2

        curr_domain_nonstim_move_act = ...
            nanmean(cell2mat(cellfun(@(mua,domain) mua(domain==curr_domain,:), ...
            striatum_nostim_move_mua(use_nostim_move_recordings & bhv.days_from_learning==curr_ld), ...
            striatum_domain_idx(use_nostim_move_recordings & bhv.days_from_learning==curr_ld),'uni',false)),1);

        curr_stim_move_act_idx = ismember(striatum_mua_moveavg_grp(:,1:2), ...
            [curr_domain,curr_ld],'rows');

        nexttile; hold on
        plot(curr_domain_nonstim_move_act,'k');
        plot(striatum_mua_moveavg(curr_stim_move_act_idx,:)');
        plot((striatum_mua_moveavg(curr_stim_move_act_idx,:)-curr_domain_nonstim_move_act)');

    end
end
linkaxes(h.Children,'xy');

% Plot striatum overlaid
ld_color = ap.colormap('BKR',5);
figure;h = tiledlayout(n_domains,3);
for curr_domain = 1:n_domains

    h1 = nexttile; hold on
    set(h1,'ColorOrder',ld_color);

    h2 = nexttile; hold on
    set(h2,'ColorOrder',ld_color);

    h3 = nexttile; hold on
    set(h3,'ColorOrder',ld_color);

    for curr_ld = -2:2
        curr_domain_nonstim_move_act = ...
            nanmean(cell2mat(cellfun(@(mua,domain) mua(domain==curr_domain,:), ...
            striatum_nostim_move_mua(use_nostim_move_recordings & bhv.days_from_learning==curr_ld), ...
            striatum_domain_idx(use_nostim_move_recordings & bhv.days_from_learning==curr_ld),'uni',false)),1);

        curr_stim_move_act_idx = ismember(striatum_mua_moveavg_grp(:,1:2), ...
            [curr_domain,curr_ld],'rows');

        plot(h1,curr_domain_nonstim_move_act,'linewidth',2);
        plot(h2,striatum_mua_moveavg(curr_stim_move_act_idx,:)','linewidth',2);
        plot(h3,striatum_mua_moveavg(curr_stim_move_act_idx,:)' - ...
            curr_domain_nonstim_move_act,'linewidth',2);
    end
end
linkaxes(h.Children,'xy');

% Plot widefield overlaid
[wf_moveavg,wf_moveavg_grp] = ...
    ap.nestgroupfun({@nanmean,@nanmean},wf_striatum_roi(:,:,:,2), ...
    wf_grp.animal,wf_grp.ld);

ld_color = ap.colormap('BKR',5);
figure;h = tiledlayout(n_domains,2);
for curr_domain = 1:n_domains

    h1 = nexttile; hold on
    set(h1,'ColorOrder',ld_color);

    h2 = nexttile; hold on
    set(h2,'ColorOrder',ld_color);

    for curr_ld = -2:2
        curr_roi_nostim_move_act = nanmean(wf_nostim_move_striatum_roi( ...
            bhv.days_from_learning(use_nostim_move_recordings)==curr_ld,:,curr_domain),1);

        curr_roi_stim_move_act = wf_moveavg(wf_moveavg_grp==curr_ld,:,curr_domain);

        plot(h1,curr_roi_nostim_move_act,'linewidth',2);
        plot(h2,curr_roi_stim_move_act,'linewidth',2);
    end
end
linkaxes(h.Children,'xy');


% Plot wheel velocity for stim/nostim
move_stim_wheel_avg = cell2mat(cellfun(@nanmean,nonstim_move.move_stim_wheel(use_nostim_move_recordings),'uni',false));
move_nostim_wheel_avg = cell2mat(cellfun(@nanmean,nonstim_move.move_nostim_wheel(use_nostim_move_recordings),'uni',false));

figure;h = tiledlayout(1,2);
h1 = nexttile; hold on
set(h1,'ColorOrder',ld_color);
h2 = nexttile; hold on
set(h2,'ColorOrder',ld_color);
for curr_ld = -2:2
    plot_idx = bhv.days_from_learning(use_nostim_move_recordings) == curr_ld;
    plot(h1,nanmean(move_stim_wheel_avg(plot_idx,:),1),'linewidth',2);
    plot(h2,nanmean(move_nostim_wheel_avg(plot_idx,:),1),'linewidth',2);
end
linkaxes(h.Children,'xy');



%% Fraction of quiescent passive trials

%%% Load non-activity data
load_dataset = 'noact';
Marica_2025.figures.load_data;
%%%

% Get which trials in passive were quiescent
animals = unique(bhv.animal,'stable');

passive_quiescent_trials = struct( ...
    'stim_x',cell(length(animals),1), ...
    'quiescent',cell(length(animals),1));

for animal_idx=1:length(animals)
   
    animal = animals{animal_idx};

    % Find passive recording days that also have task
    workflow_passive = {'lcr_passive'};
    recordings_passive = plab.find_recordings(animal, [], workflow_passive);
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);
    training_days = ismember({recordings_passive.day}, {recordings_task.day});
    train_rec_passive = recordings_passive(training_days);
    bhv_days = {train_rec_passive.day};
    ephys_days =  bhv_days([train_rec_passive.ephys]);    

    for use_rec = 1:length(ephys_days)

        rec_day = train_rec_passive(use_rec).day;
        rec_time = train_rec_passive(use_rec).recording{end};
        load_parts.behavior = true;
        ap.load_recording
   
        % Trial stim values
        trial_stim_values = vertcat(trial_events.values.TrialStimX);
        trial_stim_values = trial_stim_values(1:length(stimOn_times));

        % Quiescent trials
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');

        passive_quiescent_trials(animal_idx).stim_x{use_rec,1} = trial_stim_values;
        passive_quiescent_trials(animal_idx).quiescent{use_rec,1} = quiescent_trials;

    end
    ap.print_progress_fraction(animal_idx,length(animals));
end


% Plot quiescent trials pre/post learning
passive_quiescent_stim = cell2mat(cellfun(@(stim,quiescent) ...
    ap.groupfun(@mean,+quiescent,stim)', ...
    cat(1,passive_quiescent_trials.stim_x), ...
    cat(1,passive_quiescent_trials.quiescent),'uni',false));

[passive_quiescent_stim_ld,passive_quiescent_stim_ld_grp] = ...
    ap.groupfun(@mean,passive_quiescent_stim,bhv.days_from_learning >= 0);
passive_quiescent_stim_ld_sem = ...
    ap.groupfun(@ap.sem,passive_quiescent_stim,bhv.days_from_learning >= 0);

figure;
x_labels = ["Pre-learn","Post-learn"];
errorbar(reordercats(categorical(x_labels),x_labels), ...
    passive_quiescent_stim_ld,passive_quiescent_stim_ld_sem,'linewidth',2)
axis padded;
ylabel('Frac. quiescent trials')
ap.prettyfig;


% Plot quiescent trials binned by trial number
n_trialsplit = 9;
passive_trialsplit_idx = cell2mat(cellfun(@(stim) ...
    (floor((0:length(stim)-1)/n_trialsplit)+1)', ...
    cat(1,passive_quiescent_trials.stim_x),'uni',false));

passive_ld_idx = cell2mat(cellfun(@(rec,stim) repelem(rec,length(stim))', ...
    num2cell(bhv.days_from_learning), ...
    cat(1,passive_quiescent_trials.stim_x),'uni',false));

[passive_quiescent_stim_trialsplit,passive_quiescent_stim_trialsplit_grp] = ...
    ap.groupfun(@mean,+cell2mat(cat(1,passive_quiescent_trials.quiescent)), ...
    [passive_ld_idx,passive_trialsplit_idx,cell2mat(cat(1,passive_quiescent_trials.stim_x))]);

passive_quiescent_stim_trialsplit_grid = ...
    accumarray([grp2idx(passive_quiescent_stim_trialsplit_grp(:,1)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,2)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,3))], ...
    passive_quiescent_stim_trialsplit,[],[],nan);

passive_quiescent_stim_daysplit_x = reshape((unique(passive_quiescent_stim_trialsplit_grp(:,1)) + ...
    linspace(0,1,max(passive_quiescent_stim_trialsplit_grp(:,2)+1)))',[],1);

figure;
plot(passive_quiescent_stim_daysplit_x, ...
    reshape(permute(padarray(passive_quiescent_stim_trialsplit_grid,[0,1],nan,'post'),[2,1,3]),[],3))
xlabel('Days from learning');
ylabel('Frac. quiescent trials');










