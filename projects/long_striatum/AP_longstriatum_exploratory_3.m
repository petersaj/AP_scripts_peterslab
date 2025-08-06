%% Notes

% After starting figure code: tests for figures
% (using preprocessing in AP_longstriatum_load_data)

%% N days from learning

[n,ld] = ap.groupfun(@sum,ones(size(bhv.days_from_learning)),bhv.days_from_learning);
figure;plot(ld,n);


%% Striatal passive TAN/celltype responses by stim

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

%%%%%% SPIKE PREPROCESSING TESTING:

%%%%% don't normalize spikes, just smooth
striatum_sua = smoothdata(cell2mat(cellfun(@(x,y) x(y,:,:), ...
    ephys.unit_event_psths,striatum_units,'uni',false)),2,'gaussian',[100,0]);


%%%%% normalize to baseline and smooth
striatum_sua = cell2mat(cellfun(@(data,baseline,striatum_units) ...
    spikes_norm_smooth_reshape_fcn(data(striatum_units,:,:),baseline(striatum_units)), ...
    ephys.unit_event_psths,sua_baseline,striatum_units,'uni',false));


%%%%%%%%%%%%%%%%%%

plot_day_bins = [-inf,0,inf];
plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);



stim_t = psth_t > 0.05 & psth_t < 0.15;
n_stim = size(striatum_sua,3);

% Plot heatmap
for curr_domain = 1:n_domains
    figure;
    colormap(ap.colormap('WK',[],2));
    h = tiledlayout(n_stim,max(plot_day_grp),'TileSpacing','compact');
    for curr_stim = 1:n_stim
        for curr_day_grp = 1:length(plot_day_bins)-1
            nexttile;

            curr_units = find(striatum_sua_grp.domain_idx == curr_domain & ...
                plot_day_grp == curr_day_grp & ...
                striatum_sua_grp.tan);

            % (sort max across stim)
            [~,sort_idx] = sort(max(mean(striatum_sua(curr_units,stim_t,:),2),[],3),'descend');

            imagesc(psth_t,[],striatum_sua(curr_units(sort_idx),:,curr_stim));
            clim([0,20])
            xlim([-0.2,0.8])
            title(plot_day_bins(curr_day_grp));

            % (unused: get rec/IDs of cells)
            curr_sorted_unit_coordinate = ...
                [ephys.animal(striatum_sua_grp.rec(curr_units(sort_idx))), ...
                ephys.rec_day(striatum_sua_grp.rec(curr_units(sort_idx))), ...
                num2cell(striatum_sua_grp.unit_id(curr_units(sort_idx)))];

        end
    end
    title(h,sprintf('Striatal cluster %d',curr_domain));
    linkaxes(h.Children,'xy');
    ap.prettyfig;
end

% Plot scatter and diagonal histograms
stim_t = psth_t > 0.05 & psth_t < 0.15;
compare_stim = [2,3];
scatter_fig = figure;
h_scatter = tiledlayout(scatter_fig,n_domains,max(plot_day_grp),'TileSpacing','compact');

histogram_fig = figure;
h_histogram = tiledlayout(histogram_fig,n_domains,max(plot_day_grp),'TileSpacing','compact');
for curr_domain = 1:n_domains
    for curr_day_grp = 1:length(plot_day_bins)-1
        % (plot scatter)
        nexttile(h_scatter); axis equal; hold on;

        curr_units = find(striatum_sua_grp.domain_idx == curr_domain & ...
            plot_day_grp == curr_day_grp & ...
            striatum_sua_grp.tan);

        curr_unit_mean = squeeze(mean(striatum_sua(curr_units,stim_t,:),2));

        plot(curr_unit_mean(:,compare_stim(1)),curr_unit_mean(:,compare_stim(2)),'.k');
        xlabel(sprintf('Stim %d',compare_stim(1)));
        ylabel(sprintf('Stim %d',compare_stim(2)));

        % (~stats~)
        diff_meas = mean(diff(curr_unit_mean(:,compare_stim),[],2));
        n_shuff = 10000;
        diff_null = nan(n_shuff,1);
        for curr_shuff = 1:n_shuff
            diff_null(curr_shuff) = mean(diff(ap.shake(curr_unit_mean(:,compare_stim),2),[],2));
        end
        stat_rank = tiedrank([diff_meas;diff_null]);
        stat_p = 1-stat_rank(1)/(n_shuff+1);

        title(sprintf('%d, p = %.2g',plot_day_bins(curr_day_grp),stat_p));

        % (plot diagonal histogram)
        curr_unit_mean_stimdiff = diff(curr_unit_mean(:,compare_stim),[],2);
        hist_edges = -15:1:15;
        nexttile(h_histogram);
        histogram(-curr_unit_mean_stimdiff,hist_edges, ...
            'normalization','probability','FaceColor','k','EdgeColor','none')
        xline(0,'color',[0.5,0.5,0.5]);

    end
end
% (draw zero and guide lines)
scatter_max = 25;
hist_max = 10;
arrayfun(@(x) line(x,xlim,xlim,'color',[0.5,0.5,0.5]),h_scatter.Children);
arrayfun(@(x) line(x,xlim,xlim-hist_max,'color','r'),h_scatter.Children);
arrayfun(@(x) line(x,xlim,xlim+hist_max,'color','r'),h_scatter.Children);
linkaxes(h_scatter.Children,'xy');
[h_scatter.Children.XLim,h_scatter.Children.YLim] = deal([0,scatter_max]);

linkaxes(h_histogram.Children,'xy');
[h_histogram.Children.XLim] = deal([-1,1].*hist_max);
arrayfun(@(x) xline(x,[-1,1].*hist_max,'r'),h_histogram.Children);

figure(scatter_fig);ap.prettyfig;
figure(histogram_fig);ap.prettyfig;


% Plot average
figure;
h = tiledlayout(n_domains,n_stim*max(plot_day_grp),'TileSpacing','compact');
stim_colormap = ap.colormap('BKR',n_stim);
for curr_domain = 1:n_domains
    for curr_day_grp = 1:length(plot_day_bins)-1

        curr_units = find(striatum_sua_grp.domain_idx == curr_domain & ...
            plot_day_grp == curr_day_grp & ...
            striatum_sua_grp.tan);

        % curr_sua_mean = mean(striatum_sua(curr_units,:,:),1);
        % curr_sua_sem = AP_sem(striatum_sua(curr_units,:,:),1);

        curr_sua_mean = mean(ap.groupfun(@mean,striatum_sua(curr_units,:,:),striatum_sua_grp.animal(curr_units)),1);
        curr_sua_sem = AP_sem(ap.groupfun(@mean,striatum_sua(curr_units,:,:),striatum_sua_grp.animal(curr_units)),1);

        for curr_stim = 1:3
            nexttile;
            ap.errorfill(psth_t,curr_sua_mean(:,:,curr_stim),curr_sua_sem(:,:,curr_stim),stim_colormap(curr_stim,:));
        end

    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;
xlim([-0.2,0.8]);l


%%%%% TESTING: mean of stim diff across cell types
celltype_order = {'msn','fsi','tan'};
celltype_id = sum([striatum_sua_grp.(celltype_order{1}), ...
    striatum_sua_grp.(celltype_order{2}), ...
    striatum_sua_grp.(celltype_order{3})].*[1,2,3],2);

figure;
h = tiledlayout(n_domains,max(plot_day_grp),'TileSpacing','compact');
for curr_domain = 1:n_domains
    for curr_day_grp = 1:max(plot_day_grp)

        curr_units = find(striatum_sua_grp.domain_idx == curr_domain & ...
            plot_day_grp == curr_day_grp & celltype_id ~= 0);

        curr_unit_mean = squeeze(mean(striatum_sua(curr_units,stim_t,:),2));

        rc_diff = diff(curr_unit_mean(:,[2,3]),[],2)./ ...
            sum(curr_unit_mean(:,[2,3]),2);

        rl_diff = diff(curr_unit_mean(:,[1,3]),[],2)./ ...
            sum(curr_unit_mean(:,[1,3]),2);

        nexttile; hold on;
        violinplot(categorical(celltype_order(celltype_id(curr_units))),rc_diff,'DensityDirection','negative');
        violinplot(categorical(celltype_order(celltype_id(curr_units))),rl_diff,'DensityDirection','positive');

        rc_diff_med = ap.groupfun(@nanmedian,rc_diff,celltype_id(curr_units));
        rl_diff_med = ap.groupfun(@nanmedian,rl_diff,celltype_id(curr_units));

        plot(categorical(celltype_order),rc_diff_med,'_b','MarkerSize',10);
        plot(categorical(celltype_order),rl_diff_med,'_r','MarkerSize',10);

        yline(0)
        ylabel('Diff/Sum')

    end
end
legend({'R-C','R-L'});


%%%%% TESTING: alt version of above: exclude L response
celltype_order = {'msn','fsi','tan'};
celltype_id = sum([striatum_sua_grp.(celltype_order{1}), ...
    striatum_sua_grp.(celltype_order{2}), ...
    striatum_sua_grp.(celltype_order{3})].*[1,2,3],2);

figure;
h = tiledlayout(n_domains,1,'TileSpacing','compact');
for curr_domain = 1:n_domains

    nexttile; hold on;

    curr_units = find(striatum_sua_grp.domain_idx == curr_domain & ...
        plot_day_grp == 1 & celltype_id ~= 0);

    curr_unit_mean = squeeze(mean(striatum_sua(curr_units,stim_t,:),2));

    rc_diff = diff(curr_unit_mean(:,[2,3]),[],2)./ ...
        sum(curr_unit_mean(:,[2,3]),2);

    violinplot(categorical(celltype_order(celltype_id(curr_units))),rc_diff,'DensityDirection','negative');
    rc_diff_med = ap.groupfun(@nanmedian,rc_diff,celltype_id(curr_units));
    plot(categorical(celltype_order),rc_diff_med,'_b','MarkerSize',10);

    curr_units = find(striatum_sua_grp.domain_idx == curr_domain & ...
        plot_day_grp == 2 & celltype_id ~= 0);

    curr_unit_mean = squeeze(mean(striatum_sua(curr_units,stim_t,:),2));

    rc_diff = diff(curr_unit_mean(:,[2,3]),[],2)./ ...
        sum(curr_unit_mean(:,[2,3]),2);

    violinplot(categorical(celltype_order(celltype_id(curr_units))),rc_diff,'DensityDirection','positive');
    rc_diff_med = ap.groupfun(@nanmedian,rc_diff,celltype_id(curr_units));
    plot(categorical(celltype_order),rc_diff_med,'_r','MarkerSize',10);

    yline(0)
    ylabel('R-C')

end
legend({'Day group 1','Day group 2'});
ap.prettyfig;

%%%%% TESTING: mean of stim diff across cell types, but as sem errorbars
celltype_order = {'msn','fsi','tan'};
celltype_id = sum([striatum_sua_grp.(celltype_order{1}), ...
    striatum_sua_grp.(celltype_order{2}), ...
    striatum_sua_grp.(celltype_order{3})].*[1,2,3],2);

figure;
h = tiledlayout(n_domains,max(plot_day_grp),'TileSpacing','compact');
for curr_domain = 1:n_domains
    for curr_day_grp = 1:max(plot_day_grp)

        curr_units = find(striatum_sua_grp.domain_idx == curr_domain & ...
            plot_day_grp == curr_day_grp & celltype_id ~= 0);

        curr_unit_mean = squeeze(mean(striatum_sua(curr_units,stim_t,:),2));

        rc_diff = diff(curr_unit_mean(:,[2,3]),[],2)./ ...
            sum(curr_unit_mean(:,[2,3]),2);

        rl_diff = diff(curr_unit_mean(:,[1,3]),[],2)./ ...
            sum(curr_unit_mean(:,[1,3]),2);

        nexttile; hold on;

        rc_diff_med = ap.groupfun(@nanmedian,rc_diff,celltype_id(curr_units));
        rl_diff_med = ap.groupfun(@nanmedian,rl_diff,celltype_id(curr_units));

        rc_diff_error = ap.groupfun(@AP_sem,rc_diff,celltype_id(curr_units));
        rl_diff_error = ap.groupfun(@AP_sem,rl_diff,celltype_id(curr_units));

        errorbar(categorical(celltype_order),rc_diff_med,rc_diff_error,'.b');
        errorbar(categorical(celltype_order),rl_diff_med,rl_diff_error,'.r');

        yline(0)
        ylabel('Diff/Sum')

    end
end
linkaxes(h.Children,'xy');
legend({'R-C','R-L'});

%%%%% TESTING: mean of stim diff across cell types, but as ci errorbars
stim_t = psth_t > 0.05 & psth_t < 0.15;

celltype_order = {'msn','fsi','tan'};
celltype_id = sum([striatum_sua_grp.(celltype_order{1}), ...
    striatum_sua_grp.(celltype_order{2}), ...
    striatum_sua_grp.(celltype_order{3})].*[1,2,3],2);

figure;
h = tiledlayout(n_domains,max(plot_day_grp),'TileSpacing','compact');
for curr_domain = 1:n_domains
    for curr_day_grp = 1:max(plot_day_grp)

        curr_units = find(striatum_sua_grp.domain_idx == curr_domain & ...
            plot_day_grp == curr_day_grp & celltype_id ~= 0);

        curr_unit_mean = squeeze(mean(striatum_sua(curr_units,stim_t,:),2));

        % just diff (baseline-norm data)
        rc_diff = diff(curr_unit_mean(:,[2,3]),[],2);
        rl_diff = diff(curr_unit_mean(:,[1,3]),[],2);

        n_shuff = 1000;
        rc_shuff_med = nan(n_shuff,3);
        rl_shuff_med = nan(n_shuff,3);
        for curr_shuff = 1:n_shuff

            data_shuff = ap.shake(curr_unit_mean(:,[2,3]),2);
            shuff_diff = diff(data_shuff,[],2);
            rc_shuff_med(curr_shuff,:) = ...
                ap.groupfun(@nanmean,shuff_diff,celltype_id(curr_units));

            data_shuff = ap.shake(curr_unit_mean(:,[1,3]),2);
            shuff_diff = diff(data_shuff,[],2);
            rl_shuff_med(curr_shuff,:) = ...
                ap.groupfun(@nanmean,shuff_diff,celltype_id(curr_units));
        end
        rc_shuff_ci = prctile(rc_shuff_med,[2.5,97.5],1);
        rl_shuff_ci = prctile(rl_shuff_med,[2.5,97.5],1);

        nexttile; hold on;

        rc_diff_med = ap.groupfun(@nanmean,rc_diff,celltype_id(curr_units));
        rl_diff_med = ap.groupfun(@nanmean,rl_diff,celltype_id(curr_units));

        % plot x-categories in given order
        celltype_x = reordercats(categorical(celltype_order),celltype_order);

        plot(celltype_x,rc_diff_med,'.b','MarkerSize',20);
        plot(celltype_x,rl_diff_med,'.r','MarkerSize',20);

        plot(celltype_x,rc_shuff_ci,'_b','MarkerSize',20);
        plot(celltype_x,rl_shuff_ci,'_r','MarkerSize',20);

        yline(0)
        ylabel('Diff')

    end
end
linkaxes(h.Children,'xy');
legend({'R-C','R-L'});
ap.prettyfig;

%%%%% TESTING: mean of stim diff across cell types, but as ci errorbars
% (average across animals)
stim_t = psth_t > 0.05 & psth_t < 0.15;

celltype_order = {'msn','fsi','tan'};
celltype_id = sum([striatum_sua_grp.(celltype_order{1}), ...
    striatum_sua_grp.(celltype_order{2}), ...
    striatum_sua_grp.(celltype_order{3})].*[1,2,3],2);

figure;
h = tiledlayout(n_domains,max(plot_day_grp),'TileSpacing','compact');
for curr_domain = 1:n_domains
    for curr_day_grp = 1:max(plot_day_grp)

        curr_units = find(striatum_sua_grp.domain_idx == curr_domain & ...
            plot_day_grp == curr_day_grp & celltype_id ~= 0);

        curr_unit_mean = squeeze(mean(striatum_sua(curr_units,stim_t,:),2));

        % just diff (baseline-norm data)
        rc_diff = diff(curr_unit_mean(:,[2,3]),[],2);
        rl_diff = diff(curr_unit_mean(:,[1,3]),[],2);

        n_shuff = 1000;
        rc_shuff_med = nan(n_shuff,3);
        rl_shuff_med = nan(n_shuff,3);
        for curr_shuff = 1:n_shuff

            data_shuff = ap.shake(curr_unit_mean(:,[2,3]),2);
            shuff_diff = diff(data_shuff,[],2);
            rc_shuff_med(curr_shuff,:) = ...
                ap.nestgroupfun({@nanmean,@nanmean},shuff_diff,striatum_sua_grp.animal(curr_units),celltype_id(curr_units));

            data_shuff = ap.shake(curr_unit_mean(:,[1,3]),2);
            shuff_diff = diff(data_shuff,[],2);
            rl_shuff_med(curr_shuff,:) = ...
                ap.nestgroupfun({@nanmean,@nanmean},shuff_diff,striatum_sua_grp.animal(curr_units),celltype_id(curr_units));
        end
        rc_shuff_ci = prctile(rc_shuff_med,[5,95],1);
        rl_shuff_ci = prctile(rl_shuff_med,[5,95],1);

        nexttile; hold on;

        rc_diff_med = ap.nestgroupfun({@nanmean,@nanmean},rc_diff,striatum_sua_grp.animal(curr_units),celltype_id(curr_units));
        rl_diff_med = ap.nestgroupfun({@nanmean,@nanmean},rl_diff,striatum_sua_grp.animal(curr_units),celltype_id(curr_units));

        % plot x-categories in given order
        celltype_x = reordercats(categorical(celltype_order),celltype_order);

        plot(celltype_x,rc_diff_med,'.b','MarkerSize',20);
        plot(celltype_x,rl_diff_med,'.r','MarkerSize',20);

        plot(celltype_x,rc_shuff_ci,'_b','MarkerSize',20);
        plot(celltype_x,rl_shuff_ci,'_r','MarkerSize',20);

        yline(0)
        ylabel('Diff')

    end
end
linkaxes(h.Children,'xy');
legend({'R-C','R-L'});
ap.prettyfig;

%% (above, v2)

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

% Set striatal domains to plot (combines if multiple)
plot_domains = 1:2;

% Set days to group
plot_day_bins = [-inf,-2,0,Inf];
plot_day_grp = discretize(max(-inf,striatum_sua_grp.ld),plot_day_bins);

% Get max activity in type-relevant window
celltypes = ["msn","fsi","tan"];
striatum_sua_tavg = nan(size(striatum_sua,[1,3]));
for curr_celltype = celltypes

    curr_units = find(striatum_sua_grp.(curr_celltype));

    % Get window from mean activity across all cells
    curr_sua_mean = mean(ap.groupfun(@mean, ...
        striatum_sua(curr_units,:,3),striatum_sua_grp.animal(curr_units)),1);

    [curr_sua_mean_max,curr_sua_mean_max_idx] = max(curr_sua_mean);
    act_thresh = curr_sua_mean_max*0.2;

    act_window = ...
        psth_t([find(curr_sua_mean < act_thresh & ...
        1:length(curr_sua_mean) < curr_sua_mean_max_idx,1,'last'), ...
        find(curr_sua_mean < act_thresh & ...
        1:length(curr_sua_mean) > curr_sua_mean_max_idx,1,'first')]);

    stim_t = psth_t > act_window(1) & psth_t < act_window(2);

    % Get max activity within window
    striatum_sua_tavg(curr_units,:) = mean(striatum_sua(curr_units,stim_t,:),2);
end

n_stim = size(striatum_sua,3);

% Plot heatmap
for curr_celltype = celltypes
    figure;
    colormap(ap.colormap('BWR',[],2));
    h = tiledlayout(max(plot_day_grp),n_stim,'TileSpacing','tight');
    title(h,curr_celltype);
    for curr_day_grp = 1:length(plot_day_bins)-1

        % Get units and sorting
        curr_units = find( ...
            ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            plot_day_grp == curr_day_grp & ...
            striatum_sua_grp.(curr_celltype));

        % (sort max stim, then max within stim)
        sort_idx_cell = cell(3,1);
        [~,max_stim] = max(striatum_sua_tavg,[],2);
        for curr_stim_sort = 1:3
            curr_stim_units = find(max_stim(curr_units)==curr_stim_sort);
            [~,curr_sort_idx] = sort(max(striatum_sua_tavg(curr_units(curr_stim_units),curr_stim_sort),[],2),'descend');
            sort_idx_cell{curr_stim_sort} = curr_stim_units(curr_sort_idx);
        end
        sort_idx = cell2mat(sort_idx_cell);
        plot_units = curr_units(sort_idx);

        % (get rec/IDs of cells for single-unit PSTH plotting)
        curr_sorted_unit_coordinate = ...
            [ephys.animal(striatum_sua_grp.rec(plot_units)), ...
            ephys.rec_day(striatum_sua_grp.rec(plot_units)), ...
            num2cell(striatum_sua_grp.unit_id(plot_units))];

        % (smooth MSN/FSI because large number of units;
        max_n_cells = max(ap.groupfun(@sum, ...
            +(ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            striatum_sua_grp.(curr_celltype)),plot_day_grp));
        smooth_n = 5*max_n_cells/500;

        for curr_stim = 1:n_stim
            nexttile;
            imagesc(psth_t,[],movmean(striatum_sua(plot_units,:,curr_stim),[smooth_n,0],1));
            clim([-1,1])
            xlim([-0.2,0.8])
            yline(cumsum(cellfun(@length,sort_idx_cell)),'k');
            axis off;
        end

    end
    linkaxes(h.Children,'xy');
    ap.prettyfig;
end

% % Plot example PSTHs
% psth_unit_coordinates = { ...
%     'AM026','2024-07-25',354; ... % pre-learn
%     'AM022','2024-04-06',246};    % post-learn
%
% for curr_unit = 1:size(psth_unit_coordinates,1)
%     animal = psth_unit_coordinates{curr_unit,1};
%     rec_day = psth_unit_coordinates{curr_unit,2};
%     plot_units = psth_unit_coordinates{curr_unit,3};
%     AP_longstriatum_psth_fig;
% end

% % %%% TESTING: Plot example PSTHs
% % curr_unit = 1141;
% curr_unit = find(max_stim(plot_units)==3,1)+3;
% animal = curr_sorted_unit_coordinate{curr_unit,1};
% rec_day = curr_sorted_unit_coordinate{curr_unit,2};
% plot_units = curr_sorted_unit_coordinate{curr_unit,3};
% AP_longstriatum_psth_fig;

% Plot average
for curr_celltype = celltypes
    figure;
    h = tiledlayout(1,max(plot_day_grp),'TileSpacing','compact');
    title(h,curr_celltype);
    stim_colormap = ap.colormap('BKR',n_stim);
    for curr_day_grp = 1:length(plot_day_bins)-1
        curr_units = find(...
            ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            plot_day_grp == curr_day_grp & ...
            striatum_sua_grp.(curr_celltype));

        curr_sua_mean = squeeze(mean(ap.groupfun(@mean, ...
            striatum_sua(curr_units,:,:),striatum_sua_grp.animal(curr_units)),1));
        curr_sua_sem = squeeze(AP_sem(ap.groupfun(@mean, ...
            striatum_sua(curr_units,:,:),striatum_sua_grp.animal(curr_units)),1));

        nexttile;
        ap.errorfill(psth_t,curr_sua_mean,curr_sua_sem,stim_colormap);

    end
    linkaxes(h.Children,'xy');
    ap.prettyfig;
    xlim(h.Children(1),[-0.2,0.8]);
end

% Plot average by stim response
[~,max_stim] = max(striatum_sua_tavg,[],2);
for curr_celltype = celltypes
    figure;
    h = tiledlayout(3,max(plot_day_grp),'TileSpacing','compact','TileIndexing','ColumnMajor');
    title(h,curr_celltype);
    stim_colormap = ap.colormap('BKR',n_stim);
    for curr_day_grp = 1:length(plot_day_bins)-1
        for curr_max_stim = 1:3
            nexttile; hold on; set(gca,'ColorOrder',stim_colormap);
            curr_stim_units = ...
                ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
                plot_day_grp == curr_day_grp & ...
                striatum_sua_grp.(curr_celltype) & ...
                max_stim == curr_max_stim;
            curr_sua_mean = squeeze(mean(ap.groupfun(@mean, ...
                striatum_sua(curr_stim_units,:,:),striatum_sua_grp.animal(curr_stim_units)),1));
            plot(psth_t,curr_sua_mean);
        end
    end
    linkaxes(h.Children,'xy');
    ap.prettyfig;
    xlim(h.Children(1),[-0.2,0.8]);
end

% Plot max by stim response
[~,max_stim] = max(striatum_sua_tavg,[],2);
for curr_celltype = celltypes
    figure;
    h = tiledlayout(n_stim*max(plot_day_grp),1,'TileSpacing','none');
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
    linkaxes(h.Children,'xy');
    ap.prettyfig;
end

% % Plot max by stim response (group by day)
% [~,max_stim] = max(striatum_sua_tavg,[],2);
% for curr_celltype = celltypes
%     figure;
%     h = tiledlayout(max(plot_day_grp),1,'TileSpacing','tight');
%     title(h,curr_celltype);
%
%     stim_colormap = ap.colormap('BKR',n_stim);
%     for curr_day_grp = 1:3
%
%         curr_units = ...
%             ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
%             striatum_sua_grp.(curr_celltype) & ...
%             plot_day_grp == curr_day_grp;
%
%         curr_act_mean = ap.nestgroupfun({@nanmean,@nanmean},striatum_sua_tavg(curr_units,:), ...
%             striatum_sua_grp.animal(curr_units),max_stim(curr_units));
%
%         curr_act_sem = ap.nestgroupfun({@nanmean,@AP_sem},striatum_sua_tavg(curr_units,:), ...
%             striatum_sua_grp.animal(curr_units),max_stim(curr_units));
%
%         nexttile; hold on;colormap(stim_colormap);
%         b = bar(curr_act_mean,'FaceColor','flat');
%         stim_cdata = num2cell(1:3);
%         [b.CData] = deal(stim_cdata{:});
%         errorbar(vertcat(b.XEndPoints)',curr_act_mean,curr_act_sem, ...
%             'marker','none','linestyle','none','color','k','linewidth',1);
%
%     end
%     xlabel('Preferred stimulus');
%     linkaxes(h.Children,'xy');
%     ap.prettyfig;
% end



% Stim response by celltype
plot_celltypes = ["tan","fsi","msn"];
celltype_id = sum(cell2mat(arrayfun(@(x) ...
    striatum_sua_grp.(plot_celltypes(x)).*x, ...
    1:length(plot_celltypes),'uni',false)),2);

stim_colormap = ap.colormap('BWR',n_stim);

figure;
h = tiledlayout(1,max(plot_day_grp),'TileSpacing','compact');
for curr_day_grp = 1:max(plot_day_grp)

    curr_units = find( ...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        plot_day_grp == curr_day_grp & celltype_id ~= 0);

    curr_act_mean = ap.nestgroupfun({@nanmean,@nanmean},striatum_sua_tavg(curr_units,:), ...
        striatum_sua_grp.animal(curr_units),celltype_id(curr_units));

    curr_act_sem = ap.nestgroupfun({@nanmean,@AP_sem},striatum_sua_tavg(curr_units,:), ...
        striatum_sua_grp.animal(curr_units),celltype_id(curr_units));

    nexttile; hold on; set(gca,'ColorOrder',stim_colormap);
    celltype_x = reordercats(categorical(plot_celltypes),plot_celltypes);
    b = bar(celltype_x,curr_act_mean);
    errorbar(vertcat(b.XEndPoints)',curr_act_mean,curr_act_sem, ...
        'marker','none','linestyle','none','color','k','linewidth',1);

    % ~stats~ %
    rc_diff = ap.nestgroupfun({@nanmean,@nanmean}, ...
        diff(striatum_sua_tavg(curr_units,[2,3]),[],2), ...
        striatum_sua_grp.animal(curr_units),celltype_id(curr_units));

    n_shuff = 1000;
    rc_diff_shuff = nan(length(plot_celltypes),n_shuff);
    for curr_shuff = 1:n_shuff
        data_shuff = ap.shake(striatum_sua_tavg(curr_units,[2,3]),2);
        rc_diff_shuff(:,curr_shuff) = ...
            ap.nestgroupfun({@nanmean,@nanmean}, ...
            diff(data_shuff,[],2), ...
            striatum_sua_grp.animal(curr_units),celltype_id(curr_units));
    end
    stat_rank = tiedrank([rc_diff,rc_diff_shuff]');
    stat_p = 1-stat_rank(1,:)/(n_shuff+1);
    arrayfun(@(x) fprintf('%s daygroup %d RvC: p = %.3g\n', ...
        plot_celltypes(x),curr_day_grp,stat_p(x)),1:length(plot_celltypes))

end
linkaxes(h.Children,'xy');
legend({'L','C','R'});
ap.prettyfig;


% Stim response by celltype - grouped by celltype
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
    compare_day_grps = [1,3];
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






%% Snowflake plot

% Get max response in window
use_t = psth_t > 0 & psth_t < 0.2;
% (t max)
% unit_psth_max = permute(max(striatum_sua(:,use_t,:),[],2),[1,3,2]);
% (t avg)
% unit_psth_max = permute(mean(striatum_sua(:,use_t,:),2),[1,3,2]);
% (t abs sum)
unit_psth_max = permute(sum(abs(striatum_sua(:,use_t,:)),2),[1,3,2]);

plot_day_bins = [-Inf,0,Inf];
unit_plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

% Get responsive units
unit_responsive_p_thresh = 0.95;
unit_responsive = cell2mat(horzcat(ephys.unit_resp_p_value{:})') > unit_responsive_p_thresh;
striatum_units_responsive = unit_responsive(cell2mat(striatum_units),:);

% Scatter (all)
figure;
h = tiledlayout(1,length(plot_day_bins)-1);
plot_domains = 1:3;
for curr_day_grp = 1:length(plot_day_bins)-1
    nexttile; axis equal
    curr_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        unit_plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.tan;
    scatter(unit_psth_max(curr_units,2),unit_psth_max(curr_units,3),10, ...
        any(striatum_units_responsive(curr_units,:),2),'filled');
    line(ylim,ylim);
end
colormap([0.5,0.5,0.5;1,0,0])
linkaxes(h.Children,'xy');

% Scatter (only responsive units)
figure;
h = tiledlayout(1,length(plot_day_bins)-1);
plot_domains = 1:2;
for curr_day_grp = 1:length(plot_day_bins)-1
    nexttile; axis equal
    curr_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        unit_plot_day_grp == curr_day_grp & ...
        any(striatum_units_responsive(:,2:3),2) & ...
        striatum_sua_grp.msn;

    scatter(unit_psth_max(curr_units,2),unit_psth_max(curr_units,3),10, ...
        striatum_sua_grp.animal(curr_units),'filled');

    line(ylim,ylim);
end
colormap(ap.colormap('tube'))
linkaxes(h.Children,'xy');


% Trying out snowflake plot
unit_response_norm = unit_psth_max;
unit_response_relative = abs(unit_response_norm)./max(abs(unit_response_norm),[],2);

ternary_tform = [0,1;sqrt(3/4),-0.5;-sqrt(3/4),-0.5];
unit_response_norm_ternary = (ternary_tform\unit_response_relative')';

unit_max_response = max(abs(unit_response_norm),[],2);

response_limits = reshape([eye(3),eye(3)+circshift(eye(3),[0,1])]',3,[])';
ternary_limits = (ternary_tform\response_limits')';

figure;
h = tiledlayout(1,length(plot_day_bins)-1);
plot_domains = 1;
for curr_day_grp = 1:length(plot_day_bins)-1
    curr_data_idx = ...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        unit_plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.tan;

    % Plot each unit
    size_norm_denom = prctile(unit_max_response(ismember(striatum_sua_grp.domain_idx,plot_domains)),90);
    size_norm = min(unit_max_response(curr_data_idx)./size_norm_denom,1);
    dot_size = size_norm.*20+1;

    response_sign = sign(ap.signed_max(unit_response_norm(curr_data_idx,:),2));
    dot_color = [response_sign == 1,zeros(size(response_sign)),response_sign == -1];

    nexttile; hold on; axis image off;
    patch(ternary_limits(:,1),ternary_limits(:,2),'w','linewidth',1);
    scatter(unit_response_norm_ternary(curr_data_idx,1), ...
        unit_response_norm_ternary(curr_data_idx,2), ...
        dot_size,dot_color,'filled');
end



%% Plot task heatmap with sorting from passive

% Load task SUA
load(fullfile(data_path,'ephys_task'));

psth_t = -0.5:0.001:1;
baseline_t = psth_t < 0;
softnorm = 1;
sua_baseline = cellfun(@(sua) ...
    mean(sua(:,baseline_t,1),[2,3]), ...
    ephys.unit_event_psths,'uni',false,'ErrorHandler',@(varargin) NaN);

striatum_sua_task = cell2mat(cellfun(@(data,baseline,striatum_units) ...
    spikes_norm_smooth_reshape_fcn(data(striatum_units,:,:),baseline(striatum_units)), ...
    ephys.unit_event_psths,sua_baseline,striatum_units,'uni',false));

% Plot heatmap
figure;
colormap(ap.colormap('BWR',[],2));
h = tiledlayout(max(plot_day_grp),size(striatum_sua_task,3),'TileSpacing','compact');
for curr_day_grp = 1:length(plot_day_bins)-1

    % Get units and sorting
    curr_units = find( ...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.tan);

    % (sort max stim, then max within stim)
    sort_idx_cell = cell(3,1);
    [~,max_stim] = max(striatum_sua_tavg,[],2);
    for curr_stim_sort = 1:3
        curr_stim_units = find(max_stim(curr_units)==curr_stim_sort);
        [~,curr_sort_idx] = sort(max(striatum_sua_tavg(curr_units(curr_stim_units),curr_stim_sort),[],2),'descend');
        sort_idx_cell{curr_stim_sort} = curr_stim_units(curr_sort_idx);
    end
    sort_idx = cell2mat(sort_idx_cell);
    plot_units = curr_units(sort_idx);

    % (get rec/IDs of cells for single-unit PSTH plotting)
    curr_sorted_unit_coordinate = ...
        [ephys.animal(striatum_sua_grp.rec(plot_units)), ...
        ephys.rec_day(striatum_sua_grp.rec(plot_units)), ...
        num2cell(striatum_sua_grp.unit_id(plot_units))];

    % (smooth MSN/FSI because large number of units;
    max_n_cells = max(ap.groupfun(@sum, ...
        +(ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        striatum_sua_grp.(curr_celltype)),plot_day_grp));
    smooth_n = 5*max_n_cells/500;

    for curr_stim = 1:n_stim
        nexttile;
        imagesc(psth_t,[],movmean(striatum_sua_task(plot_units,:,curr_stim),[smooth_n,0],1));
        clim([-1,1])
        xlim([-0.2,0.8])
        yline(cumsum(cellfun(@length,sort_idx_cell)),'k');
        axis off;
    end

end
linkaxes(h.Children,'xy');
ap.prettyfig;

% %%%%%%%% TESTING
% curr_unit = 77;
% animal = curr_sorted_unit_coordinate{curr_unit,1};
% rec_day = curr_sorted_unit_coordinate{curr_unit,2};
% plot_units = curr_sorted_unit_coordinate{curr_unit,3};
% AP_longstriatum_psth_fig;
% %%%%%%%%%%%%

% Plot average
figure;
h = tiledlayout(max(plot_day_grp),size(striatum_sua_task,3),'TileSpacing','compact');
stim_colormap = ap.colormap('BKR',n_stim);
for curr_day_grp = 1:length(plot_day_bins)-1

    curr_units = find(...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.tan);

    for curr_align = 1:size(striatum_sua_task,3)
        curr_sua_mean = mean(ap.groupfun(@mean, ...
            striatum_sua_task(curr_units,:,curr_align),striatum_sua_grp.animal(curr_units)),1);
        curr_sua_sem = AP_sem(ap.groupfun(@mean, ...
            striatum_sua_task(curr_units,:,curr_align),striatum_sua_grp.animal(curr_units)),1);
        nexttile;
        ap.errorfill(psth_t,curr_sua_mean,curr_sua_sem,[0.7,0,0]);
    end

end
linkaxes(h.Children,'xy');
ap.prettyfig;
xlim(h.Children(1),[-0.2,0.8]);


% Plot average overlaid across days
align_xlim = {[-0.2,0.15],[-0.05,0.4],[-0.1,0.5]};

figure;
h = tiledlayout(1,size(striatum_sua_task,3),'TileSpacing','compact');
stim_colormap = ap.colormap('BKR',length(plot_day_bins)-1);
for curr_align = 1:size(striatum_sua_task,3)
    nexttile; hold on; set(gca,'ColorOrder',stim_colormap);
    for curr_day_grp = 1:length(plot_day_bins)-1
        curr_units = find(...
            ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            plot_day_grp == curr_day_grp & ...
            striatum_sua_grp.tan);

        % curr_sua_mean = mean(ap.groupfun(@mean, ...
        %     striatum_sua_task(curr_units,:,curr_align),striatum_sua_grp.animal(curr_units)),1);
        % curr_sua_sem = AP_sem(ap.groupfun(@mean, ...
        %     striatum_sua_task(curr_units,:,curr_align),striatum_sua_grp.animal(curr_units)),1);
        % ap.errorfill(psth_t,curr_sua_mean,curr_sua_sem);

        plot(psth_t,median(striatum_sua_task(curr_units,:,curr_align),1),'linewidth',2);
        xlim(align_xlim{curr_align});
    end

end
linkaxes(h.Children,'y');
[h.Children.DataAspectRatio] = deal(min(vertcat(h.Children.DataAspectRatio),[],1));
ap.prettyfig;


% Plot stim rise vs reward dip
stim_t = [0,0.2];
reward_t = [0,0.2];

stim_use_t = isbetween(psth_t,stim_t(1),stim_t(2));
reward_use_t = isbetween(psth_t,reward_t(1),reward_t(2));

stim_act = mean(striatum_sua_task(:,stim_use_t,1),2);
reward_act = mean(striatum_sua_task(:,reward_use_t,3),2);

figure;
h = tiledlayout(1,length(plot_day_bins)-1);
for curr_day_grp = 1:length(plot_day_bins)-1
    % Get units and sorting
    curr_units = find( ...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.tan);

    nexttile; hold on;
    plot(stim_act(curr_units),reward_act(curr_units),'.k','MarkerSize',15);
    xlabel('Stim activity');
    ylabel('Reward activity');
    xline(0);
    yline(0);
end
linkaxes(h.Children,'xy');

figure;
use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & striatum_sua_grp.tan;
stim_act_mean = ap.groupfun(@median,stim_act(use_units),plot_day_grp(use_units));
reward_act_mean = ap.groupfun(@median,reward_act(use_units),plot_day_grp(use_units));
plot(stim_act_mean,reward_act_mean,'linewidth',2);

xline(0);
yline(0);
ap.prettyfig;


use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & striatum_sua_grp.tan;
stim_act_mean = ap.groupfun(@median,stim_act(use_units),plot_day_grp(use_units));
reward_act_mean = ap.groupfun(@median,reward_act(use_units),plot_day_grp(use_units));

figure; hold on;
plot(stim_act_mean,reward_act_mean,'linewidth',2);

figure;
bar([stim_act_mean,reward_act_mean])
legend({'Stim','Reward'})


% mean w/i animal, mean across animals
wi_animal_fcn = @median;

use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & striatum_sua_grp.tan;
[stim_act_mean,grp] = ap.groupfun(wi_animal_fcn,stim_act(use_units),[striatum_sua_grp.animal(use_units),plot_day_grp(use_units)]);
reward_act_mean = ap.groupfun(wi_animal_fcn,reward_act(use_units),[striatum_sua_grp.animal(use_units),plot_day_grp(use_units)]);

figure; hold on;
arrayfun(@(x) plot(stim_act_mean(grp(:,2)==x),reward_act_mean(grp(:,2)==x),'.','MarkerSize',15),1:max(grp(:,2)));

use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & striatum_sua_grp.tan;
stim_act_mean = ap.nestgroupfun({wi_animal_fcn,@mean},stim_act(use_units),striatum_sua_grp.animal(use_units),plot_day_grp(use_units));
reward_act_mean = ap.nestgroupfun({wi_animal_fcn,@mean},reward_act(use_units),striatum_sua_grp.animal(use_units),plot_day_grp(use_units));
stim_act_sem = ap.nestgroupfun({wi_animal_fcn,@AP_sem},stim_act(use_units),striatum_sua_grp.animal(use_units),plot_day_grp(use_units));
reward_act_sem = ap.nestgroupfun({wi_animal_fcn,@AP_sem},reward_act(use_units),striatum_sua_grp.animal(use_units),plot_day_grp(use_units));

errorbar(stim_act_mean,reward_act_mean,stim_act_sem,stim_act_sem,reward_act_sem,reward_act_sem,'linewidth',2);
xline(0);
yline(0);


% plotting reward by grouped stim act
use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & striatum_sua_grp.tan;

figure; hold on;
for curr_daygrp = 1:max(plot_day_grp)
    curr_stim_act = stim_act(use_units & plot_day_grp==curr_daygrp);
    curr_reward_act = reward_act(use_units & plot_day_grp==curr_daygrp);

    curr_stim_bins = discretize(curr_stim_act,prctile(curr_stim_act,linspace(0,100,4)));
    plot(ap.groupfun(@median,curr_stim_act,curr_stim_bins),ap.groupfun(@median,curr_reward_act,curr_stim_bins));
end
xline(0);
yline(0);
xlabel('Stim activity');
ylabel('Reward activity');



%% No-stim TAN rewards

load_dataset = 'task';
loaded_dataset = load_dataset;

data_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Marica_2025\data\nostim';

% Load behavior
load(fullfile(data_path,'bhv'));

% Set "learned_days" and "days_from_learning"
bhv.learned_days = bhv.stimwheel_pval<0.05;
for curr_animal = unique(bhv.animal)'
    curr_rows = strcmp(bhv.animal,curr_animal);
    bhv.days_from_learning(curr_rows) = ...
        (1:sum(curr_rows))' - ...
        max([NaN,find(bhv.learned_days(curr_rows),1)]);
end

%%%%%%%%% run from load manual here

% Set striatal domains to plot (combines if multiple)
plot_domains = 1:2;

% Set days to group
plot_day_bins = [-inf,0,inf];
plot_day_grp = discretize(max(-inf,striatum_sua_grp.ld),plot_day_bins);

% Plot heatmap
figure;
colormap(ap.colormap('BWR',[],2));
h = tiledlayout(max(plot_day_grp),size(striatum_sua,3),'TileSpacing','compact');
for curr_day_grp = 1:length(plot_day_bins)-1

    % Get units and sorting
    curr_units = find( ...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.tan);

    % (sort max stim, then max within stim)
    stim_t = [0.05,0.2];
    stim_use_t = isbetween(psth_t,stim_t(1),stim_t(2));
    striatum_sua_tavg = squeeze(nanmean(striatum_sua(:,stim_use_t,1),2));

    [~,sort_idx] = sort(max(striatum_sua_tavg(curr_units),[],2),'descend');
    plot_units = curr_units(sort_idx);

    % (get rec/IDs of cells for single-unit PSTH plotting)
    curr_sorted_unit_coordinate = ...
        [ephys.animal(striatum_sua_grp.rec(plot_units)), ...
        ephys.rec_day(striatum_sua_grp.rec(plot_units)), ...
        num2cell(striatum_sua_grp.unit_id(plot_units))];

    for curr_align = 1:size(striatum_sua,3)
        nexttile;
        imagesc(psth_t,[],striatum_sua(plot_units,:,curr_align));
        clim([-1,1])
        xlim([-0.2,0.8])
    end

end
linkaxes(h.Children,'xy');
ap.prettyfig;

% Plot average
figure;
h = tiledlayout(max(plot_day_grp),size(striatum_sua,3),'TileSpacing','compact');
for curr_day_grp = 1:length(plot_day_bins)-1

    curr_units = find(...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.tan);

    for curr_align = 1:size(striatum_sua,3)
        curr_sua_mean = mean(ap.groupfun(@mean, ...
            striatum_sua(curr_units,:,curr_align),striatum_sua_grp.animal(curr_units)),1);
        curr_sua_sem = AP_sem(ap.groupfun(@mean, ...
            striatum_sua(curr_units,:,curr_align),striatum_sua_grp.animal(curr_units)),1);
        nexttile;
        ap.errorfill(psth_t,curr_sua_mean,curr_sua_sem,[0.7,0,0]);
    end

end
linkaxes(h.Children,'xy');
ap.prettyfig;
xlim(h.Children(1),[-0.2,0.8]);


% %%%%%%%% TESTING: single cell PSTH
% curr_unit = 10;
% animal = curr_sorted_unit_coordinate{curr_unit,1};
% rec_day = curr_sorted_unit_coordinate{curr_unit,2};
% plot_units = curr_sorted_unit_coordinate{curr_unit,3};
% AP_longstriatum_psth_fig;
%%%%%%%%%%%%


% Plot average overlaid across stim/nostim
align_xlim = {[-0.2,0.15],[-0.05,0.4],[-0.1,0.5],[-0.1,0.5]};

figure;
h = tiledlayout(1,size(striatum_sua,3)/2,'TileSpacing','compact');
for curr_align = 1:size(striatum_sua,3)/2
    nexttile; hold on; set(gca,'ColorOrder',[0,0,0;1,0,0]);

    curr_units = find(...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        plot_day_grp == length(plot_day_bins)-1 & ...
        striatum_sua_grp.tan);

    plot(psth_t,median(striatum_sua(curr_units,:,curr_align),1),'linewidth',2);
    plot(psth_t,median(striatum_sua(curr_units,:,curr_align+size(striatum_sua,3)/2-1),1),'linewidth',2);
    xlim(align_xlim{curr_align});

end
linkaxes(h.Children,'y');
[h.Children.DataAspectRatio] = deal(min(vertcat(h.Children.DataAspectRatio),[],1));
legend({'stim','no-stim'})
ap.prettyfig;


% Plot stim rise vs reward dip
stim_t = [0,0.2];
reward_t = [0,0.2];

stim_use_t = isbetween(psth_t,stim_t(1),stim_t(2));
reward_use_t = isbetween(psth_t,reward_t(1),reward_t(2));

stim_act = squeeze(nanmean(striatum_sua(:,stim_use_t,[1,4]),2));
reward_act = squeeze(nanmean(striatum_sua(:,reward_use_t,[3,6]),2));

figure; hold on;
curr_units = find( ...
    ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    plot_day_grp == 2 & ...
    striatum_sua_grp.tan);

plot(stim_act(curr_units,1),reward_act(curr_units,1),'.k','MarkerSize',15);
plot(stim_act(curr_units,2),reward_act(curr_units,2),'.r','MarkerSize',15);
xlabel('Stim activity');
ylabel('Reward activity');
xline(0);
yline(0);

plot(stim_act(curr_units,:)',reward_act(curr_units,:)','color',[0.5,0.5,0.5]);

figure; hold on;
plot(reward_act(curr_units,:)','k');
axis padded



use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & plot_day_grp==2 & striatum_sua_grp.tan;
stim_act_mean = median(stim_act(use_units,:),1);
reward_act_mean = median(reward_act(use_units,:),1);

figure;
bar([stim_act_mean',reward_act_mean'])
legend({'Stim','Reward'})


%% Fraction responsive units (to match above)

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

% Set days to group
plot_day_bins = [-Inf,-2,0,Inf];
plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

plot_domains = 1:2;

% Grab "responsive" striatal units
unit_responsive_p_thresh = 0.99;
unit_responsive = cell2mat(horzcat(ephys.unit_resp_p_value{:})') > unit_responsive_p_thresh;
striatum_units_responsive = unit_responsive(cell2mat(striatum_units),:);

% Get fraction of responsive unit
celltypes = ["msn","fsi","tan"];

figure;
h = tiledlayout(length(celltypes),1);
stim_colormap = ap.colormap('BKR',3);
for curr_celltype = celltypes
    use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & striatum_sua_grp.(curr_celltype);

    [unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
        +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
        plot_day_grp(use_units));

    unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
        +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
        plot_day_grp(use_units));

    nexttile; hold on;
    set(gca,'ColorOrder',stim_colormap);
    errorbar(unit_responsive_mean,unit_responsive_sem,'linewidth',2);
    ylabel('Frac. responsive units');
    axis padded
    title(curr_celltype);
end
ap.prettyfig;

% Get fraction of responsive unit as barplot
celltypes = ["msn","fsi","tan"];

figure;
h = tiledlayout(1,length(celltypes));
stim_colormap = ap.colormap('BKR',3);
for curr_celltype = celltypes
    use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & striatum_sua_grp.(curr_celltype);

    [unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
        +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
        plot_day_grp(use_units));

    unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
        +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
        plot_day_grp(use_units));

    nexttile; hold on;
    set(gca,'ColorOrder',stim_colormap);
    b = bar(unit_responsive_mean','FaceColor','flat');
    [b.CData] = deal(stim_colormap);

    errorbar(vertcat(b.XEndPoints)',unit_responsive_mean',unit_responsive_sem', ...
        'marker','none','linestyle','none','color','k','linewidth',1);

    ylabel('Frac units');
    axis padded
    title(curr_celltype);

end
ap.prettyfig;


% Get fraction of responsive unit as barplot (by responsive to stim)
celltypes = ["msn","fsi","tan"];

figure;
h = tiledlayout(n_stim,length(celltypes));
stim_colormap = ap.colormap('BKR',3);
for curr_celltype = celltypes
    for curr_stim = 1:n_stim
        use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            striatum_sua_grp.(curr_celltype) & ...
            striatum_units_responsive(:,curr_stim);

        [unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
            +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
            plot_day_grp(use_units));

        unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
            +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
            plot_day_grp(use_units));

        nexttile; hold on;
        set(gca,'ColorOrder',stim_colormap);
        b = bar(unit_responsive_mean','FaceColor','flat');
        [b.CData] = deal(stim_colormap);

        errorbar(vertcat(b.XEndPoints)',unit_responsive_mean',unit_responsive_sem', ...
            'marker','none','linestyle','none','color','k','linewidth',1);

        ylabel('Frac units');
        axis padded
        title(curr_celltype);
    end
end
ap.prettyfig;


%% Cortex movie (task or passive)

plot_day_bins = [-Inf,-2,0,2,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% (task - use rxn > thresh)
use_trials = wf_grp.rxn >= 0.3;
[wf_avg,wf_avg_grp] = ap.nestgroupfun({@mean,@mean}, ...
    cell2mat(cellfun(@(data,use_trials) data(use_trials,:,:,:), ...
    wf.V_event_align,mat2cell(use_trials, ...
    cellfun(@(x) size(x,1),wf.V_event_align)),'uni',false)), ...
    wf_grp.animal(use_trials),plot_day_grp(use_trials));
wf_avg_px = plab.wf.svd2px(U_master,permute(wf_avg,[3,2,1,4]));

% (passive)
[wf_avg,wf_avg_grp] = ap.nestgroupfun({@mean,@mean},cell2mat(wf.V_event_align), ...
    wf_grp.animal,[plot_day_grp,cell2mat(wf.trial_stim_values)]);
wf_avg_px = plab.wf.svd2px(U_master,permute(wf_avg(wf_avg_grp(:,2)==90,:,:),[3,2,1]));

% Plot baseline-subtracted
wf_baseline_t = wf_t < 0;
ap.imscroll(wf_avg_px(:,:,:,:,1) - mean(wf_avg_px(:,:,wf_baseline_t,:,1),3));
colormap(ap.colormap('PWG',[],1.5));
clim(max(abs(clim)).*[-1,1]);
axis image off;
ap.wf_draw('ccf',[0.5,0.5,0.5]);


%% mPFC - troubleshooting per animal

wf_striatum_roi_rec = cellfun(@(act,stim) ap.groupfun(@mean,act,stim), ...
    mat2cell(wf_striatum_roi,cellfun(@length,wf.trial_stim_values)), ...
    wf.trial_stim_values,'uni',false);

figure; tiledlayout;
animals = unique(bhv.animal);
for curr_animal = animals'
    curr_data = cat(4,wf_striatum_roi_rec{strcmp(bhv.animal,curr_animal)});

    nexttile;
    imagesc(permute(curr_data(3,:,2,:),[4,2,1,3]));

    clim([-1,1].*3e-3);
    colormap(ap.colormap('PWG',[],1.5));
    title(curr_animal);
end


% movies for single animal
plot_animal = 'AM015';

wf_avg_rec = cellfun(@(act,stim) ap.groupfun(@mean,act,stim), ...
    wf.V_event_align,wf.trial_stim_values,'uni',false);

wf_animal = plab.wf.svd2px(U_master, ...
    permute(cat(4,wf_avg_rec{strcmp(wf.animal,plot_animal)}),[3,2,4,1]));

ap.imscroll(wf_animal(:,:,:,:,3));
colormap(ap.colormap('PWG',[],1.5));
clim(max(abs(clim)).*[-1,1]);
axis image off;
ap.wf_draw('ccf',[0.5,0.5,0.5]);
set(gcf,'name',plot_animal);


%% Reinhold behavior analysis: P(cue|move)

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
        p_cue_given_move(curr_recording,:) = ap.groupfun(@mean,+move_poststim,move_session_bins);

        % Get association stat
        [rxn_stat_p(curr_recording), ...
            rxn_stat(curr_recording),rxn_null_stat(curr_recording)] = ...
            AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    bhv(curr_animal_idx).move_rate = move_rate;
    bhv(curr_animal_idx).p_cue_given_move = p_cue_given_move;
    bhv(curr_animal_idx).rxn_stat_p = rxn_stat_p;

    % Clear vars except pre-load for next loop
    clearvars('-except',preload_vars{:});
    AP_print_progress_fraction(curr_recording,length(recordings));

end


ld = cellfun(@(x) ((1:size(x,1)) - find(x<0.05,1))',{bhv.rxn_stat_p}','uni',false);

use_animals = ~cellfun(@isempty,ld);

ld_cat = cell2mat(ld(use_animals));

move_rate_cat = cell2mat({bhv(use_animals).move_rate}');
p_c2m_cat = cell2mat({bhv(use_animals).p_cue_given_move}');

ld_split = ld_cat + linspace(0,(n_bins-1)/n_bins,n_bins);

[move_rate_avg,move_rate_grp] = ap.groupfun(@mean,move_rate_cat(:),ld_split(:));
move_rate_sem = ap.groupfun(@AP_sem,move_rate_cat(:),ld_split(:));

[p_c2m_avg,p_c2m_avg_grp] = ap.groupfun(@mean,p_c2m_cat(:),ld_split(:));
p_c2m_sem = ap.groupfun(@AP_sem,p_c2m_cat(:),ld_split(:));

figure; tiledlayout(2,1);

nexttile; 
errorbar( ...
    padarray(reshape(move_rate_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(move_rate_sem,n_bins,[]),[1,0],NaN,'post'),'k','linewidth',2);
xlabel('Learned day')
ylabel(' Move n');
xline(0,'r');

nexttile;
errorbar( ...
    padarray(reshape(p_c2m_avg_grp,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(p_c2m_avg,n_bins,[]),[1,0],NaN,'post'), ...
    padarray(reshape(p_c2m_sem,n_bins,[]),[1,0],NaN,'post'),'k','linewidth',2);
xlabel('Learned day')
ylabel('P(stim|move)');
xline(0,'r');

%% Behavior analysis: wheel velocity per ITI time

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

% Grab learning day for each mouse
bhv = struct;

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};

    use_workflow = 'stim_wheel_right_stage\d';
    recordings = plab.find_recordings(animal,[],use_workflow);

    rxn = nan(length(recordings),1);
    iti_move = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get movement during whole session, and when stim off
        stim_present = interp1([0;photodiode_times],[0;photodiode_values], ...
            timelite.timestamps,'previous','extrap');

        rxn_mean(curr_recording) = mean(stim_to_move);
        move_total(curr_recording) = mean(abs(wheel_velocity));
        move_nostim(curr_recording) = mean(abs(wheel_velocity(~stim_present)));

    end

    bhv(curr_animal_idx).rxn_mean = rxn_mean;
    bhv(curr_animal_idx).move_total = move_total;
    bhv(curr_animal_idx).move_nostim = move_nostim;

    % Clear vars except pre-load for next loop
    clearvars('-except',preload_vars{:});
    AP_print_progress_fraction(curr_recording,length(recordings));

end

rxn_mean_cat = AP_padcatcell(cellfun(@transpose,{bhv.rxn_mean},'uni',false));
move_total_cat = AP_padcatcell(cellfun(@transpose,{bhv.move_total},'uni',false));
move_nostim_cat = AP_padcatcell(cellfun(@transpose,{bhv.move_nostim},'uni',false));

figure; tiledlayout(3,1);
nexttile; errorbar(nanmean(rxn_mean_cat,2),AP_sem(rxn_mean_cat,2),'k','linewidth',2);
ylabel('Reaction');
nexttile; errorbar(nanmean(move_total_cat,2),AP_sem(move_total_cat,2),'k','linewidth',2);
ylabel('Total movement');
nexttile; errorbar(nanmean(move_nostim_cat,2),AP_sem(move_nostim_cat,2),'k','linewidth',2);
ylabel('No-stim movement');

%% Single-unit heatmap (C/R only, all days concat)

% Set striatal domains to plot (combines if multiple)
plot_domains = 1:2;
plot_stim = 2:3;

% Get mean activity in window after stim onset
stim_t = [0,0.2];
stim_use_t = isbetween(psth_t,stim_t(1),stim_t(2));
striatum_sua_tavg = permute(mean(striatum_sua(:,stim_use_t,:),2),[1,3,2]);

% Plot heatmap
celltypes = ["msn","fsi","tan"];

figure;
h = tiledlayout(1,length(celltypes),'TileSpacing','tight');
for curr_celltype = celltypes
    h_sub = tiledlayout(h,1,length(plot_stim),'TileSpacing','tight');
    [~,curr_celltype_idx] = ismember(curr_celltype,celltypes);
    h_sub.Layout.Tile = curr_celltype_idx;

    colormap(ap.colormap('BWR',[],2));
    title(h_sub,curr_celltype);

    % Get units and sorting
    curr_units = find( ...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        any(striatum_units_responsive(:,plot_stim),2) & ...
        striatum_sua_grp.(curr_celltype));

    % (sort max stim, then max within stim)
    sort_idx_cell = cell(length(plot_stim),1);
    [~,max_stim] = max(striatum_sua_tavg(:,plot_stim),[],2);
    for curr_stim_sort = 1:length(plot_stim)
        curr_stim_units = find(max_stim(curr_units)==curr_stim_sort);
        [~,curr_sort_idx] = sort(max(striatum_sua_tavg(curr_units(curr_stim_units), ...
            plot_stim(curr_stim_sort)),[],2),'descend');
        sort_idx_cell{curr_stim_sort} = curr_stim_units(curr_sort_idx);
    end
    sort_idx = cell2mat(sort_idx_cell);
    plot_units = curr_units(sort_idx);

    % (get rec/IDs of cells for single-unit PSTH plotting)
    curr_sorted_unit_coordinate = ...
        [ephys.animal(striatum_sua_grp.rec(plot_units)), ...
        ephys.rec_day(striatum_sua_grp.rec(plot_units)), ...
        num2cell(striatum_sua_grp.unit_id(plot_units))];

    % (y-smooth heatmaps with large number of units)
    smooth_n = 3*length(plot_units)/500;

    for curr_stim = plot_stim
        nexttile(h_sub);
        imagesc(psth_t,[],movmean(striatum_sua(plot_units,:,curr_stim),[smooth_n,0],1));
        clim([-1,1])
        yline(cumsum(cellfun(@length,sort_idx_cell(1:end-1))),'k');
        axis off;
    end

end
linkaxes(vertcat(h.Children.Children),'x');
xlim(vertcat(h.Children.Children),[0,0.5]);
ap.prettyfig;


% % Plot example PSTHs (temp - just uses last run sorting)
% plot_units_n = 5;
% figure;
% h_units = tiledlayout(1,plot_units_n);
% for curr_unit = sum(max_stim(curr_units)==1)+1:sum(max_stim(curr_units)==1)+plot_units_n
%     animal = curr_sorted_unit_coordinate{curr_unit,1};
%     rec_day = curr_sorted_unit_coordinate{curr_unit,2};
%     unit_id = curr_sorted_unit_coordinate{curr_unit,3};
%     AP_longstriatum_psth_fig(animal,rec_day,unit_id,false,h_units);    
% end
% linkaxes(vertcat(h_units.Children.Children),'x');
% ap.prettyfig;


%% Single unit scatter/selectivity plots

plot_celltype = "tan";
plot_domains = 1:2;

% Get max response in window
use_t = psth_t > 0.05 & psth_t < 0.15;
% (t max)
unit_psth_max = permute(max(striatum_sua(:,use_t,:),[],2),[1,3,2]);
% % (t avg)
% unit_psth_max = permute(mean(striatum_sua(:,use_t,:),2),[1,3,2]);
% % (t abs sum)
% % unit_psth_max = permute(sum(abs(striatum_sua(:,use_t,:)),2),[1,3,2]);

plot_day_bins = [-Inf,-2,0,Inf];
unit_plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

% Get responsive units
unit_responsive_p_thresh = 0.95;
unit_responsive = cell2mat(horzcat(ephys.unit_resp_p_value{:})') > unit_responsive_p_thresh;
striatum_units_responsive = unit_responsive(cell2mat(striatum_units),:);


% % %%%%%%%%%% SHUFFLE TEST: C/R TO TEST FOR CHANCE OVERLAP
%
% % shuffle within animal/day group/plot domains/cell type/C or R responsive
% [~,~,shuff_grp] = unique([striatum_sua_grp.animal, ...
%     unit_plot_day_grp, ...
%     ismember(striatum_sua_grp.domain_idx,plot_domains), ...
%     striatum_sua_grp.(plot_celltype), ...
%     any(striatum_units_responsive(:,2:3),2)],'rows');
%
% % striatum_units_responsive = ap.shake(striatum_units_responsive,1,shuff_grp);
% unit_psth_max = ap.shake(unit_psth_max,1,shuff_grp);
%
% %%%%%%%%%%

% Scatter
[~,striatum_units_responsive_category] = ...
    ismember(striatum_units_responsive(:,2:3), ...
    [0,0; ...         % 1 = neither C nor R
    1,0; ...          % 2 = C
    0,1; ...          % 3 = R
    1,1],'rows');     % 4 = C+R
striatum_units_responsive_category_col = ...
    [0.5,0.5,0.5; ...
    0,0,0; ...
    1,0,0; ...
    1,0,1];

figure;
h = tiledlayout(1,length(plot_day_bins)-1,'TileSpacing','tight');
for curr_day_grp = 1:length(plot_day_bins)-1
    curr_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        unit_plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.(plot_celltype);

    nexttile; axis equal
    scatter(unit_psth_max(curr_units,2),unit_psth_max(curr_units,3), ...
        3+20*(striatum_units_responsive_category(curr_units)~=1), ...
        striatum_units_responsive_category(curr_units),'filled');
    line(ylim,ylim);
    clim([1,4]);
end
colormap(striatum_units_responsive_category_col);

ap.prettyfig;
linkaxes(h.Children,'xy');
title(h,plot_celltype);

figure;
h = tiledlayout(1,length(plot_day_bins)-1,'TileSpacing','tight');
for curr_day_grp = 1:length(plot_day_bins)-1
    curr_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        unit_plot_day_grp == curr_day_grp & ...
        any(striatum_units_responsive(:,2:3),2) & ...
        striatum_sua_grp.(plot_celltype);

    nexttile; axis equal
    scatter(unit_psth_max(curr_units,2),unit_psth_max(curr_units,3), ...
        10,[0.5,0.5,0.5],'filled');
    line(ylim,ylim);
    clim([1,4]);
end
colormap(striatum_units_responsive_category_col);

ap.prettyfig;
linkaxes(h.Children,'xy');
title(h,plot_celltype);


% Average response to stimuli by C/R responsiveness
figure;
h = tiledlayout(2,length(plot_day_bins)-1,'TileSpacing','tight');
title(h,plot_celltype);

n_stim = size(striatum_sua,3);
stim_colormap = ap.colormap('BKR',n_stim);
for curr_stim = 2:3
    for curr_day_grp = 1:length(plot_day_bins)-1

        curr_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            unit_plot_day_grp == curr_day_grp & ...
            striatum_units_responsive(:,curr_stim) & ...
            striatum_sua_grp.(plot_celltype);

        curr_sua_mean = permute(mean(ap.groupfun(@mean,striatum_sua(curr_units,:,:), ...
            striatum_sua_grp.animal(curr_units)),1),[2,3,1]);
        curr_sua_sem = permute(AP_sem(ap.groupfun(@mean,striatum_sua(curr_units,:,:), ...
            striatum_sua_grp.animal(curr_units)),1),[2,3,1]);

        nexttile; set(gca,'ColorOrder',stim_colormap);
        ap.errorfill(psth_t,curr_sua_mean,curr_sua_sem);
    end
end
xlim([-0.2,0.8]);
ap.prettyfig;
linkaxes(h.Children,'xy');

% Average response to stimuli by C/R/C+R responsiveness
figure;
h = tiledlayout(3,length(plot_day_bins)-1,'TileSpacing','tight');
title(h,plot_celltype);

n_stim = size(striatum_sua,3);
stim_colormap = ap.colormap('BKR',n_stim);
for curr_responsive_category = 2:4
    for curr_day_grp = 1:length(plot_day_bins)-1

        nexttile; set(gca,'ColorOrder',stim_colormap);

        curr_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            unit_plot_day_grp == curr_day_grp & ...
            striatum_units_responsive_category == curr_responsive_category & ...
            striatum_sua_grp.(plot_celltype);
        if ~any(curr_units)
            continue
        end

        curr_sua_mean = permute(mean(ap.groupfun(@mean,striatum_sua(curr_units,:,:), ...
            striatum_sua_grp.animal(curr_units)),1),[2,3,1]);
        curr_sua_sem = permute(AP_sem(ap.groupfun(@mean,striatum_sua(curr_units,:,:), ...
            striatum_sua_grp.animal(curr_units)),1),[2,3,1]);

        ap.errorfill(psth_t,curr_sua_mean,curr_sua_sem);
    end
end
xlim([-0.2,0.8]);
ap.prettyfig;
linkaxes(h.Children,'xy');

% Average response to stimuli by C/R/C+R responsiveness (all days combined)
figure;
h = tiledlayout(1,3,'TileSpacing','tight');
title(h,plot_celltype);

n_stim = size(striatum_sua,3);
stim_colormap = ap.colormap('BKR',n_stim);
for curr_responsive_category = 2:4

    nexttile; set(gca,'ColorOrder',stim_colormap);

    curr_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        striatum_units_responsive_category == curr_responsive_category & ...
        striatum_sua_grp.(plot_celltype);
    if ~any(curr_units)
        continue
    end

    curr_sua_mean = permute(mean(ap.groupfun(@mean,striatum_sua(curr_units,:,:), ...
        striatum_sua_grp.animal(curr_units)),1),[2,3,1]);
    curr_sua_sem = permute(AP_sem(ap.groupfun(@mean,striatum_sua(curr_units,:,:), ...
        striatum_sua_grp.animal(curr_units)),1),[2,3,1]);

    ap.errorfill(psth_t,curr_sua_mean,curr_sua_sem);
end
xlim([-0.2,0.8]);
ap.prettyfig;
linkaxes(h.Children,'xy');

% Plot proportion overlap
figure;
h = tiledlayout(1,2);

% (as fraction of responsive cells)
curr_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    striatum_sua_grp.(plot_celltype);

frac_responsive = nan(3,length(plot_day_bins)-1);
for curr_day = 1:length(plot_day_bins)-1
    for curr_grp = 2:4

        frac_responsive(curr_grp-1,curr_day) = ...
            sum(curr_units & unit_plot_day_grp == curr_day & ...
            striatum_units_responsive_category == curr_grp)./ ...
            sum(curr_units & unit_plot_day_grp == curr_day & ...
            ismember(striatum_units_responsive_category,2:4));

    end
end
nexttile;
bar(frac_responsive');
legend(["C","R","C+R"]);
ylabel('# this responsive / # any responsive');

% (as fraction of total cells)
curr_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    striatum_sua_grp.(plot_celltype);

frac_responsive = nan(3,length(plot_day_bins)-1);
for curr_day = 1:length(plot_day_bins)-1
    for curr_grp = 2:4

        frac_responsive(curr_grp-1,curr_day) = ...
            sum(curr_units & unit_plot_day_grp == curr_day & ...
            striatum_units_responsive_category == curr_grp)./ ...
            sum(curr_units & unit_plot_day_grp == curr_day);

    end
end
nexttile;
bar(frac_responsive([1,3,2],:)','stacked');
legend(["C","R","C+R"]);
ylabel('# this responsive / # any responsive');


% Fraction of responsive units
figure;
h = tiledlayout(1,3);

% (L/C/R, including double-counts for overlap)
use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    striatum_sua_grp.(plot_celltype);

[unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
    +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
    [unit_plot_day_grp(use_units)]);

unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    +striatum_units_responsive(use_units,:),striatum_sua_grp.animal(use_units), ...
    [unit_plot_day_grp(use_units)]);

nexttile; hold on;
stim_colormap = ap.colormap('BKR',3);
set(gca,'ColorOrder',stim_colormap);

binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');

errorbar(binned_days_x(unit_responsive_mean_group), ...
    unit_responsive_mean,unit_responsive_sem,'linewidth',2);
ylabel('Frac. responsive units');
ap.prettyfig;


% (R/C/R+C, splitting overlap)
use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    striatum_sua_grp.(plot_celltype);

[unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
    +(striatum_units_responsive_category(use_units,:) == [2,3,4]),striatum_sua_grp.animal(use_units), ...
    [unit_plot_day_grp(use_units)]);

unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    +(striatum_units_responsive_category(use_units,:) == [2,3,4]),striatum_sua_grp.animal(use_units), ...
    [unit_plot_day_grp(use_units)]);

nexttile; hold on;
stim_colormap = striatum_units_responsive_category_col(2:end,:);
set(gca,'ColorOrder',stim_colormap);

binned_days_x = interp1(find(~isinf(plot_day_bins)),...
    plot_day_bins(~isinf(plot_day_bins)),1:length(plot_day_bins)-1,'linear','extrap');

errorbar(binned_days_x(unit_responsive_mean_group), ...
    unit_responsive_mean,unit_responsive_sem,'linewidth',2);
ylabel('Frac. responsive units');
ap.prettyfig;


%%%% TESTING: histograms
unit_psth_max_rc = diff(unit_psth_max(:,2:3),[],2);

figure; hold on
for curr_day_grp = 1:length(plot_day_bins)-1
    use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        any(striatum_units_responsive(:,2:3),2) & ...
        unit_plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.(plot_celltype);

    histogram(unit_psth_max_rc(use_units),100,'normalization','probability');
end

figure; hold on
use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    any(striatum_units_responsive(:,2:3),2) & ...
    striatum_sua_grp.(plot_celltype);
histogram(unit_psth_max_rc(use_units),100,'normalization','probability');


unit_psth_max_rc_null = diff(ap.shake(unit_psth_max(:,2:3),2),[],2);
figure; hold on
use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    any(striatum_units_responsive(:,2:3),2) & ...
    striatum_sua_grp.(plot_celltype);
histogram(unit_psth_max_rc_null(use_units),100,'normalization','probability');



% just one stim response
figure;
h = tiledlayout(length(plot_day_bins)-1,1);
for curr_day_grp = 1:length(plot_day_bins)-1
    use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        any(striatum_units_responsive(:,2:3),2) & ...
        unit_plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.(plot_celltype);

    nexttile; hold on;
    histogram(unit_psth_max(use_units,3),linspace(0,10,30),'normalization','probability');
end
linkaxes(h.Children,'xy');

% Frac R/C/R+C as stacked barplot
use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    striatum_sua_grp.(plot_celltype);

[unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
    +(striatum_units_responsive_category(use_units,:) == [2,4,3]),striatum_sua_grp.animal(use_units), ...
    unit_plot_day_grp(use_units));

unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
    +(striatum_units_responsive_category(use_units,:) == [2,4,3]),striatum_sua_grp.animal(use_units), ...
    unit_plot_day_grp(use_units));

figure; hold on;
set(gca,'ColorOrder',striatum_units_responsive_category_col([2,4,3],:));
bar(unit_responsive_mean,'stacked');
ylabel('Frac. responsive units');
legend(["C","C+R","R"]);
ap.prettyfig;

% stats: testing amount of overlap by chance

% shuffle within animal/day group/plot domains/cell type/C or R responsive
[~,~,shuff_grp] = unique([striatum_sua_grp.animal, ...
    unit_plot_day_grp, ...
    ismember(striatum_sua_grp.domain_idx,plot_domains), ...
    striatum_sua_grp.(plot_celltype)],'rows');

responsive_overlap_meas = ap.nestgroupfun({@mean,@mean}, ...
    +(all(striatum_units_responsive(use_units,2:3),2)),striatum_sua_grp.animal(use_units), ...
    unit_plot_day_grp(use_units));

n_shuff = 1000;
responsive_overlap_shuff = nan(length(plot_day_bins)-1,n_shuff);
for curr_shuff = 1:n_shuff

    striatum_units_responsive_shuff = ap.shake(striatum_units_responsive(use_units,:),1,shuff_grp(use_units));
    responsive_overlap_shuff(:,curr_shuff) = ap.nestgroupfun({@mean,@mean}, ...
        +(all(striatum_units_responsive_shuff(:,2:3),2)),striatum_sua_grp.animal(use_units), ...
        unit_plot_day_grp(use_units));
    disp(curr_shuff);

end

stat_rank = tiedrank([responsive_overlap_meas,responsive_overlap_shuff]')';
stat_p = 1-stat_rank(:,1)/(n_shuff+1);


%% R v C unit plots (v2, from above)

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_domains = 1:2;

% Set stim to compare, and colors
compare_stim = [2,3]; % C, R stim
stim_1_color = [0,0,0.8];
stim_2_color = [0.8,0,0];

% Get max response in window
use_t = psth_t > 0.05 & psth_t < 0.15;
unit_psth_max = permute(max(striatum_sua(:,use_t,:),[],2),[1,3,2]);

plot_day_bins = [-Inf,-2,0,Inf];
unit_plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

% Get responsive units
unit_responsive_p_thresh = 0.99;
unit_responsive = cell2mat(horzcat(ephys.unit_resp_p_value{:})') > unit_responsive_p_thresh;
striatum_units_responsive = unit_responsive(cell2mat(striatum_units),:);

% Loop through celltypes, plot scatter and stim responsive overlap
striatum_celltypes = ["msn","fsi","tan"];

figure;
h = tiledlayout(length(striatum_celltypes),1,'TileSpacing','tight');
example_units = nan(length(striatum_celltypes),2);
for curr_celltype = striatum_celltypes

    h_sub = tiledlayout(h,1,length(plot_day_bins),'TileSpacing','tight');
    [~,curr_celltype_idx] = ismember(curr_celltype,striatum_celltypes);
    h_sub.Layout.Tile = curr_celltype_idx;
    title(h_sub,curr_celltype);

    % Scatter
    h_scatter = gobjects(3,1);
    for curr_day_grp = 1:length(plot_day_bins)-1
        curr_units = ...
            ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            unit_plot_day_grp == curr_day_grp & ...
            any(striatum_units_responsive(:,compare_stim),2) & ...
            striatum_sua_grp.(curr_celltype);

        curr_dot_size = 15;
        curr_dot_color = ...
            stim_1_color.*striatum_units_responsive(curr_units,compare_stim(1)) + ...
            stim_2_color.*striatum_units_responsive(curr_units,compare_stim(2));

        h_scatter(curr_day_grp) = nexttile(h_sub); hold on;
        scatter(unit_psth_max(curr_units,2),unit_psth_max(curr_units,3), ...
            curr_dot_size,curr_dot_color,'filled');
    end

    linkaxes(h_scatter,'xy');
    ax_range = [min([xlim,ylim]),max([xlim,ylim])];
    xlim(ax_range);ylim(ax_range);
    axis(h_scatter,'square')
    arrayfun(@(x) line(h_scatter(x),ax_range,ax_range,'color','k'),1:length(h_scatter));

    % Grab example units to plot
    plot_unit_prctile = 90;
    plot_day_grp = 2;
    for curr_stim = compare_stim

        curr_units = ...
            ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            unit_plot_day_grp == plot_day_grp & ...
            striatum_units_responsive(:,curr_stim) & ...
            ~striatum_units_responsive(:,setxor(curr_stim,compare_stim)) & ...
            striatum_sua_grp.(curr_celltype);

        curr_unit_idx = find(curr_units);
        [~,sort_idx] = sort(unit_psth_max(curr_units,curr_stim));
        
        curr_example_unit = curr_unit_idx(sort_idx(ceil(sum(curr_units).*plot_unit_prctile/100)));

        % (circle example units on scatter)
        scatter(h_scatter(plot_day_grp), ...
            unit_psth_max(curr_example_unit,2),unit_psth_max(curr_example_unit,3), ...
            45,[0.5,0.5,0.5],'linewidth',2);

        [~,curr_stim_idx] = ismember(curr_stim,compare_stim);
        example_units(curr_celltype_idx,curr_stim_idx) = curr_example_unit;

    end

    % Frac R/C/R+C as stacked barplot
    use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        striatum_sua_grp.(curr_celltype);

    striatum_units_responsive_1_12_2 = permute(all(striatum_units_responsive(:,compare_stim) == ...
        cat(3,[1,0],[1,1],[0,1]),2),[1,3,2]);

    [unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
        +striatum_units_responsive_1_12_2(use_units,:), ...
        striatum_sua_grp.animal(use_units), ...
        unit_plot_day_grp(use_units));

    unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
        +striatum_units_responsive_1_12_2(use_units,:), ...
        striatum_sua_grp.animal(use_units), ...
        unit_plot_day_grp(use_units));

    nexttile(h_sub); hold on;
    set(gca,'ColorOrder',[stim_1_color;stim_1_color+stim_2_color;stim_2_color]);
    bar(unit_responsive_mean,'stacked');
    ylabel('Frac. responsive units');
    legend(["C","C+R","R"]);

    % % ~~~ STATS ~~~
    % fprintf('---- STATS ----\n')
    % n_shuff = 10000;
    % 
    % % Compare R+C-responsive overlap to shuffling R/C responsiveness
    % use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    %     striatum_sua_grp.(curr_celltype);
    % 
    % [~,~,overlap_shuff_grp] = unique([striatum_sua_grp.animal, ...
    %     unit_plot_day_grp, ...
    %     ismember(striatum_sua_grp.domain_idx,plot_domains), ...
    %     any(striatum_units_responsive(:,compare_stim),2) & ...
    %     striatum_sua_grp.(curr_celltype)],'rows');
    % 
    % responsive_overlap_meas = ap.nestgroupfun({@mean,@mean}, ...
    %     +(all(striatum_units_responsive(use_units,compare_stim),2)), ...
    %     striatum_sua_grp.animal(use_units),unit_plot_day_grp(use_units));
    % 
    % responsive_overlap_shuff = nan(length(plot_day_bins)-1,n_shuff);
    % for curr_shuff = 1:n_shuff
    %     striatum_units_responsive_shuff = ap.shake(striatum_units_responsive(use_units,:),1,overlap_shuff_grp(use_units));
    %     responsive_overlap_shuff(:,curr_shuff) = ap.nestgroupfun({@mean,@mean}, ...
    %         +(all(striatum_units_responsive_shuff(:,compare_stim),2)),striatum_sua_grp.animal(use_units), ...
    %         unit_plot_day_grp(use_units));
    % end
    % 
    % stat_rank = tiedrank([responsive_overlap_meas,responsive_overlap_shuff]')';
    % stat_p = stat_rank(:,1)/(n_shuff+1);
    % 
    % fprintf('Overlap stim %d + %d < shuffle:\n',compare_stim);
    % stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
    % for curr_day = 1:length(plot_day_bins)-1
    %     fprintf('%s day grp %d, p = %.2g%s\n',curr_celltype,curr_day, ...
    %         stat_p(curr_day),stat_sig(curr_day));
    % end
    % 
    % 
    % % Compare R-responsive fraction across days
    % for curr_compare_day = 1:length(plot_day_bins)-2
    % 
    %     compare_day_grps = curr_compare_day+[0,1];
    % 
    %     use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
    %         ismember(unit_plot_day_grp,compare_day_grps) & ...
    %         striatum_sua_grp.(curr_celltype);
    % 
    %     [~,~,r_shuff_grp] = unique([striatum_sua_grp.animal, ...
    %         ismember(striatum_sua_grp.domain_idx,plot_domains), ...
    %         striatum_sua_grp.(curr_celltype)],'rows');
    % 
    %     r_frac_meas = diff(ap.nestgroupfun({@mean,@mean}, ...
    %         +(all(striatum_units_responsive(use_units,3),2)), ...
    %         striatum_sua_grp.animal(use_units),unit_plot_day_grp(use_units)));
    % 
    %     r_frac_shuff = nan(n_shuff,1);
    %     for curr_shuff = 1:n_shuff
    %         unit_plot_day_grp_shuff = ap.shake(unit_plot_day_grp(use_units,:),1,r_shuff_grp(use_units));
    %         r_frac_shuff(curr_shuff) = diff(ap.nestgroupfun({@mean,@mean}, ...
    %             +(all(striatum_units_responsive(use_units,3),2)), ...
    %             striatum_sua_grp.animal(use_units),unit_plot_day_grp_shuff));
    %     end
    % 
    %     stat_rank = tiedrank([r_frac_meas;r_frac_shuff]);
    %     stat_p = 1-stat_rank(1)/(n_shuff+1);
    % 
    %     stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
    %     fprintf('%s R-frac day %d vs %d: p = %.2g%s\n',curr_celltype,compare_day_grps,stat_p,stat_sig);
    % 
    % end
end
ap.prettyfig;


% Plot example PSTHs
figure;
h_units = tiledlayout(1,numel(example_units));
for curr_unit = reshape(example_units',1,[])
    animal = ephys.animal{striatum_sua_grp.rec(curr_unit)};
    rec_day = ephys.rec_day{striatum_sua_grp.rec(curr_unit)};
    unit_id = striatum_sua_grp.unit_id(curr_unit);
    AP_longstriatum_psth_fig(animal,rec_day,unit_id,true,h_units);    
end
linkaxes(vertcat(h_units.Children.Children),'x');
ap.prettyfig;


%% R v C unit plots (v3, from above)
% (scatter from all days, barplots in day groups)

%%% Load data for figure
load_dataset = 'passive';
AP_longstriatum_load_data;
%%%

plot_domains = 1:2;

% Set stim to compare, and colors
compare_stim = [2,3]; % C, R stim
stim_1_color = [0,0,0.8];
stim_2_color = [0.8,0,0];

% Get max response in window
use_t = psth_t > 0.05 & psth_t < 0.15;
unit_psth_max = permute(max(striatum_sua(:,use_t,:),[],2),[1,3,2]);

plot_day_bins = [-Inf,-2,0,Inf];
unit_plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

% Get responsive units
unit_responsive_p_thresh = 0.99;
unit_responsive = cell2mat(horzcat(ephys.unit_resp_p_value{:})') > unit_responsive_p_thresh;
striatum_units_responsive = unit_responsive(cell2mat(striatum_units),:);

% Set celltypes to loop through
striatum_celltypes = ["msn","fsi","tan"];

figure;
h = tiledlayout(length(striatum_celltypes),2, ...
    'TileSpacing','tight','TileIndexing','columnmajor');

% Scatter plots R vs C (days combined)
example_units = nan(length(striatum_celltypes),2);
h_scatter = gobjects(size(striatum_celltypes));
for curr_celltype = striatum_celltypes

    curr_units = ...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        any(striatum_units_responsive(:,compare_stim),2) & ...
        striatum_sua_grp.(curr_celltype);

    curr_dot_size = 15;
    curr_dot_color = ...
        stim_1_color.*striatum_units_responsive(curr_units,compare_stim(1)) + ...
        stim_2_color.*striatum_units_responsive(curr_units,compare_stim(2));

    [~,curr_celltype_idx] = ismember(curr_celltype,striatum_celltypes);
    h_scatter(curr_celltype_idx) = nexttile; hold on;
    scatter(unit_psth_max(curr_units,2),unit_psth_max(curr_units,3), ...
        curr_dot_size,curr_dot_color,'filled');

    ax_range = [min([xlim,ylim]),max([xlim,ylim])];
    xlim(ax_range);ylim(ax_range);
    axis(h_scatter(curr_celltype_idx),'square')
    title(curr_celltype);

end

% Fraction [R,R+C] as stacked barplot over days
for curr_celltype = striatum_celltypes

    use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        striatum_sua_grp.(curr_celltype);

    striatum_units_responsive_2_12 = permute(all(striatum_units_responsive(:,compare_stim) == ...
        cat(3,[0,1],[1,1]),2),[1,3,2]);

    [unit_responsive_mean,unit_responsive_mean_group] = ap.nestgroupfun({@mean,@mean}, ...
        +striatum_units_responsive_2_12(use_units,:), ...
        striatum_sua_grp.animal(use_units), ...
        unit_plot_day_grp(use_units));

    unit_responsive_sem = ap.nestgroupfun({@mean,@AP_sem}, ...
        +striatum_units_responsive_2_12(use_units,:), ...
        striatum_sua_grp.animal(use_units), ...
        unit_plot_day_grp(use_units));

    nexttile; hold on;
    set(gca,'ColorOrder',[stim_2_color;stim_1_color+stim_2_color]);
    bar(unit_responsive_mean,'stacked');
    ylabel('Frac. responsive units');
    legend(["R","C+R"]);
    title(curr_celltype);

end

ap.prettyfig;

% Plot example PSTHs

% (grab example units to plot)
plot_unit_prctile = 80;
for curr_celltype = striatum_celltypes
    for curr_stim = compare_stim
        curr_units = ...
            ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            striatum_units_responsive(:,curr_stim) & ...
            ~striatum_units_responsive(:,setxor(curr_stim,compare_stim)) & ...
            striatum_sua_grp.(curr_celltype);

        curr_unit_idx = find(curr_units);
        [~,sort_idx] = sort(unit_psth_max(curr_units,curr_stim));

        curr_example_unit = curr_unit_idx(sort_idx(ceil(sum(curr_units).*plot_unit_prctile/100)));

        % (circle example units on scatter)
        [~,curr_celltype_idx] = ismember(curr_celltype,striatum_celltypes);

        scatter(h_scatter(curr_celltype_idx), ...
            unit_psth_max(curr_example_unit,2),unit_psth_max(curr_example_unit,3), ...
            45,[0.5,0.5,0.5],'linewidth',2);

        [~,curr_stim_idx] = ismember(curr_stim,compare_stim);
        example_units(curr_celltype_idx,curr_stim_idx) = curr_example_unit;
    end
end

% (load and plot example units)
figure('Name',sprintf('%d%%ile',plot_unit_prctile));
h_units = tiledlayout(1,numel(example_units));
for curr_unit = reshape(example_units',1,[])
    animal = ephys.animal{striatum_sua_grp.rec(curr_unit)};
    rec_day = ephys.rec_day{striatum_sua_grp.rec(curr_unit)};
    unit_id = striatum_sua_grp.unit_id(curr_unit);

    task_flag = true; % (plot task activity)
    quiescence_flag = false; % (plot all trials in passive);
    AP_longstriatum_psth_fig(animal,rec_day,unit_id,task_flag,quiescence_flag,h_units);    
end
xlim(vertcat(h_units.Children.Children),[-0.1,0.5])
ap.prettyfig;

% ~~~ STATS ~~~
fprintf('---- STATS ----\n')
n_shuff = 10000;

for curr_celltype = striatum_celltypes

    % Compare R+C overlap to shuffling R/C responsiveness
    use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        striatum_sua_grp.(curr_celltype);

    [~,~,overlap_shuff_grp] = unique([ ...
        striatum_sua_grp.animal, ...
        unit_plot_day_grp, ...
        striatum_sua_grp.domain_idx, ...
        any(striatum_units_responsive(:,compare_stim),2) & ...
        striatum_sua_grp.(curr_celltype)],'rows');

    responsive_overlap_meas = mean(ap.groupfun(@mean, ...
        +(all(striatum_units_responsive(use_units,compare_stim),2)), ...
        striatum_sua_grp.animal(use_units)));

    responsive_overlap_shuff = nan(n_shuff,1);
    for curr_shuff = 1:n_shuff
        striatum_units_responsive_shuff = ap.shake(striatum_units_responsive(use_units,:),1,overlap_shuff_grp(use_units));
        responsive_overlap_shuff(curr_shuff) = mean(ap.groupfun(@mean, ...
            +(all(striatum_units_responsive_shuff(:,compare_stim),2)), ...
            striatum_sua_grp.animal(use_units)));
    end

    stat_rank = tiedrank([responsive_overlap_meas;responsive_overlap_shuff]);
    stat_p = stat_rank(1)/(n_shuff+1);

    stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
    fprintf('%s R+C overlap vs shuffle: p = %.2g%s\n',curr_celltype, ...
        stat_p,stat_sig);

    % Compare R fraction (R-only and R+C) across days
    for curr_compare_day = 1:length(plot_day_bins)-2

        compare_day_grps = curr_compare_day+[0,1];

        use_units = ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
            ismember(unit_plot_day_grp,compare_day_grps) & ...
            striatum_sua_grp.(curr_celltype);

        [~,~,r_shuff_grp] = unique([striatum_sua_grp.animal, ...
            striatum_sua_grp.domain_idx, ...
            striatum_sua_grp.(curr_celltype)],'rows');

        r_frac_meas = diff(ap.nestgroupfun({@mean,@mean}, ...
            +(all(striatum_units_responsive(use_units,3),2)), ...
            striatum_sua_grp.animal(use_units),unit_plot_day_grp(use_units)));

        r_frac_shuff = nan(n_shuff,1);
        for curr_shuff = 1:n_shuff
            unit_plot_day_grp_shuff = ap.shake(unit_plot_day_grp(use_units,:),1,r_shuff_grp(use_units));
            r_frac_shuff(curr_shuff) = diff(ap.nestgroupfun({@mean,@mean}, ...
                +(all(striatum_units_responsive(use_units,3),2)), ...
                striatum_sua_grp.animal(use_units),unit_plot_day_grp_shuff));
        end

        stat_rank = tiedrank([r_frac_meas;r_frac_shuff]);
        stat_p = 1-stat_rank(1)/(n_shuff+1);

        stat_sig = discretize(stat_p < 0.05,[0,1,Inf],["","*"]);
        fprintf('%s R-frac day %d vs %d: p = %.2g%s\n',curr_celltype,compare_day_grps,stat_p,stat_sig);
    end

end

