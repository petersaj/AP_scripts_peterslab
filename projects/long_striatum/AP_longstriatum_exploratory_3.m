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
unit_psth_max = permute(max(striatum_sua(:,use_t,:),[],2),[1,3,2]);

plot_day_bins = [-Inf,-2,0,Inf];
unit_plot_day_grp = discretize(striatum_sua_grp.ld,plot_day_bins);

% Scatter (raw)
figure;
h = tiledlayout(1,length(plot_day_bins)-1);

curr_domain_idx = 1;
for curr_day_grp = 1:length(plot_day_bins)-1
    nexttile; axis equal
    curr_units = ismember([striatum_sua_grp.domain_idx,unit_plot_day_grp],[curr_domain_idx,curr_day_grp],'rows') & ...
        striatum_sua_grp.msn;
    plot(unit_psth_max(curr_units,2),unit_psth_max(curr_units,3),'.k');
end
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
    size_norm_denom = prctile(unit_max_response(striatum_sua_grp.domain_idx == curr_domain_idx),90);
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

plot_day_bins = [-Inf,-1:1,Inf];
plot_day_grp = discretize(max(wf_grp.ld,-inf),plot_day_bins);

% (task)
[wf_avg,wf_avg_grp] = ap.nestgroupfun({@mean,@mean},cell2mat(wf.V_event_align), ...
    wf_grp.animal,plot_day_grp);
wf_avg_px = plab.wf.svd2px(U_master,permute(wf_avg,[3,2,1,4]));

% (passive)
[wf_avg,wf_avg_grp] = ap.nestgroupfun({@mean,@mean},cell2mat(wf.V_event_align), ...
    wf_grp.animal,[plot_day_grp,cell2mat(wf.trial_stim_values)]);
wf_avg_px = plab.wf.svd2px(U_master,permute(wf_avg(wf_avg_grp(:,2)==90,:,:),[3,2,1]));


ap.imscroll(wf_avg_px(:,:,:,:,1));
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







