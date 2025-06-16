%% Notes

% After starting figure code: tests for figures


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
for curr_k = 1:n_k
    figure;
    colormap(ap.colormap('WK',[],2));
    h = tiledlayout(n_stim,max(plot_day_grp),'TileSpacing','compact');
    for curr_stim = 1:n_stim
        for curr_day_grp = 1:length(plot_day_bins)-1
            nexttile;

            curr_units = find(striatum_sua_grp.domain_idx == curr_k & ...
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
    title(h,sprintf('Striatal cluster %d',curr_k));
    linkaxes(h.Children,'xy');
    ap.prettyfig;
end

% Plot scatter and diagonal histograms
stim_t = psth_t > 0.05 & psth_t < 0.15;
compare_stim = [2,3];
scatter_fig = figure;
h_scatter = tiledlayout(scatter_fig,n_k,max(plot_day_grp),'TileSpacing','compact');

histogram_fig = figure;
h_histogram = tiledlayout(histogram_fig,n_k,max(plot_day_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day_grp = 1:length(plot_day_bins)-1
        % (plot scatter)
        nexttile(h_scatter); axis equal; hold on;
         
        curr_units = find(striatum_sua_grp.domain_idx == curr_k & ...
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
h = tiledlayout(n_k,n_stim*max(plot_day_grp),'TileSpacing','compact');
stim_colormap = ap.colormap('BKR',n_stim);
for curr_k = 1:n_k
    for curr_day_grp = 1:length(plot_day_bins)-1

        curr_units = find(striatum_sua_grp.domain_idx == curr_k & ...
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
h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day_grp = 1:max(plot_day_grp)

        curr_units = find(striatum_sua_grp.domain_idx == curr_k & ...
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
h = tiledlayout(n_k,1,'TileSpacing','compact');
for curr_k = 1:n_k

    nexttile; hold on;

    curr_units = find(striatum_sua_grp.domain_idx == curr_k & ...
        plot_day_grp == 1 & celltype_id ~= 0);

    curr_unit_mean = squeeze(mean(striatum_sua(curr_units,stim_t,:),2));

    rc_diff = diff(curr_unit_mean(:,[2,3]),[],2)./ ...
        sum(curr_unit_mean(:,[2,3]),2);

    violinplot(categorical(celltype_order(celltype_id(curr_units))),rc_diff,'DensityDirection','negative');
    rc_diff_med = ap.groupfun(@nanmedian,rc_diff,celltype_id(curr_units));
    plot(categorical(celltype_order),rc_diff_med,'_b','MarkerSize',10);

    curr_units = find(striatum_sua_grp.domain_idx == curr_k & ...
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
h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day_grp = 1:max(plot_day_grp)

        curr_units = find(striatum_sua_grp.domain_idx == curr_k & ...
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
h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day_grp = 1:max(plot_day_grp)

        curr_units = find(striatum_sua_grp.domain_idx == curr_k & ...
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
h = tiledlayout(n_k,max(plot_day_grp),'TileSpacing','compact');
for curr_k = 1:n_k
    for curr_day_grp = 1:max(plot_day_grp)

        curr_units = find(striatum_sua_grp.domain_idx == curr_k & ...
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
plot_domains = 1;

% Set days to group
plot_day_bins = [-inf,-2,0,inf];
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
figure;
colormap(ap.colormap('BWR',[],2));
h = tiledlayout(max(plot_day_grp),n_stim,'TileSpacing','compact');
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

    for curr_stim = 1:n_stim
        nexttile;
        imagesc(psth_t,[],striatum_sua(plot_units,:,curr_stim));
        clim([-1,1])
        xlim([-0.2,0.8])
        yline(cumsum(cellfun(@length,sort_idx_cell)),'k');
    end

end
linkaxes(h.Children,'xy');
ap.prettyfig;

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

%%%% TESTING
% curr_unit = 249+1;
% animal = curr_sorted_unit_coordinate{curr_unit,1};
% rec_day = curr_sorted_unit_coordinate{curr_unit,2};
% plot_units = curr_sorted_unit_coordinate{curr_unit,3};
% AP_longstriatum_psth_fig;


% Plot average
figure;
h = tiledlayout(1,n_stim*max(plot_day_grp),'TileSpacing','compact');
stim_colormap = ap.colormap('BKR',n_stim);
for curr_day_grp = 1:length(plot_day_bins)-1
    curr_units = find(...
        ismember(striatum_sua_grp.domain_idx,plot_domains) & ...
        plot_day_grp == curr_day_grp & ...
        striatum_sua_grp.tan);

    curr_sua_mean = mean(ap.groupfun(@mean, ...
        striatum_sua(curr_units,:,:),striatum_sua_grp.animal(curr_units)),1);
    curr_sua_sem = AP_sem(ap.groupfun(@mean, ...
        striatum_sua(curr_units,:,:),striatum_sua_grp.animal(curr_units)),1);

    for curr_stim = 1:3
        nexttile;
        ap.errorfill(psth_t,curr_sua_mean(:,:,curr_stim),curr_sua_sem(:,:,curr_stim),stim_colormap(curr_stim,:));
    end
end
linkaxes(h.Children,'xy');
ap.prettyfig;
xlim(h.Children(1),[-0.2,0.8]);

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


%% Plot task heatmap with same sorting as above

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
colormap(ap.colormap('BWR',[],1));
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

    for curr_align = 1:size(striatum_sua_task,3)
        nexttile;
        imagesc(psth_t,[],striatum_sua_task(plot_units,:,curr_align));
        clim([-2,2])
        xlim([-0.2,0.8])
        yline(cumsum(cellfun(@length,sort_idx_cell)),'k');
    end

end
linkaxes(h.Children,'xy');
ap.prettyfig;


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



