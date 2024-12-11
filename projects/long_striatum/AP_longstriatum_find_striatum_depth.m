
% Find striatum boundaries on probe

%% Find gaps in units

unit_gap_thresh = 150;
template_depths_sorted = sort(template_depths);
template_gap = template_depths_sorted(find(diff(template_depths_sorted) > ...
    unit_gap_thresh)+1)-1; % back up 1um

% If there's no gap - just use the first unit
if isempty(template_gap)
    template_gap = min(template_depths)-1;
end

%% Start of striatum: last gap before halfway

probe_halfway_point = 3840/2;

striatum_start = max(template_gap(template_gap < probe_halfway_point));

%% End of striatum: end of probe, or before correlated MUA block near tip

% Get correlation of MUA in sliding windows
depth_corr_window = 150; % MUA window in microns
depth_corr_window_spacing = 50; % MUA window spacing in microns

max_depths = 3840; % (hardcode, sometimes kilosort drops channels)

depth_corr_bins = [0:depth_corr_window_spacing:(max_depths-depth_corr_window); ...
    (0:depth_corr_window_spacing:(max_depths-depth_corr_window))+depth_corr_window];
depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

spike_binning_t = 0.1; % seconds
spike_binning_t_edges = nanmin(spike_times_timelite):spike_binning_t:nanmax(spike_times_timelite);

binned_spikes_depth = zeros(size(depth_corr_bins,2),length(spike_binning_t_edges)-1);
for curr_depth = 1:size(depth_corr_bins,2)
    curr_depth_templates_idx = ...
        find(template_depths >= depth_corr_bins(1,curr_depth) & ...
        template_depths < depth_corr_bins(2,curr_depth));

    binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timelite( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

mua_corr = corrcoef(binned_spikes_depth');

% Look at average correlation backwards from tip on both dimensions
depth_back = 1000;
groups_back = depth_back/depth_corr_window_spacing;
mua_corr_end = medfilt2(mua_corr(end-groups_back+1:end,end-groups_back+1:end),[3,3]);
mua_corr_end(triu(true(length(mua_corr_end)),0)) = nan;
mean_corr_dim1 = nanmean(mua_corr_end,2);
mean_corr_dim2 = nanmean(mua_corr_end,1)';


if max(mean_corr_dim2) > max(mean_corr_dim1)*0.5 && ...
    min(mean_corr_dim1-mean_corr_dim2) < -max(mean_corr_dim1)*0.5

    % If there's a correlated block at end: end of probe is correlated, and
    % end of probe is not correlated with middle of probe
    [~,mua_corr_min_idx] = min(max(mean_corr_dim1,mean_corr_dim2,'includenan'));
    depths_back_centers = depth_corr_bin_centers(end-groups_back+1:end);
    striatum_end = depths_back_centers(mua_corr_min_idx);

else

    % If no correlated block at end, use last unit
    striatum_end = max(template_depths)+1;

end


%% Set striatum depth

if ~isempty(striatum_start) && ~isempty(striatum_end)
    striatum_depth = [striatum_start;striatum_end];
else
    striatum_depth = NaN(2,1);
end


%% Get number of units in striatum

n_striatum_units = sum(template_depths >= striatum_depth(1) & template_depths <= striatum_depth(2));

% if there < N units, exclude this striatal recording
if n_striatum_units < 50
    striatum_depth = NaN(2,1);
end

%% Drop recording: if GPe included

% Compare number of spikes in pre-striatum cortex to striatum. If
% "striatum" has higher-firing units than cortex, it's GPe.
% (only use if there are decent number of cortical units(

n_cortex_units = sum(template_depths >= 0 & template_depths < striatum_depth(1));

if n_cortex_units > 20
    spike_n = accumarray(findgroups(spike_templates),1);

    ctx_v_str_thresh = 1.5;
    gpe_flag = mean(spike_n(template_depths < striatum_start))*ctx_v_str_thresh < ...
        mean(spike_n(template_depths >= striatum_start & template_depths <= striatum_end));

    if gpe_flag
        striatum_depth = NaN(2,1);
    end
end

%% Drop recording: if there's a gap in units more than halfway

if template_gap(end) > probe_halfway_point
    striatum_depth = NaN(2,1);
end

%% Get MUA correlation and plot

if false

    % Plot units and MUA correlation with striatum boundaries
    figure('name',sprintf('%s %s',animal,rec_day));
    h = tiledlayout(1,2);

    ha = nexttile;
    ap.plot_unit_depthrate(spike_templates,template_depths,probe_areas,ha);
    yline(striatum_depth,'r','linewidth',3);

    nexttile;
    imagesc(depth_corr_bin_centers,depth_corr_bin_centers,mua_corr);
    clim([-1,1].*0.5);
    colormap(ap.colormap('BWR'))
    xline(striatum_depth,'g','linewidth',3);
    yline(striatum_depth,'g','linewidth',3);

    linkaxes(h.Children,'y');

end



