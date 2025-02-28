
% Find striatum boundaries on probe

%% Find gaps in units

% Based on number of units in a window less than a threshold

unit_gap_window = 200; % window to count units
unit_gap_window_slide = 10; % sliding time for window
unit_gap_thresh = 3; % number of units in window to be called a "gap"

unit_gap_bins = (0:unit_gap_window_slide:(max(channel_positions(:,2))-unit_gap_window)) + ...
    [0;unit_gap_window];
unit_gap_bin_ends = unit_gap_bins(2,:);

unit_bin_counts = cellfun(@(x) histcounts(template_depths,x),num2cell(unit_gap_bins,1));
unit_bin_thresh = unit_bin_counts <= unit_gap_thresh;

unit_gap_ends = unit_gap_bin_ends(diff(unit_bin_thresh) == -1);

% If there's no gap, just call the start of the probe an end of a gap
if isempty(unit_gap_ends)
    unit_gap_ends = min(template_depths)-1;
end

%% Start of striatum: last gap in the top portion of probe

late_gap_threshold = max(channel_positions(:,2))*0.6;

striatum_start = max(unit_gap_ends(unit_gap_ends < late_gap_threshold));

%% End of striatum: end of probe, or before correlated MUA block near tip

% Get correlation of MUA in sliding windows
depth_corr_window = 150; % MUA window in microns
depth_corr_window_spacing = 50; % MUA window spacing in microns

max_depths = 3840; % (hardcode, sometimes kilosort drops channels)

depth_corr_bins = [0:depth_corr_window_spacing:(max_depths-depth_corr_window); ...
    (0:depth_corr_window_spacing:(max_depths-depth_corr_window))+depth_corr_window];
depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

spike_binning_t = 0.1; % seconds
spike_binning_t_edges = nanmin(spike_times_openephys):spike_binning_t:nanmax(spike_times_openephys);

binned_spikes_depth = zeros(size(depth_corr_bins,2),length(spike_binning_t_edges)-1);
for curr_depth = 1:size(depth_corr_bins,2)
    curr_depth_templates_idx = ...
        find(template_depths >= depth_corr_bins(1,curr_depth) & ...
        template_depths < depth_corr_bins(2,curr_depth));

    binned_spikes_depth(curr_depth,:) = histcounts(spike_times_openephys( ...
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
    diff(prctile(max(mean_corr_dim1,mean_corr_dim2,'includenan'),[0,100])) > max(mean_corr_dim1)*0.2

    % If there's a correlated block at end: end of probe is correlated, and
    % there's a dip in correlation between the middle and end
    [~,mua_corr_min_idx] = min(max(mean_corr_dim1,mean_corr_dim2,'includenan'));
    depths_back_centers = depth_corr_bin_centers(end-groups_back+1:end);
    striatum_end = depths_back_centers(mua_corr_min_idx-1);

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

if unit_gap_ends(end) > late_gap_threshold
    striatum_depth = NaN(2,1);
end

%% Get MUA correlation and plot

if false

    % Plot units and MUA correlation with striatum boundaries
    figure('name',sprintf('%s %s',animal,rec_day));
    mua_plot = tiledlayout(1,2);

    ha = nexttile;
    ap.plot_unit_depthrate(spike_templates,template_depths,probe_areas,ha);
    yline(striatum_depth,'r','linewidth',3);

    nexttile;
    imagesc(depth_corr_bin_centers,depth_corr_bin_centers,mua_corr);
    clim([-1,1].*0.5);
    colormap(ap.colormap('BWR'))
    xline(striatum_depth,'g','linewidth',3);
    yline(striatum_depth,'g','linewidth',3);

    linkaxes(mua_plot.Children,'y');

end



