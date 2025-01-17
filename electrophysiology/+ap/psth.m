function [psth,raster,t] = psth(event_times,align_times,event_groups,opts)

arguments

    % Discrete data
    event_times
    align_times
    event_groups = ones(size(event_times))

    % Raster window options
    opts.window (2,1) = [-0.5,1]
    opts.bin_size = 0.001

    % PSTH post-processing options
    opts.smoothing {mustBeNonnegative} = 0
    opts.norm_window (2,1) = [NaN,NaN]
    opts.softnorm {mustBeNonnegative} = 0

end

% [psth,raster,t] = psth(event_times,align_times,event_groups,opts)
% 
% Make PSTHs for discrete data (e.g. event times)
%
% INPUTS
% event_times = vector of event times
% align_times = vector or cell array of times to align events
% event_groups (optional) = grouping variable for event times (ignore
% NaNs). Currently: must be indicies (e.g. 1,2,3), not categorical.
% 
% Name-Value arguments: 
% window = window around align times to gather (default = [-0.5,1])
% bin_size = bin size (default = 0.001)
% 
% smoothing = gaussian smoothing window size (units in bins) (default = off)
% norm_window = averaging window to define "baseline" for normalizing (default = off)
% softnorm = softening factor for normalizing (default = 0)
%
% OUTPUTS
%
% psth = group x t x align
% raster = trial x t x align
% t = timepoints (bin centers)

% If 'align_times' is a vector, convert it into a cell
if ~iscell(align_times)
    align_times = {align_times};
end

% Set time bins
t = opts.window(1):opts.bin_size:opts.window(2);
t_bins = [t-(opts.bin_size/2),t(end)+(opts.bin_size/2)];

% Set event groups
event_groups_unique = 1:max(event_groups);

% Get raster [group x t x align]
% (Use 2D histogram because it's much faster than loop. This requires
% monotonic bins: if there is overlap in windows, split the alignments into
% groups until they're not overlapping)
raster = cell(length(align_times),1);
for curr_align = 1:length(align_times)

    % Get max number of align times that fall within overlapping window
    align_time_diff = align_times{curr_align} - align_times{curr_align}';
    max_t_diff = max(t_bins)-min(t_bins);
    n_split = max(sum(abs(align_time_diff)<max_t_diff));

    % Split align times by the max number of overlaps above
    curr_raster = nan(length(align_times{curr_align}),length(t),length(event_groups_unique));

    for curr_split = 1:n_split

        % Choose align times (split if necessary)
        curr_align_idx = curr_split:n_split:length(align_times{curr_align});
        curr_align_times = align_times{curr_align}(curr_align_idx);

        % Get time bins around event (ensure column vector align times)
        align_bins = reshape(curr_align_times,[],1) + t_bins;
        align_bins_vector = reshape(align_bins',[],1);

        % Make event group "bins" (centered on integers)
        event_group_bins = 0.5:1:length(event_groups_unique)+1;

        % Set events to use (non-NaN group, within time range)
        use_events = ...
            event_times >= min(align_bins_vector) & ...
            event_times <= max(align_bins_vector) & ...
            ~isnan(event_groups);

        % Bin events
        events_binned_continuous = histcounts2( ...
            event_times(use_events), ...
            event_groups(use_events), ...
            align_bins_vector, ...
            event_group_bins);

        % Throw away between-event bins and reshape to expected size
        use_continuous_bins = reshape(padarray( ...
            true(size(align_bins(:,1:end-1)')), ...
            [1,0],false,'post'),[],1);

        events_binned_aligned = ...
            permute(reshape(events_binned_continuous(use_continuous_bins,:), ...
            size(align_bins,2)-1,size(align_bins,1),length(event_groups_unique)), ...
            [2,1,3]);

        events_binned_aligned_rate = events_binned_aligned./opts.bin_size;

        curr_raster(curr_align_idx,:,:) = events_binned_aligned_rate;

    end

    % Check for NaNs (should be none if all filled)
    if any(isnan(curr_raster(:)))
        error('Raster contains NaN: problem in filling all bins')
    end

    % Store full raster
    raster{curr_align} = curr_raster;

end

% Get PSTH
psth = permute(cell2mat(reshape(cellfun(@(x) nanmean(x,1),raster,'uni',false),[],1)),[3,2,1]);

% PSTH post-processing (if selected)

% Smooth
if opts.smoothing > 0
    psth = smoothdata(psth,2,'gaussian',opts.smoothing);
end

% Normalize
if ~all(isnan(opts.norm_window))
    t_baseline = t >= opts.norm_window(1) & ...
        t <= opts.norm_window(2);
    % (compute baseline as average across all alignments)
    psth_baseline = nanmean(psth(:,t_baseline,:),2);
    psth = (psth - psth_baseline)./(psth_baseline + opts.softnorm);
end

% If one alignment group, set raster as matrix
% (if multiple groups, leave as cell)
if length(align_times) == 1
    raster = cell2mat(raster);
end







