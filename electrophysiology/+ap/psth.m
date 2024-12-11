function [psth_avg,psth_trial,psth_t] = psth(event_times,align_times,event_groups,psth_opts)

arguments

    % Discrete data
    event_times
    align_times
    event_groups = ones(size(event_times))

    % PSTH window options
    psth_opts.window (2,1) = [-0.5,1]
    psth_opts.bin_size = 0.001

    % PSTH post-processing options
    psth_opts.smoothing {mustBeNonnegative} = 0
    psth_opts.norm_window (2,1) = [NaN,NaN]
    psth_opts.softnorm {mustBeNonnegative} = 0

end

% [psth_avg,psth_trial,psth_t] = psth(event_times,align_times,event_groups,psth_opts)
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
% psth_avg = group x t x align (averaged across instances within alignment)
% psth_trial = trial x t x align
% psth_t = timepoints for psth (bin centers)

% If 'align_times' is a vector, convert it into a cell
if ~iscell(align_times)
    align_times = {align_times};
end

% Set time bins
t_centers = psth_opts.window(1):psth_opts.bin_size:psth_opts.window(2);
t_bins = [t_centers-(psth_opts.bin_size/2),t_centers(end)+(psth_opts.bin_size/2)];

% Set event groups
event_groups_unique = 1:max(event_groups);

% Get PSTH unit x t x align
% (Use 2D histogram because it's much faster than loop. This requires
% monotonic bins: if there is overlap in windows, split the alignments into
% groups until they're not overlapping)
psth = cell(length(align_times),1);
for curr_align = 1:length(align_times)

    % Get max number of align times that fall within overlapping window
    align_time_diff = align_times{curr_align} - align_times{curr_align}';
    max_t_diff = max(t_bins)-min(t_bins);
    n_split = max(sum(abs(align_time_diff)<max_t_diff));

    % Split align times by the max number of overlaps above
    curr_psth = nan(length(align_times{curr_align}),length(t_centers),length(event_groups_unique));

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

        events_binned_aligned_rate = events_binned_aligned./psth_opts.bin_size;

        curr_psth(curr_align_idx,:,:) = events_binned_aligned_rate;

    end

    % Check for NaNs (should be none if all filled)
    if any(isnan(curr_psth(:)))
        error('PSTH contains NaN: problem in filling all bins')
    end

    % Store full PSTH
    psth{curr_align} = curr_psth;

end

% Post-processing (if selected)

% Smooth
if psth_opts.smoothing > 0
    psth = cellfun(@(x) smoothdata(x,2,'gaussian', ...
        psth_opts.smoothing),psth,'uni',false);
end

% Normalize
if ~all(isnan(psth_opts.norm_window))
    t_baseline = t_centers >= psth_opts.norm_window(1) & ...
        t_centers <= psth_opts.norm_window(2);
    % (compute baseline as average across all alignments)
    psth_baseline = cellfun(@(x) nanmean(x(:,t_baseline,:),[1,2]),psth,'uni',false);
    psth = cellfun(@(x,bl) (x - bl)./(bl + psth_opts.softnorm),psth,psth_baseline,'uni',false);
end

% Set outputs
% (one alignment group: psth = matrix, multiple alignments: psth = cell)
if length(align_times) == 1
    psth_trial = cell2mat(psth);
else
    psth_trial = psth;
end
psth_avg = permute(cell2mat(reshape(cellfun(@(x) nanmean(x,1),psth,'uni',false),[],1)),[3,2,1]);
psth_t = t_centers;






