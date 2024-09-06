function [psth_avg,psth_trial,psth_t] = ephys_psth(spike_times,align_times,spike_groups,psth_opts)

arguments

    % Spike data
    spike_times
    align_times
    spike_groups = ones(size(spike_times))

    % PSTH window options
    psth_opts.window (2,1) = [-0.5,1]
    psth_opts.bin_size = 0.001

    % PSTH post-processing options
    psth_opts.smoothing {mustBeNonnegative} = 0
    psth_opts.norm_window (2,1) = [NaN,NaN]
    psth_opts.softnorm {mustBeNonnegative} = 0

end

% [psth_avg,psth_trial] = ephys_psth(spike_times,align_times,spike_groups,psth_opts)
% 
% Make PSTHs for spike data
%
% INPUTS
% align_times = cell array of times to align
% spike_groups (optional) = grouping variable for spike times (ignore
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
t_bins = psth_opts.window(1):psth_opts.bin_size:psth_opts.window(2);
t_centers = t_bins(1:end-1) + diff(t_bins)./2;

% Set spike groups
spike_groups_unique = 1:max(spike_groups);

% Get PSTH unit x t x align
psth = cell(length(align_times),1);
for curr_align = 1:length(align_times)

    % Get time bins around event (ensure column vector align times)
    align_bins = reshape(align_times{curr_align},[],1) + t_bins;
    align_bins_vector = reshape(align_bins',[],1);

    % Make spike group "bins" (centered on integers)
    spike_group_bins = 0.5:1:length(spike_groups_unique)+1;

    % Set spikes to use (non-NaN group, within time range)
    use_spikes = ...
        spike_times >= min(align_bins_vector) & ...
        spike_times <= max(align_bins_vector) & ...
        ~isnan(spike_groups);

    % Bin spikes 
    spikes_binned_continuous = histcounts2( ...
        spike_times(use_spikes), ...
        spike_groups(use_spikes), ...
        align_bins_vector, ...
        spike_group_bins);

    % Throw away between-event bins and reshape to expected size
    use_continuous_bins = reshape(padarray( ...
        true(size(align_bins(:,1:end-1)')), ...
        [1,0],false,'post'),[],1);

    spikes_binned_aligned = ...
        permute(reshape(spikes_binned_continuous(use_continuous_bins,:), ...
        size(align_bins,2)-1,size(align_bins,1),length(spike_groups_unique)), ...
        [2,1,3]);

    spikes_binned_aligned_rate = spikes_binned_aligned./psth_opts.bin_size;

    psth{curr_align} = spikes_binned_aligned_rate;

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






