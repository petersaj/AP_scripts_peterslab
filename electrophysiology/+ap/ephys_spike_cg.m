function spike_xcorr = ephys_spike_cg(template_1,template_2,plot_flag)
% spike_xcorr = ephys_spike_cg(template_1,template_2,plot_flag)
%
% Calculate and plot spike Auto/crosscorrelogram

% Grab data from workspace
spike_times_timelite = evalin('base','spike_times_timelite');
spike_templates = evalin('base','spike_templates');

% Parameters
bin_size_expt = 100; % experiment size chunk to get ACG from (s)
max_lag = 500; % max correlation lag (ms)
t_smooth = 50; % gaussian smoothing window size (ms)

if nargin < 3 || isempty(plot_flag)
    plot_flag = false;
end

t_bins_expt = min(spike_times_timelite):bin_size_expt:max(spike_times_timelite);

if nargin == 1 || isempty(template_2) || template_1 == template_2
    % Autocorrelogram

    % Find window with largest number of spikes
    binned_spikes_expt = histcounts(spike_times_timelite(spike_templates == template_1),t_bins_expt);
    [~,max_idx] = max(binned_spikes_expt);
    t_bins = t_bins_expt(max_idx):0.001:t_bins_expt(max_idx+1);
    
    % Bin spikes, autocorrelate, plot
    binned_spikes = histcounts(spike_times_timelite(spike_templates == template_1),t_bins);
    [binned_spikes_xcorr,lags] = xcorr(binned_spikes,max_lag,'coeff');
    binned_spikes_xcorr(lags == 0) = nan;

else
    % Crosscorrelogram

    % Find window with largest number of coincident spikes
    binned_spikes_1_expt = histcounts(spike_times_timelite(spike_templates == template_1),t_bins_expt);
    binned_spikes_2_expt = histcounts(spike_times_timelite(spike_templates == template_2),t_bins_expt);

    [~,max_idx] = max(normalize(binned_spikes_1_expt,'range').*normalize(binned_spikes_2_expt,'range'));
    t_bins = t_bins_expt(max_idx):0.001:t_bins_expt(max_idx+1);

    % Bin spikes, cross-correlate, plot
    binned_spikes_1 = histcounts(spike_times_timelite(spike_templates == template_1),t_bins);
    binned_spikes_2 = histcounts(spike_times_timelite(spike_templates == template_2),t_bins);
    [binned_spikes_xcorr,lags] = xcorr(binned_spikes_1,binned_spikes_2,max_lag,'coeff');

end

binned_spikes_xcorr_smoothed = smoothdata(binned_spikes_xcorr,'gaussian',t_smooth);

if plot_flag
    figure; plot(lags,binned_spikes_xcorr_smoothed,'k');
    xline(0,'color','r');
    xlabel('Lag (ms)');
    title(sprintf('Unit %d, lagged unit %d',template_1,template_2));
end

spike_xcorr = binned_spikes_xcorr_smoothed;







