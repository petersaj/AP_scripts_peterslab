function spike_xcorr = ephys_spike_cg(template_1,template_2,plot_flag)
% spike_xcorr = ephys_spike_cg(template_1,template_2,plot_flag)
%
% Calculate and plot spike Auto/crosscorrelogram

% Grab data from workspace
spike_times_timelite = evalin('base','spike_times_timelite');
spike_templates = evalin('base','spike_templates');

% Parameters
bin_size_expt = 60*10; % experiment size chunk to get ACG from (s)
max_lag = 1000; % max correlation lag (ms)
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

spike_xcorr = binned_spikes_xcorr_smoothed;

% ALTERNATE METHOD: get ACG from many sections and average
% (didn't look that different, and took longer)

% % Get percent section of recording with the most spikes for this unit
% bin_size_expt = 60*10; % experiment size chunk to get ACG from (s)
% use_recording_percent = 0.5; % fraction of recording to use for stats
% sliding_window_t = 60; % seconds for window sliding to find most spikes
% acg_window = 60*1000; % ms to get each acg (then average)
% 
% recording_times = prctile(spike_times_timelite(spike_templates == template_1),[0,100]);
% acg_window_t = diff(recording_times).*use_recording_percent;
% 
% % (if this spike was not recorded at least for ACG window s, skip)
% if acg_window_t< acg_window/1000
%     spike_xcorr = nan(1,max_lag*2+1);
%     return
% end
% 
% recording_bins = num2cell((recording_times(1):sliding_window_t: ...
%     (recording_times(end)-acg_window_t))+[0;acg_window_t],1);
% 
% curr_spike_counts = cellfun(@(x) histcounts(spike_times_timelite(spike_templates == ...
%     template_1),x),recording_bins);
% 
% [~,max_window] = max(curr_spike_counts);
% recording_window_t = recording_bins{max_window};
% 
% % Within recording section, get ACG in windows
% t_bins = recording_window_t(1):0.001:recording_window_t(2);
% 
% binned_spikes_grouped = reshape(histcounts( ...
%     spike_times_timelite(spike_templates == template_1), ...
%     t_bins(1:end+1-mod(length(t_bins),acg_window))),acg_window,[]);
% 
% binned_spikes_xcorr = cell2mat(cellfun(@(x) xcorr(x,max_lag,'coeff'), ...
%     num2cell(binned_spikes_grouped,1),'uni',false));
% binned_spikes_xcorr(ceil(size(binned_spikes_xcorr,1)/2),:) = NaN;
% 
% % Average and smooth windowed ACG for final
% binned_spikes_xcorr_smoothed = smoothdata(nanmean(binned_spikes_xcorr,2),'gaussian',t_smooth);
% 
% spike_xcorr = binned_spikes_xcorr_smoothed';

% Plot (if flag true)
if plot_flag
    figure; plot(lags,binned_spikes_xcorr_smoothed,'k');
    xline(0,'color','r');
    xlabel('Lag (ms)');
    title(sprintf('Unit %d, lagged unit %d',template_1,template_2));
end















