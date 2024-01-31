function deconvolved_activity = wf_deconv(activity,frame_rate)
% deconvolved_activity = wf_deconv(activity,frame_rate)
%
% Deconvolve fluorescence from tetO-GC6s widefield to match MUA spikes
%
% Kernel created from spike/fluorescence (35Hz) recording and calculated in
% AP_deconv_wf_kernelfit
%
% activity = ND array, time in 2nd dim
% sample_rate = sample rate of signal to deconvolve (default = 35Hz)

%% Set default frame rate
if ~exist('frame_rate','var')
    frame_rate = 35;
end

%% Load GCaMP6s widefield kernel

% Use the kernel within the directory with the function
% (tried regression from supersample - didn't change shape)
kernel_path = fileparts(mfilename('fullpath'));
kernel_fn = [kernel_path filesep 'gcamp6s_kernel.mat'];
load(kernel_fn);

kernel_sample_rate = 1/mean(diff(gcamp6s_kernel.regression_t));

% Choose kernel
kernel_cat = vertcat(gcamp6s_kernel.regression{:});

% If lower target frame rate: pad, filter kernel at nyquist of frame rate
if frame_rate < kernel_sample_rate
    kernel_cat_pad = padarray(kernel_cat,[0,size(kernel_cat,2)],0);
    kernel_cat_filt_pad = lowpass(kernel_cat_pad',frame_rate/2,kernel_sample_rate)';
    kernel_cat_filt = kernel_cat_filt_pad(:,size(kernel_cat,2)+1:end-size(kernel_cat,2));
else
    kernel_cat_filt = kernel_cat;
end

% Normalize and average
kernel_mean = nanmean(kernel_cat_filt./vecnorm(kernel_cat_filt,2,2),1);

% Resample kernel to given sampling rate, normalize
% (get range as nearest integer sample)
kernel_resample_t_range = ...
    floor(frame_rate*max(abs(gcamp6s_kernel.regression_t)))/frame_rate;
kernel_resample_t = -kernel_resample_t_range:1/frame_rate:kernel_resample_t_range;
kernel_resample = interp1(gcamp6s_kernel.regression_t,kernel_mean,kernel_resample_t);
kernel = kernel_resample./norm(kernel_resample);

%% Convolve with deconvolution filter
% (zero nans and pad edges)

% Remove NaNs for convolution (put them back after)
activity_nonan = activity;
activity_nonan(isnan(activity_nonan)) = 0;

% Deconvolve to same size with replication padding
deconvolved_activity = convn(padarray(activity_nonan, ...
    [0,floor(length(kernel)/2)],'replicate','both'), ...
    kernel,'valid');

% Put original NaNs back in
deconvolved_activity(isnan(activity)) = NaN;








