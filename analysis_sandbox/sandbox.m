%% SANDBOX
% 
% Temporary code

%% Align widefield to event

% align_times = photodiode_times(1:2:end);
% align_category = vertcat(trial_events.values.TrialStimX);

align_times = reward_times;
align_category = ones(size(align_times));

surround_window = [-0.2,4];

surround_samplerate = 35;
t = surround_window(1):1/surround_samplerate:surround_window(2);
peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

use_U = wf_U;
use_V = wf_V;
use_wf_t = wf_times;

aligned_v = reshape(interp1(use_wf_t,use_V',peri_event_t,'previous'), ...
    length(align_times),length(t),[]);

align_id = findgroups(align_category);
aligned_v_avg = permute(splitapply(@nanmean,aligned_v,align_id),[3,2,1]);
aligned_v_avg_baselined = aligned_v_avg - nanmean(aligned_v_avg(:,t < 0,:),2);

aligned_px_avg = plab.wf.svd2px(use_U,aligned_v_avg_baselined);

AP_imscroll(aligned_px_avg,t);
colormap(AP_colormap('PWG'));
clim(prctile(abs(aligned_px_avg(:)),99.5).*[-1,1]);
axis image;














