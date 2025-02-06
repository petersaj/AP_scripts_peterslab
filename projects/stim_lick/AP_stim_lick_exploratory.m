
%% Lick raster

[lick_psth,lick_raster,lick_t] = ap.psth(lick_times,stimOn_times,...
    'window',[-4,4],'bin_size',0.01,'smoothing',10);

[lick_trial,lick_t_raster_idx] = find(lick_raster);
lick_t_raster = lick_t(lick_t_raster_idx);

figure;
h = tiledlayout(4,1);

nexttile(1)
plot(lick_t,lick_psth,'k');
xline(0,'r');

nexttile([3,1])
plot(lick_t_raster,lick_trial,'.k');
xline(0,'r');
set(gca,'YDir','reverse');

linkaxes(h.Children,'x');


%% Lick raster (stim move conditioning)

% stimOn_times = [photodiode_on_times(1);interp1(photodiode_on_times,photodiode_on_times,reward_times(1:n_trials-1)+1,'next')];
% stim_move_times = [photodiode_off_times(1);interp1(photodiode_off_times,photodiode_off_times,stimOn_times(2:n_trials),'next')];

stimOn_times = photodiode_on_times(1:2:end);
stim_move_times = photodiode_off_times(1:2:end);
stim_center_times = photodiode_off_times(2:2:end);

[lick_psth,lick_raster,lick_t] = ap.psth(lick_times,stimOn_times,...
    'window',[-8,8],'bin_size',0.01,'smoothing',10);

figure;
h = tiledlayout(4,1);

nexttile(1)
plot(lick_t,lick_psth,'k');
xline(0,'r');

nexttile([3,1]);

% stim_static_t = vertcat(trial_events.values(1:n_trials).TrialStimStaticTime);
% [~,sort_idx] = sort(stim_static_t);
% 
% [~,sort_idx] = sort(reward_times - stimOn_times(1:n_trials));

[lick_trial,lick_t_raster_idx] = find(lick_raster(sort_idx,:));
lick_t_raster = lick_t(lick_t_raster_idx);

plot(lick_t_raster,lick_trial,'.k');
hold on;
plot(stim_move_times(sort_idx)-stimOn_times(sort_idx),1:n_trials,'b')
plot(stim_center_times(sort_idx)-stimOn_times(sort_idx),1:n_trials,'c')

xline(0,'r');
set(gca,'YDir','reverse');

linkaxes(h.Children,'x');

% (hacked at the moment)
figure;plot(reward_times - stim_move_times(1:n_trials) - 1,'.k')
ylabel('Stim to reward time')










