
%% Lick raster

[lick_psth,lick_raster] = ap.psth(lick_times,stimOn_times,...
    'window',[-4,4],'bin_size',0.01,'smoothing',10);

[lick_trial,lick_t] = find(lick_raster);

figure;
h = tiledlayout(4,1);

nexttile(1)
plot(lick_psth,'k');

nexttile([3,1])
plot(lick_t,lick_trial,'.k');

linkaxes(h.Children,'x');





