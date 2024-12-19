%% ~~~~~~~~~~~ INDIVIDUAL

%% Get bursts of DMS activity

ap.plot_unit_depthrate(spike_templates,template_depths,probe_areas);

% Get bursts
use_depth = [2000,2500];
yline(use_depth,'r','linewidth',2)

spike_binning_t = 0.005; % seconds
spike_binning_t_edges = (min(timelite.timestamps):spike_binning_t:max(timelite.timestamps))';
spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;

use_spikes = spike_depths > use_depth(1) & spike_depths < use_depth(2);
binned_spikes = histcounts(spike_times_timelite(use_spikes), ...
    spike_binning_t_edges);

binned_spikes_smoothed = smoothdata(binned_spikes,2,'gaussian',0.1/spike_binning_t);

binned_spikes_smoothed_std = binned_spikes_smoothed./std(binned_spikes_smoothed);

burst_thresh = 4;
burst_times = spike_binning_t_centers(find(diff(binned_spikes_smoothed_std > burst_thresh)==1)+1);

% (trim the burst times)
burst_times(burst_times < 30) = [];
burst_times(find(diff(burst_times) < 2)+1) = [];

figure;plot(spike_binning_t_centers,binned_spikes_smoothed_std,'k');
xline(burst_times,'r');

ap.cellraster(burst_times);

% Get widefield aligned to bursts
surround_window = [-0.5,1];
t = surround_window(1):1/wf_framerate:surround_window(2);
peri_event_t = reshape(burst_times,[],1) + reshape(t,1,[]);

aligned_v = reshape(interp1(wf_t,wf_V',peri_event_t,'previous'), ...
    length(burst_times),length(t),[]);

aligned_px_avg = plab.wf.svd2px(wf_U,permute(nanmean(aligned_v - ...
    nanmean(aligned_v(:,t<0,:),2),1),[3,2,1]));

AP_imscroll(aligned_px_avg,t);
colormap(AP_colormap('PWG'));
clim(prctile(abs(aligned_px_avg(:)),100).*[-1,1]);
axis image;
set(gcf,'name',sprintf('%s %s %s',animal,rec_day,bonsai_workflow));



%% Get striatal and cortex ROI trial activity aligned to event

% Set upsample value for regression
upsample_factor = 1;
sample_rate = (1/mean(diff(wf_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% (group multiunit by size bins)

% Find striatum boundaries
AP_longstriatum_find_striatum_depth

depth_size = 400;
depth_group_edges = striatum_depth(1):depth_size:striatum_depth(2);
[depth_group_n,~,depth_group] = histcounts(spike_depths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

ap.plot_unit_depthrate(spike_templates,template_depths,probe_areas);
yline(depth_group_edges,'r','linewidth',2)

% Bin spikes
n_depths = length(depth_group_edges) - 1;
depth_group = discretize(spike_depths,depth_group_edges);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timelite(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

use_svs = 1:100;
kernel_t = [0,0];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 10;
zs = [false,false];
cvfold = 5;
return_constant = false;
use_constant = true;

fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

% Regress cortex to spikes
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample, ...
    binned_spikes_std,kernel_frames, ...
    lambda,zs,cvfold,return_constant,use_constant);

% Convert kernel to pixel space
r_px = plab.wf.svd2px(wf_U(:,:,use_svs),k);

AP_imscroll(r_px,kernel_frames/wf_framerate);
clim([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
colormap(AP_colormap('BWR'));
axis image;


% Regress cortex to spikes
[k,predicted_spikes,explained_var] = ...
    ap.regresskernel(fVdf_deconv_resample, ...
    binned_spikes_std,kernel_frames, ...
    lambda,zs,cvfold,return_constant,use_constant);

% Convert kernel to pixel space
r_px = permute(plab.wf.svd2px(wf_U(:,:,use_svs),k),[1,2,4,3]);
k_roi = r_px > std(r_px(:))*2;
k_roi(:,round(size(k_roi,2)/2):end,:) = false;

% Align data
if contains(bonsai_workflow,'passive')
    % LCR passive: align to quiescent stim onset
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        1:length(stimOn_times))';

    stim_x = vertcat(trial_events.values.TrialStimX);
    align_times = stimOn_times(quiescent_trials & stim_x == 90);
    sort_idx = 1:length(align_times);

elseif contains(bonsai_workflow,'stim_wheel')
    % (task)
    align_times = stimOn_times;
    [~,sort_idx] = sort(stim_to_move);

end

% Align data to align times

psth_sample_rate = 200;

% (striatum)
% (note - missing depth group handled weirdly atm)
[~,str_psth_all,psth_t] = ap.ephys_psth(spike_times_timelite,align_times,depth_group, ...
    'norm_window',[-0.5,0],'smoothing',10,'bin_size',1/psth_sample_rate);
str_psth = nan(size(align_times,1),size(str_psth_all,2),n_depths);
str_psth(:,:,unique(depth_group(~isnan(depth_group)))) = str_psth_all;

% (cortex)
[roi_trace,roi_mask] = ap.wf_roi(wf_U,wf_V,[],[],k_roi);

surround_window = [-0.5,1];
t = surround_window(1):1/psth_sample_rate:surround_window(2);
peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

aligned_trace = reshape(interp1(wf_t,roi_trace',peri_event_t,'previous'), ...
    length(align_times),length(t),[]);
aligned_trace_baselinesub = aligned_trace - mean(aligned_trace(:,t < 0,:),2);

% Plot all
figure;
h = tiledlayout(4,n_depths,'TileIndexing','ColumnMajor');
for curr_str = 1:n_depths
    nexttile;
    imagesc(k_roi(:,:,curr_str));
    axis image off;
    set(gca,'colormap',ap.colormap('WK'));
    ap.wf_draw('ccf','r');

    nexttile;
    imagesc(t,[],aligned_trace_baselinesub(sort_idx,:,curr_str))
    clim(0.8*max(max(aligned_trace_baselinesub,[],'all')).*[-1,1]);
    xline(0);
    set(gca,'colormap',AP_colormap('PWG'))
    title('Cortex');

    nexttile;
    imagesc(psth_t,[],str_psth(sort_idx,:,curr_str))
    clim(0.2*max(max(str_psth,[],'all')).*[-1,1]);
    xline(0);
    set(gca,'colormap',AP_colormap('PWG'))
    title('Striatum');

    aligned_trace_baselinesub_norm = aligned_trace_baselinesub(sort_idx,:,curr_str)./ ...
        std(aligned_trace_baselinesub(:,:,curr_str),[],'all');
    str_psth_norm = str_psth(sort_idx,:,curr_str)./std(str_psth(:,:,curr_str),[],'all');
    nexttile;
    imagesc(aligned_trace_baselinesub_norm - str_psth_norm);
    xline(0);
    clim(max(abs(clim)).*[-1,1]);
    set(gca,'colormap',AP_colormap('BWR'))
    title('Cortex - striatum')

end
title(h,sprintf('%s %s',animal,rec_day))






%% Cortical activity split by striatal activity amount

use_depth = [1000,2000];
ap.plot_unit_depthrate(spike_templates,template_depths,probe_areas);
yline(use_depth,'r');

use_spikes = spike_depths >= use_depth(1) & spike_depths <= use_depth(2);

if contains(bonsai_workflow,'passive')
    % (passive)
    % PSTH for quiescent right-stim trials
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        (1:length(stimOn_times))');

    stim_x = vertcat(trial_events.values.TrialStimX);
    align_times = stimOn_times(stim_x == 90 & quiescent_trials);
elseif contains(bonsai_workflow,'stim_wheel')
    % (task)
    align_times = stimOn_times;
end

[~,psth,psth_t] = ap.ephys_psth(spike_times_timelite(use_spikes),align_times,'smoothing',20);

use_t = [0,0.3];
psth_tavg = nanmean(psth(:,psth_t>=use_t(1) & psth_t<=use_t(2)),2);
[~,sort_idx] = sort(psth_tavg);

n_grps = 3;
trial_grp = discretize(psth_tavg,prctile(psth_tavg,linspace(0,100,n_grps+1)));
psth_grp_avg = ap.groupfun(@nanmean,psth,trial_grp,[]);
figure;
nexttile;
imagesc(psth(sort_idx,:));
nexttile;
colororder(copper(n_grps));
plot(psth_grp_avg');
title('PSTH grouped by activity');

% Widefield by striatal percentile
align_category = trial_grp;
baseline_times = align_times;

surround_window = [-1,4];
baseline_window = [-0.5,-0.1];

surround_samplerate = 35;

t = surround_window(1):1/surround_samplerate:surround_window(2);
baseline_t = baseline_window(1):1/surround_samplerate:baseline_window(2);

peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);
baseline_event_t = reshape(baseline_times,[],1) + reshape(baseline_t,1,[]);

aligned_v = reshape(interp1(wf_t,wf_V',peri_event_t,'previous'), ...
    length(align_times),length(t),[]);
aligned_baseline_v = nanmean(reshape(interp1(wf_t,wf_V',baseline_event_t,'previous'), ...
    length(baseline_times),length(baseline_t),[]),2);

aligned_v_baselinesub = aligned_v - aligned_baseline_v;

align_id = findgroups(reshape(align_category,[],1));

aligned_v_avg = permute(splitapply(@nanmean,aligned_v_baselinesub,align_id),[3,2,1]);
aligned_px_avg = plab.wf.svd2px(wf_U,aligned_v_avg);

AP_imscroll(aligned_px_avg,t);
colormap(AP_colormap('PWG'));
clim(prctile(abs(aligned_px_avg(:)),100).*[-1,1]);
axis image;
set(gcf,'name',sprintf('%s %s %s',animal,rec_day,bonsai_workflow));

% (testing: correlate grouped striatal/cortical activity)
str_grp_tavg = ap.groupfun(@nanmean,psth_tavg,trial_grp,[]);

str_mua_cat = reshape(1-pdist2(str_grp_tavg',reshape(aligned_px_avg,[],n_grps), ...
    'correlation'),size(aligned_px_avg,[1,2,3]));

ap.imscroll(str_mua_cat,t); axis image;
clim([-1,1]);
colormap(AP_colormap('BWR',[],5));

%% MUA task regression

mua_length = 200;

% Find striatum boundaries
AP_longstriatum_find_striatum_depth

% Discretize spikes by depth
depth_group_edges = striatum_depth(1):mua_length:striatum_depth(2);
depth_group_edges(end) = striatum_depth(2);
depth_group = discretize(spike_depths,depth_group_edges);

% Get MUA binned by widefield
sample_rate = (1/mean(diff(wf_t)));
skip_seconds = 60;
time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

binned_spikes = zeros(max(depth_group),length(time_bins)-1);
for curr_depth = 1:max(depth_group)
    curr_spike_times = spike_times_timelite(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

% Task regressors
stim_regressors = histcounts(stimOn_times,time_bins);
reward_regressors = histcounts(reward_times,time_bins);

stim_move_regressors = histcounts(stim_move_time,time_bins);
nonstim_move_times = ...
    setdiff(timelite.timestamps(find(diff(wheel_move) == 1)+1), ...
    stim_move_regressors);
nonstim_move_regressors = histcounts(nonstim_move_times,time_bins);

% Concatenate selected regressors, set parameters
task_regressors = {stim_regressors;reward_regressors;stim_move_regressors;nonstim_move_regressors};
task_regressor_labels = {'Stim','Reward','Stim move','Nonstim move'};

task_t_shifts = { ...
    [-0.2,1]; ... % stim
    [-0.2,1];  ... % outcome
    [-0.2,1];  ... % nonstim move
    [-0.2,1]};    % stim move

task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(wf_framerate)): ...
    round(x(2)*(wf_framerate)),task_t_shifts,'uni',false);
lambda = 0;
zs = [false,false];
cvfold = 5;
use_constant = true;
return_constant = false;

[mua_task_k,mua_task_long,mua_task_expl_var,mua_task_reduced_long] = ...
    AP_regresskernel(task_regressors,binned_spikes,task_regressor_sample_shifts, ...
    lambda,zs,cvfold,return_constant,use_constant);

mua_task_k_permute = cellfun(@(x) permute(x,[3,2,1]),mua_task_k,'uni',false);

%% No-stim task: behavior
% (fixed quiescence, % no stimuli)

animals = {'AP024','AP026'};

% Create master tiled layout
figure;
t = tiledlayout(1,length(animals),'TileSpacing','tight');

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};

    use_workflow = {'stim_wheel_right_fixedquiescence*'};
    recordings = plab.find_recordings(animal,[],use_workflow);

    surround_time = [-5,5];
    surround_sample_rate = 100;
    surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

    n_trials_success = nan(length(recordings),2);
    frac_move_day = nan(length(recordings),1);
    rxn_med = nan(length(recordings),2);
    frac_move_stimalign = nan(length(recordings),length(surround_time_points),2);
    rxn_stat_p = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get total trials/water
        n_trials_success(curr_recording,:) = ...
            [length([trial_events.values.Outcome]), ...
            sum([trial_events.values.Outcome])];

        % Get median stim-outcome time
        % (for stim and no-stim trials)
        n_trials = length([trial_events.values.Outcome]);
        trial_opacity = 1-[trial_events.values(1:n_trials).TrialOpacity]';
        rxn_med(curr_recording,:) = ap.groupfun(@median,stim_to_move,trial_opacity,[]);

        % Get statistical learning (faster to stim than no-stim)
        n_shuff = 1000;
        rxn_diff_null = arrayfun(@(x) diff(ap.groupfun(@median,stim_to_move,ap.shake(trial_opacity),[])),1:n_shuff);
        rxn_diff_real = diff(ap.groupfun(@median,stim_to_move,trial_opacity,[]));
        rxn_stat_rank = tiedrank(horzcat(rxn_diff_real,rxn_diff_null));
        rxn_stat_p = 1-(rxn_stat_rank(1)./(n_shuff+2));

        % Align wheel movement to stim onset
        align_times = stimOn_times(1:n_trials);
        pull_times = align_times + surround_time_points;

        frac_move_day(curr_recording) = nanmean(wheel_move);

        event_aligned_wheel_vel = interp1(timelite.timestamps, ...
            wheel_velocity,pull_times);
        event_aligned_wheel_move = interp1(timelite.timestamps, ...
            +wheel_move,pull_times,'previous');

        frac_move_stimalign(curr_recording,:,:) = ...
            permute(ap.groupfun(@nanmean,event_aligned_wheel_move,trial_opacity,[]),[3,2,1]);

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});
        AP_print_progress_fraction(curr_recording,length(recordings));

    end

    % Draw in tiled layout nested in master
    t_animal = tiledlayout(t,5,1);
    t_animal.Layout.Tile = curr_animal_idx;
    title(t_animal,animal);

    relative_day = days(datetime({recordings.day}) - datetime({recordings(1).day}))+1;
    nonrecorded_day = setdiff(1:length(recordings),relative_day);

    nexttile(t_animal);
    yyaxis left; plot(relative_day,n_trials_success);
    ylabel('# trials');
    yyaxis right; plot(relative_day,frac_move_day);
    ylabel('Fraction time moving');
    xlabel('Day');

    nexttile(t_animal);
    yyaxis left
    plot(relative_day,rxn_med)
    set(gca,'YScale','log');
    ylabel('Med. rxn');
    xlabel('Day');

    yyaxis right
    prestim_max = max(frac_move_stimalign(:,surround_time_points < 0,:),[],2);
    poststim_max = max(frac_move_stimalign(:,surround_time_points > 0,:),[],2);
    stim_move_frac_ratio = (poststim_max-prestim_max)./(poststim_max+prestim_max);
    plot(relative_day,permute(stim_move_frac_ratio,[1,3,2]));
    yline(0);
    ylabel('pre/post move idx');
    xlabel('Day');
    legend({'Stim','No stim'},'location','best')

    nexttile(t_animal); hold on
    set(gca,'ColorOrder',copper(length(recordings)));
    plot(surround_time_points,frac_move_stimalign(:,:,1)','linewidth',2);
    xline(0,'color','k');
    ylabel('Fraction moving');
    xlabel('Time from stim');
    title('Stim');
    ylim([0,1]);

    nexttile(t_animal); hold on
    set(gca,'ColorOrder',copper(length(recordings)));
    plot(surround_time_points,frac_move_stimalign(:,:,2)','linewidth',2);
    xline(0,'color','k');
    ylabel('Fraction moving');
    xlabel('Time from stim');
    title('No stim');
    ylim([0,1]);

    nexttile(t_animal); hold on
    set(gca,'ColorOrder',copper(length(recordings)));
    plot(surround_time_points,(frac_move_stimalign(:,:,1) - ...
        frac_move_stimalign(:,:,2))','linewidth',2);
    xline(0,'color','k');
    ylabel('Fraction moving');
    xlabel('Time from stim');
    title('Stim - no stim');
    ylim([-1,1]);

end





%% No-stim task: cell raster

trial_opacity = [trial_events.values(1:n_trials).TrialOpacity]';

% (regular task)
wheel_starts = timelite.timestamps(diff([0;wheel_move]) == 1);
wheel_stops = timelite.timestamps(diff([0;wheel_move]) == -1);

% (get wheel starts when no stim on screen: not sure this works yet)
iti_move_idx = interp1(photodiode_times, ...
    photodiode_values,wheel_starts,'previous') == 0;
[~,iti_move_sortidx] = sort(wheel_stops(iti_move_idx) - ...
    wheel_starts(iti_move_idx));

[~,rxn_sort_idx] = sort(stim_to_move);

reward_sort_idx = 1:length(reward_times);
[~,opacity_rxn_sort_idx] = sortrows([trial_opacity,stim_to_move]);
trial_sort_idx = [trial_opacity,opacity_rxn_sort_idx,rxn_sort_idx];
rewarded_trials = [trial_events.values(1:n_trials).Outcome] == 1;

ap.cellraster({stimOn_times(1:n_trials),stim_move_time,wheel_starts(iti_move_idx),reward_times_task}, ...
    {trial_sort_idx,trial_sort_idx,iti_move_sortidx,trial_opacity(rewarded_trials)});

set(gcf,'Name',sprintf('%s, %s',animal,rec_day));


%% No-stim task: single day widefield

trial_opacity = [trial_events.values(1:n_trials).TrialOpacity]';

align_times = [ ...
    stimOn_times(trial_opacity == 1); ...
    stimOn_times(trial_opacity == 0); ...
    stim_move_time(trial_opacity == 1); ...
    stim_move_time(trial_opacity == 0) ...
    ];
align_category = vertcat( ...
    1*ones(sum(trial_opacity == 1),1), ...
    2*ones(sum(trial_opacity == 0),1), ...
    3*ones(sum(trial_opacity == 1),1), ...
    4*ones(sum(trial_opacity == 0),1) ...
    );
baseline_times = vertcat(...
    stimOn_times(trial_opacity == 1), ...
    stimOn_times(trial_opacity == 0), ...
    stimOn_times(trial_opacity == 1), ...
    stimOn_times(trial_opacity == 0) ...
    );

surround_window = [-1,4];
baseline_window = [-0.5,-0.1];

surround_samplerate = 35;

t = surround_window(1):1/surround_samplerate:surround_window(2);
baseline_t = baseline_window(1):1/surround_samplerate:baseline_window(2);

peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);
baseline_event_t = reshape(baseline_times,[],1) + reshape(baseline_t,1,[]);

use_U = wf_U;
use_V = wf_V;
use_wf_t = wf_t;

aligned_v = reshape(interp1(use_wf_t,use_V',peri_event_t,'previous'), ...
    length(align_times),length(t),[]);
aligned_baseline_v = nanmean(reshape(interp1(use_wf_t,use_V',baseline_event_t,'previous'), ...
    length(baseline_times),length(baseline_t),[]),2);

aligned_v_baselinesub = aligned_v - aligned_baseline_v;

align_id = findgroups(reshape(align_category,[],1));

aligned_v_avg = permute(splitapply(@nanmean,aligned_v_baselinesub,align_id),[3,2,1]);
aligned_px_avg = plab.wf.svd2px(use_U,aligned_v_avg);

AP_imscroll(aligned_px_avg,t);
colormap(AP_colormap('PWG'));
clim(prctile(abs(aligned_px_avg(:)),100).*[-1,1]);
axis image;
set(gcf,'name',sprintf('%s %s %s',animal,rec_day,bonsai_workflow));



%% ~~~~~~~~~~~ SET BATCH ANIMALS

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

% Excluding:
% AM008, AP009, AP009 - ephys-only
% AM027 - very anterior / in SM-striatum

% Save
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
save_fn = fullfile(data_path,'animals');
save(save_fn,'animals')


%% ~~~~~~~~~~~ BATCH (averaged data)


%% Set animals (fixed quiescence task)

animals = {'AP024','AP026'};

% Save
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data\fixed_quiescence';
save_fn = fullfile(data_path,'animals');
save(save_fn,'animals')

%% Behavior

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

% Create master tiled layout
figure;
t = tiledlayout(1,length(animals),'TileSpacing','tight');

% Grab learning day for each mouse
bhv = struct;

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};

    use_workflow = {'stim_wheel*'};
    recordings = plab.find_recordings(animal,[],use_workflow);

    surround_time = [-5,5];
    surround_sample_rate = 100;
    surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

    n_trials_water = nan(length(recordings),2);
    frac_move_day = nan(length(recordings),1);
    stim_to_move_med = nan(length(recordings),1);
    stim_to_outcome_med = nan(length(recordings),1);
    frac_move_stimalign = nan(length(recordings),length(surround_time_points));
    rxn_stat_p = nan(length(recordings),1);
    rxn_idx =  nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get total trials/water
        n_trials_water(curr_recording,:) = ...
            [length(trial_events.timestamps), ...
            sum(([trial_events.values.Outcome] == 1)*6)];

        % Get median stim-move/outcome time
        stim_to_move_med(curr_recording) = median(stim_to_move);
        stim_to_outcome_med(curr_recording) = ...
            median(seconds([trial_events.timestamps(1:n_trials).Outcome] - ...
            cellfun(@(x) x(1),{trial_events.timestamps(1:n_trials).StimOn})));

        % Align wheel movement to stim onset
        align_times = stimOn_times;
        pull_times = align_times + surround_time_points;

        frac_move_day(curr_recording) = nanmean(wheel_move);

        event_aligned_wheel_vel = interp1(timelite.timestamps, ...
            wheel_velocity,pull_times);
        event_aligned_wheel_move = interp1(timelite.timestamps, ...
            +wheel_move,pull_times,'previous');

        frac_move_stimalign(curr_recording,:) = nanmean(event_aligned_wheel_move,1);

        % Get association stat
        % (skip if only a few trials)
        if n_trials < 10
            continue
        end

        [rxn_stat_p(curr_recording),rxn_med,rxn_null_med] = ...
            AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move);

        % Get "performance index" from 2022 Fig3C
        rxn_idx(curr_recording) = (rxn_med-rxn_null_med)./(rxn_med+rxn_null_med);

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});
        AP_print_progress_fraction(curr_recording,length(recordings));

    end

    % Define learned day from reaction stat p-value and reaction time
    learned_day = rxn_stat_p < 0.05 & stim_to_outcome_med < 5;

    relative_day = days(datetime({recordings.day}) - datetime({recordings(1).day}))+1;
    nonrecorded_day = setdiff(1:length(recordings),relative_day);

    % Draw in tiled layout nested in master
    t_animal = tiledlayout(t,4,1);
    t_animal.Layout.Tile = curr_animal_idx;
    title(t_animal,animal);

    nexttile(t_animal);
    yyaxis left; plot(relative_day,n_trials_water(:,1));
    ylabel('# trials');
    yyaxis right; plot(relative_day,frac_move_day);
    ylabel('Fraction time moving');
    xlabel('Day');
    if any(nonrecorded_day)
        xline(nonrecorded_day,'--k');
    end
    if any(learned_day)
        xline(relative_day(learned_day),'g');
    end

    nexttile(t_animal);
    yyaxis left
    plot(relative_day,stim_to_outcome_med)
    set(gca,'YScale','log');
    ylabel('Med. rxn');
    xlabel('Day');
    if any(nonrecorded_day)
        xline(nonrecorded_day,'--k');
    end
    if any(learned_day)
        xline(relative_day(learned_day),'g');
    end

    yyaxis right
    prestim_max = max(frac_move_stimalign(:,surround_time_points < 0),[],2);
    poststim_max = max(frac_move_stimalign(:,surround_time_points > 0),[],2);
    stim_move_frac_ratio = (poststim_max-prestim_max)./(poststim_max+prestim_max);
    plot(relative_day,stim_move_frac_ratio);
    yline(0);
    ylabel('pre/post move idx');
    xlabel('Day');

    nexttile(t_animal);
    imagesc(surround_time_points,[],frac_move_stimalign); hold on;
    clim([0,1]);
    colormap(gca,AP_colormap('WK'));
    set(gca,'YTick',1:length(recordings),'YTickLabel', ...
        cellfun(@(day,num) sprintf('%d (%s)',num,day(6:end)), ...
        {recordings.day},num2cell(1:length(recordings)),'uni',false));
    xlabel('Time from stim');
    if any(learned_day)
        plot(0,find(learned_day),'.g')
    end

    nexttile(t_animal); hold on
    set(gca,'ColorOrder',copper(length(recordings)));
    plot(surround_time_points,frac_move_stimalign','linewidth',2);
    xline(0,'color','k');
    ylabel('Fraction moving');
    xlabel('Time from stim');
    if any(learned_day)
        AP_errorfill(surround_time_points,frac_move_stimalign(learned_day,:)', ...
            0.02,[0,1,0],0.1,false);

        % Store behavior across animals
        bhv(curr_animal_idx).stim_to_move_med = stim_to_move_med;
        bhv(curr_animal_idx).stim_to_outcome_med = stim_to_outcome_med;
        bhv(curr_animal_idx).rxn_idx = rxn_idx;
        bhv(curr_animal_idx).frac_move_stimalign = frac_move_stimalign;
        bhv(curr_animal_idx).stim_move_frac_ratio = stim_move_frac_ratio;
        bhv(curr_animal_idx).learned_day = learned_day;

    end

    drawnow;

end

% Save
save_fn = fullfile(data_path,'bhv');
save(save_fn,'bhv')
fprintf('Saved %s\n',save_fn);

%% Behavior (fixed quiescence task)

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals_fixedquiescence'));

% Create master tiled layout
figure;
t = tiledlayout(1,length(animals),'TileSpacing','tight');

% Grab learning day for each mouse
bhv = struct;

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};

    use_workflow = {'stim_wheel*'};
    recordings = plab.find_recordings(animal,[],use_workflow);

    surround_time = [-5,5];
    surround_sample_rate = 100;
    surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

    n_trials_success = nan(length(recordings),2);
    frac_move_day = nan(length(recordings),1);
    rxn_med = nan(length(recordings),2);
    frac_move_stimalign = nan(length(recordings),length(surround_time_points),2);
    rxn_stat_p = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get total trials/water
        n_trials_success(curr_recording,:) = ...
            [length([trial_events.values.Outcome]), ...
            sum([trial_events.values.Outcome])];

        % Get median stim-outcome time
        % (for stim and no-stim trials)
        n_trials = length([trial_events.values.Outcome]);
        trial_opacity = 1-[trial_events.values(1:n_trials).TrialOpacity]';
        rxn_med(curr_recording,:) = ap.groupfun(@median,stim_to_move,trial_opacity,[]);

        % Get statistical learning (faster to stim than no-stim)
        n_shuff = 10000;
        rxn_diff_null = arrayfun(@(x) diff(ap.groupfun(@median,stim_to_move,ap.shake(trial_opacity),[])),1:n_shuff);
        rxn_diff_real = diff(ap.groupfun(@median,stim_to_move,trial_opacity,[]));
        rxn_stat_rank = tiedrank(horzcat(rxn_diff_real,rxn_diff_null));
        rxn_stat_p(curr_recording) = 1-(rxn_stat_rank(1)./(n_shuff+2));

        % Align wheel movement to stim onset
        align_times = stimOn_times(1:n_trials);
        pull_times = align_times + surround_time_points;

        frac_move_day(curr_recording) = nanmean(wheel_move);

        event_aligned_wheel_vel = interp1(timelite.timestamps, ...
            wheel_velocity,pull_times);
        event_aligned_wheel_move = interp1(timelite.timestamps, ...
            +wheel_move,pull_times,'previous');

        frac_move_stimalign(curr_recording,:,:) = ...
            permute(ap.groupfun(@nanmean,event_aligned_wheel_move,trial_opacity,[]),[3,2,1]);

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});
        AP_print_progress_fraction(curr_recording,length(recordings));

    end

    % Define learned day from reaction stat p-value and reaction time
    learned_day = rxn_stat_p < 0.05;

    % Draw in tiled layout nested in master
    t_animal = tiledlayout(t,5,1);
    t_animal.Layout.Tile = curr_animal_idx;
    title(t_animal,animal);

    relative_day = days(datetime({recordings.day}) - datetime({recordings(1).day}))+1;
    nonrecorded_day = setdiff(1:length(recordings),relative_day);

    nexttile(t_animal);
    yyaxis left; plot(relative_day,n_trials_success);
    ylabel('# trials');
    yyaxis right; plot(relative_day,frac_move_day);
    ylabel('Fraction time moving');
    xlabel('Day');

    nexttile(t_animal);
    yyaxis left
    plot(relative_day,rxn_med)
    set(gca,'YScale','log');
    ylabel('Med. rxn');
    xlabel('Day');

    yyaxis right
    prestim_max = max(frac_move_stimalign(:,surround_time_points < 0,:),[],2);
    poststim_max = max(frac_move_stimalign(:,surround_time_points > 0,:),[],2);
    stim_move_frac_ratio = (poststim_max-prestim_max)./(poststim_max+prestim_max);
    plot(relative_day,permute(stim_move_frac_ratio,[1,3,2]));
    yline(0);
    ylabel('pre/post move idx');
    xlabel('Day');
    legend({'Stim','No stim'},'location','best')

    nexttile(t_animal); hold on
    set(gca,'ColorOrder',copper(length(recordings)));
    plot(surround_time_points,frac_move_stimalign(:,:,1)','linewidth',2);
    xline(0,'color','k');
    ylabel('Fraction moving');
    xlabel('Time from stim');
    title('Stim');
    ylim([0,1]);

    nexttile(t_animal); hold on
    set(gca,'ColorOrder',copper(length(recordings)));
    plot(surround_time_points,frac_move_stimalign(:,:,2)','linewidth',2);
    xline(0,'color','k');
    ylabel('Fraction moving');
    xlabel('Time from stim');
    title('No stim');
    ylim([0,1]);

    nexttile(t_animal); hold on
    set(gca,'ColorOrder',copper(length(recordings)));
    plot(surround_time_points,(frac_move_stimalign(:,:,1) - ...
        frac_move_stimalign(:,:,2))','linewidth',2);
    xline(0,'color','k');
    ylabel('Fraction moving');
    xlabel('Time from stim');
    title('Stim - no stim');
    ylim([-1,1]);

    % Store behavior across animals
    bhv(curr_animal_idx).learned_day = learned_day;

    drawnow;

end

% Save
save_fn = fullfile(data_path,'bhv');
save(save_fn,'bhv')
fprintf('Saved %s\n',save_fn);


%% Widefield passive (master U)

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.03;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

day_V_all = cell(length(animals),1);

for curr_animal = 1:length(animals)

    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    n_components = 2000;
    n_align = 3; % (just hardcoded for now)
    day_V = nan(n_components,length(t_centers),n_align,length(recordings));

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % If widefield isn't aligned, skip
        if ~load_parts.widefield_align
            continue
        end

        % Stim times (quiescent trials only)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        stim_x = vertcat(trial_events.values.TrialStimX);
        % sometimes not all trials shown?
        stim_x = stim_x(1:length(stimOn_times));
        align_times = cellfun(@(x) stimOn_times(stim_x == x & quiescent_trials), ...
            num2cell(unique(stim_x)),'uni',false);

        % Align ROI trace to align times
        aligned_V_mean = nan(n_components,length(t_centers),length(align_times));
        for curr_align = 1:length(align_times)
            peri_event_t = align_times{curr_align} + t_centers;

            aligned_V = interp1(wf_t,wf_V',peri_event_t);

            aligned_V_baselinesub = aligned_V - ...
                mean(aligned_V(:,t_centers < 0,:),2);

            aligned_V_mean(:,:,curr_align) = ...
                permute(nanmean(aligned_V_baselinesub,1),[3,2,1]);
        end

        % Store mean aligned ROI trace
        day_V(:,:,:,curr_recording) = aligned_V_mean;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_V_all{curr_animal} = day_V;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Save
save_fn = fullfile(data_path,'wf_passive');
save(save_fn,'day_V_all')
fprintf('Saved %s\n',save_fn);

%% Widefield task (master U)

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.03;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

day_V_all = cell(length(animals),1);

for curr_animal = 1:length(animals)

    animal = animals{curr_animal};

    use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    n_components = 2000;
    n_align = 3; % (just hardcoded for now)
    day_V = nan(n_components,length(t_centers),n_align,length(recordings));

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % If widefield isn't aligned, skip
        if ~load_parts.widefield_align
            continue
        end

        % Stim times (task)
        align_times = {stimOn_times,stim_move_time,reward_times};

        % Align ROI trace to align times
        aligned_V_mean = nan(n_components,length(t_centers),length(align_times));
        for curr_align = 1:length(align_times)
            peri_event_t = align_times{curr_align} + t_centers;

            aligned_V = interp1(wf_t,wf_V',peri_event_t);

            aligned_V_baselinesub = aligned_V - ...
                mean(aligned_V(:,t_centers < 0,:),2);

            aligned_V_mean(:,:,curr_align) = ...
                permute(nanmean(aligned_V_baselinesub,1),[3,2,1]);
        end

        % Store mean aligned ROI trace
        day_V(:,:,:,curr_recording) = aligned_V_mean;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_V_all{curr_animal} = day_V;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Save
save_fn = fullfile(data_path,'wf_task');
save(save_fn,'day_V_all')
fprintf('Saved %s\n',save_fn);


%% MUA passive by depth

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

plot_depths = false;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set MUA length (microns)
mua_length = 200;

day_mua_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_mua = cell(length(recordings),1);

    if plot_depths
        figure;
        h = tiledlayout(1,length(recordings));
        title(h,animal);
    end

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Find striatum boundaries (skip if no marked striatum)
        AP_longstriatum_find_striatum_depth
        if any(isnan(striatum_depth))
            continue
        end

        % Discretize spikes by depth
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(2);
        depth_group_edges(end) = striatum_depth(2);
        spike_depth_group = discretize(spike_depths,depth_group_edges);
        n_depths = max(spike_depth_group);

        % Plot units and depth groupings
        if plot_depths
            nexttile(h);
            set(gca,'YDir','reverse'); hold on
            norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
            unit_dots = scatter(norm_spike_n,template_depths( ...
                unique(spike_templates)),20,'k','filled');
            yline(depth_group_edges,'color','r');
            drawnow;
        end

        % Stim times (quiescent trials only)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        if isfield(trial_events.values,'TrialStimX')
            stim_x = vertcat(trial_events.values.TrialStimX);
            % sometimes not all trials shown?
            stim_x = stim_x(1:length(stimOn_times));
            align_times = cellfun(@(x) stimOn_times(stim_x == x & quiescent_trials), ...
                num2cell(unique(stim_x)),'uni',false);
        else
            % Stim times (task)
            % (skip first trial - somes very short iti?)
            % (also skip stim_to_move negative - bad quiescence))
            use_trials = true(n_trials,1);
            use_trials(1) = false;
            use_trials(stim_to_move < 0) = false;

            align_times = {stimOn_times(use_trials),stim_move_time(use_trials),reward_times_task};
        end

        % Get PSTH depth x t x align
        depth_psth = nan(n_depths,length(t_bins)-1,length(align_times));
        for curr_align = 1:length(align_times)

            t_peri_event = align_times{curr_align} + t_bins;

            use_spikes = spike_times_timelite >= min(t_peri_event,[],'all') & ...
                spike_times_timelite <= max(t_peri_event,[],'all') & ...
                ~isnan(spike_depth_group);

            spikes_binned_continuous = histcounts2(spike_times_timelite(use_spikes),spike_depth_group(use_spikes), ...
                reshape(t_peri_event',[],1),1:n_depths+1)./psth_bin_size;

            use_continuous_bins = reshape(padarray(true(size(t_peri_event(:,1:end-1)')),[1,0],false,'post'),[],1);
            spikes_binned = spikes_binned_continuous(use_continuous_bins,:);

            depth_psth(:,:,curr_align) = ...
                nanmean(permute(reshape(spikes_binned, ...
                size(t_peri_event,2)-1,[],n_depths),[3,1,2]),3);
        end

        smooth_size = 100;
        depth_psth_smooth = smoothdata(depth_psth,2,'gaussian',smooth_size);

        % Store MUA (raw - to combine and normalize later)
        day_mua{curr_recording} = depth_psth_smooth;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_mua_all{curr_animal} = day_mua;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Save
save_fn = fullfile(data_path,'mua_passive');
save(save_fn,'day_mua_all')
fprintf('Saved %s\n',save_fn);

%% MUA passive by depth (excluding TANs)

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

plot_depths = false;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set MUA length (microns)
mua_length = 200;

day_mua_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
    %     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_mua = cell(length(recordings),1);

    if plot_depths
        figure;
        h = tiledlayout(1,length(recordings));
        title(h,animal);
    end

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Find striatum boundaries (skip if no marked striatum)
        AP_longstriatum_find_striatum_depth
        if any(isnan(striatum_depth))
            continue
        end

        % Discretize spikes by depth
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(2);
        depth_group_edges(end) = striatum_depth(2);
        spike_depth_group = discretize(spike_depths,depth_group_edges);
        n_depths = max(spike_depth_group);

        % Plot units and depth groupings
        if plot_depths
            nexttile(h);
            set(gca,'YDir','reverse'); hold on
            norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
            unit_dots = scatter(norm_spike_n,template_depths( ...
                unique(spike_templates)),20,'k','filled');
            yline(depth_group_edges,'color','r');
            drawnow;
        end

        % Stim times (quiescent trials only)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        if isfield(trial_events.values,'TrialStimX')
            stim_x = vertcat(trial_events.values.TrialStimX);
            % sometimes not all trials shown?
            stim_x = stim_x(1:length(stimOn_times));
            align_times = cellfun(@(x) stimOn_times(stim_x == x & quiescent_trials), ...
                num2cell(unique(stim_x)),'uni',false);
        else
            % Stim times (task)
            % (skip first trial - somes very short iti?)
            % (also skip stim_to_move negative - bad quiescence))
            use_trials = true(n_trials,1);
            use_trials(1) = false;
            use_trials(stim_to_move < 0) = false;

            align_times = {stimOn_times(use_trials),stim_move_time(use_trials),reward_times_task};
        end

         % Classify tans (by spike rate and ACG time to steady-state)
        spike_acg = cell2mat(arrayfun(@(x) ap.ephys_spike_cg(x),(1:size(waveforms,1))','uni',false));
        acg_isi = arrayfun(@(x) ...
            find(spike_acg(x,ceil(size(spike_acg,2)/2):end) > ...
            mean(spike_acg(x,end-100:end),2),1,'first'),(1:size(templates,1))');
        spike_rate = accumarray(spike_templates,1)/diff(prctile(spike_times_timelite,[0,100]));
        tans = find(spike_rate >= 4 & spike_rate <= 12 & acg_isi >= 40);

        % Get PSTH depth x t x align
        depth_psth = nan(n_depths,length(t_bins)-1,length(align_times));
        for curr_align = 1:length(align_times)

            t_peri_event = align_times{curr_align} + t_bins;

            use_spikes = spike_times_timelite >= min(t_peri_event,[],'all') & ...
                spike_times_timelite <= max(t_peri_event,[],'all') & ...
                ~isnan(spike_depth_group) & ...
                ~ismember(spike_templates,tans); % Exclude TANs

            spikes_binned_continuous = histcounts2(spike_times_timelite(use_spikes),spike_depth_group(use_spikes), ...
                reshape(t_peri_event',[],1),1:n_depths+1)./psth_bin_size;

            use_continuous_bins = reshape(padarray(true(size(t_peri_event(:,1:end-1)')),[1,0],false,'post'),[],1);
            spikes_binned = spikes_binned_continuous(use_continuous_bins,:);

            depth_psth(:,:,curr_align) = ...
                nanmean(permute(reshape(spikes_binned, ...
                size(t_peri_event,2)-1,[],n_depths),[3,1,2]),3);
        end

        smooth_size = 100;
        depth_psth_smooth = smoothdata(depth_psth,2,'gaussian',smooth_size);

        % Store MUA (raw - to combine and normalize later)
        day_mua{curr_recording} = depth_psth_smooth;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_mua_all{curr_animal} = day_mua;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Save
save_fn = fullfile(data_path,'mua_passive_notan');
save(save_fn,'day_mua_all')
fprintf('Saved %s\n',save_fn);


%% MUA task by depth

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

plot_depths = false;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set MUA length (microns)
mua_length = 200;

day_mua_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_mua = cell(length(recordings),1);

    if plot_depths
        figure;
        h = tiledlayout(1,length(recordings));
        title(h,animal);
    end

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Find striatum boundaries (skip if no marked striatum)
        AP_longstriatum_find_striatum_depth
        if any(isnan(striatum_depth))
            continue
        end

        % Discretize spikes by depth
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(2);
        depth_group_edges(end) = striatum_depth(2);
        spike_depth_group = discretize(spike_depths,depth_group_edges);
        n_depths = max(spike_depth_group);

        % Plot units and depth groupings
        if plot_depths
            nexttile(h);
            set(gca,'YDir','reverse'); hold on
            norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
            unit_dots = scatter(norm_spike_n,template_depths( ...
                unique(spike_templates)),20,'k','filled');
            yline(depth_group_edges,'color','r');
            drawnow;
        end

        % Stim times (quiescent trials only)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        if isfield(trial_events.values,'TrialStimX')
            stim_x = vertcat(trial_events.values.TrialStimX);
            % sometimes not all trials shown?
            stim_x = stim_x(1:length(stimOn_times));
            align_times = cellfun(@(x) stimOn_times(stim_x == x & quiescent_trials), ...
                num2cell(unique(stim_x)),'uni',false);
        else
            % Stim times (task)
            % (skip first trial - somes very short iti?)
            % (also skip stim_to_move negative - bad quiescence))
            use_trials = true(n_trials,1);
            use_trials(1) = false;
            use_trials(stim_to_move < 0) = false;

            align_times = {stimOn_times(use_trials),stim_move_time(use_trials),reward_times_task};
        end

        % Get PSTH depth x t x align
        depth_psth = nan(n_depths,length(t_bins)-1,length(align_times));
        for curr_align = 1:length(align_times)

            t_peri_event = align_times{curr_align} + t_bins;

            use_spikes = spike_times_timelite >= min(t_peri_event,[],'all') & ...
                spike_times_timelite <= max(t_peri_event,[],'all') & ...
                ~isnan(spike_depth_group);

            spikes_binned_continuous = histcounts2(spike_times_timelite(use_spikes),spike_depth_group(use_spikes), ...
                reshape(t_peri_event',[],1),1:n_depths+1)./psth_bin_size;

            use_continuous_bins = reshape(padarray(true(size(t_peri_event(:,1:end-1)')),[1,0],false,'post'),[],1);
            spikes_binned = spikes_binned_continuous(use_continuous_bins,:);

            depth_psth(:,:,curr_align) = ...
                nanmean(permute(reshape(spikes_binned, ...
                size(t_peri_event,2)-1,[],n_depths),[3,1,2]),3);
        end

        smooth_size = 100;
        depth_psth_smooth = smoothdata(depth_psth,2,'gaussian',smooth_size);

        % Store MUA (raw - to combine and normalize later)
        day_mua{curr_recording} = depth_psth_smooth;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_mua_all{curr_animal} = day_mua;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Save
save_fn = fullfile(data_path,'mua_task');
save(save_fn,'day_mua_all')
fprintf('Saved %s\n',save_fn);

%% MUA task regression by depth

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

plot_depths = false;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set MUA length (microns)
mua_length = 200;

mua_task_k_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    %     use_workflow = 'lcr_passive';
    use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_mua_task_k = cell(length(recordings),1);

    if plot_depths
        figure;
        h = tiledlayout(1,length(recordings));
        title(h,animal);
    end

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Find striatum boundaries (skip if no marked striatum)
        AP_longstriatum_find_striatum_depth
        if any(isnan(striatum_depth))
            continue
        end

        % Discretize spikes by depth
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(2);
        depth_group_edges(end) = striatum_depth(2);
        spike_depth_group = discretize(spike_depths,depth_group_edges);
        n_depths = max(spike_depth_group);

        % Get MUA binned by widefield framerate
        % (NOT DEFINED IN THIS, SO DEFINING HERE)
%         sample_rate = (1/mean(diff(wf_t)));
        sample_rate = 35;
        skip_seconds = 60;
        time_bins = skip_seconds:1/sample_rate:timelite.timestamps(end)-skip_seconds;
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

        binned_spikes = zeros(max(spike_depth_group),length(time_bins)-1);
        for curr_depth = 1:max(spike_depth_group)
            curr_spike_times = spike_times_timelite(spike_depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end

        % Baseline binned spikes by pre-stim periods
        baseline_periods = stimOn_times + [-0.5,0];
        spike_baselines = nanmean(cell2mat(cellfun(@(x) nanmean(binned_spikes(:,x),2), ...
            arrayfun(@(x) time_bin_centers > baseline_periods(x,1) & ...
            time_bin_centers < baseline_periods(x,2),1:length(stimOn_times), ...
            'uni', false),'uni',false)),2);
        binned_spikes_baselined = binned_spikes - spike_baselines;

        % Task regressors
        stim_regressors = histcounts(stimOn_times,time_bins);
        reward_regressors = histcounts(reward_times,time_bins);

        stim_move_regressors = histcounts(stim_move_time,time_bins);
        nonstim_move_times = ...
            setdiff(timelite.timestamps(find(diff(wheel_move) == 1)+1), ...
            stim_move_regressors);
        nonstim_move_regressors = histcounts(nonstim_move_times,time_bins);

        % Concatenate selected regressors, set parameters
        task_regressors = {stim_regressors;reward_regressors;stim_move_regressors;nonstim_move_regressors};
        task_regressor_labels = {'Stim','Reward','Stim move','Nonstim move'};

        task_t_shifts = { ...
            [-0.2,1]; ... % stim
            [-0.2,1];  ... % outcome
            [-0.2,1];  ... % nonstim move
            [-0.2,1]};    % stim move

        task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),task_t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;

        [mua_task_k,mua_task_long,mua_task_expl_var,mua_task_reduced_long] = ...
            AP_regresskernel(task_regressors,binned_spikes_baselined, ...
            task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);

        mua_task_k_permute = cellfun(@(x) permute(x,[3,2,1]),mua_task_k,'uni',false);

        % Store MUA task kernel (raw - to combine and normalize later)
        day_mua_task_k{curr_recording} = mua_task_k_permute;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    mua_task_k_all{curr_animal} = day_mua_task_k;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Save
save_fn = fullfile(data_path,'mua_task_k');
save(save_fn,'mua_task_k_all')
fprintf('Saved %s\n',save_fn);


%% Visually responsive units

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

plot_depths = false;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set MUA length (microns)
mua_length = 200;

stim_responsive_cells = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    stim_responsive_cells_day = cell(length(recordings),1);

    if plot_depths
        figure;
        h = tiledlayout(1,length(recordings));
        title(h,animal);
    end

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Find striatum boundaries (skip if no marked striatum)
        AP_longstriatum_find_striatum_depth
        if any(isnan(striatum_depth))
            continue
        end

        % Discretize spikes by depth
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(2);       
        depth_group_edges(end) = striatum_depth(2);
        template_depth_group = discretize(template_depths,depth_group_edges);
        n_depths = max(template_depth_group);

        % Get responsive units
        % (get quiescent trials)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');

        % (for vis passive: right-side stim)
        stim_type = vertcat(trial_events.values.TrialStimX);
        use_align = stimOn_times(stim_type(1:length(stimOn_times)) == 90 & quiescent_trials);

        baseline_t = [-0.2,0];
        response_t = [0,0.2];

        baseline_bins = use_align + baseline_t;
        response_bins = use_align + response_t;

        event_bins = [baseline_bins,response_bins];
        spikes_binned_continuous = histcounts2(spike_times_timelite,spike_templates, ...
            reshape([baseline_bins,response_bins]',[],1),1:size(templates,1)+1);

        event_spikes = permute(reshape(spikes_binned_continuous(1:2:end,:),2, ...
            size(event_bins,1),[]),[2,1,3]);

        event_response = squeeze(mean(diff(event_spikes,[],2),1));

        n_shuff = 10000;
        event_response_shuff = cell2mat(arrayfun(@(shuff) ...
            squeeze(mean(diff(ap.shake(event_spikes,2),[],2),1)), ...
            1:n_shuff,'uni',false));

        event_response_rank = tiedrank(horzcat(event_response,event_response_shuff)')';
        event_response_p = event_response_rank(:,1)./(n_shuff+1);

        responsive_cells = event_response_p > 0.95;

        depth_responsive_n = accumarray(template_depth_group(~isnan(template_depth_group)), ...
            responsive_cells(~isnan(template_depth_group)),[n_depths,1],@sum);

        depth_total_n = accumarray(template_depth_group(~isnan(template_depth_group)), ...
            responsive_cells(~isnan(template_depth_group)),[n_depths,1],@length);

        % Store responsive cells
        stim_responsive_cells_day{curr_recording} = [depth_responsive_n,depth_total_n];

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    stim_responsive_cells{curr_animal} = stim_responsive_cells_day;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Save
save_fn = fullfile(data_path,'stim_cells');
save(save_fn,'stim_responsive_cells')
fprintf('Saved %s\n',save_fn);


%% TAN PSTHs (passive)

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

% Set MUA length (microns)
mua_length = 200;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

tan_psth_all = cell(length(animals),1);
tan_info_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    disp(animal);

    use_workflow = 'lcr_passive';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    tan_psth_day = cell(length(recordings),1);
    tan_info_day = cell(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Find striatum boundaries (skip if no marked striatum)
        AP_longstriatum_find_striatum_depth
        if any(isnan(striatum_depth))
            continue
        end

        % Discretize spikes by depth
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(2);
        depth_group_edges(end) = striatum_depth(2);
        n_depths = length(depth_group_edges)-1;

        % Initialize TAN PSTH by depth
        tan_info_day{curr_recording} = cell(n_depths,3);
        tan_psth_day{curr_recording} = cell(n_depths,1);

        % Classify tans
        AP_longstriatum_classify_striatal_units;        

        % Use TANs in "striatum" and depth-sort
        use_tans = find(striatum_celltypes.tan);

        % Align times (quiescent stim)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        stim_x = vertcat(trial_events.values.TrialStimX);
        % sometimes not all trials shown?
        stim_x = stim_x(1:length(stimOn_times));
        align_times = cellfun(@(x) stimOn_times(stim_x == x & quiescent_trials), ...
            num2cell(unique(stim_x)),'uni',false);

        % If TANs, store PSTHs
        if any(use_tans)

            % Get TAN PSTH
            [use_spikes,use_spike_groups] = ismember(spike_templates,use_tans);
            tan_psth = ...
                ap.psth(spike_times_timelite(use_spikes), ...
                align_times,use_spike_groups(use_spikes), ...
                'smoothing',50);

            % Store TAN PSTH and info by depth group
            tan_depth_group = discretize(template_depths(use_tans),depth_group_edges);
            for curr_depth = unique(tan_depth_group)'
                % PSTH
                tan_psth_day{curr_recording}{curr_depth,1} = ...
                    tan_psth(tan_depth_group == curr_depth,:,:);

                % Info
                curr_depth_tans = use_tans(tan_depth_group == curr_depth);
                tan_info_day{curr_recording}{curr_depth,1} = repmat(animal,size(curr_depth_tans));
                tan_info_day{curr_recording}{curr_depth,2} = repmat(rec_day,size(curr_depth_tans));
                tan_info_day{curr_recording}{curr_depth,3} = curr_depth_tans;
            end            

        end

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});
        ap.print_progress_fraction(curr_recording,length(recordings));
    end

    % Store animal data
    tan_psth_all{curr_animal} = tan_psth_day;
    tan_info_all{curr_animal} = tan_info_day;

end

% Save
save_fn = fullfile(data_path,'tan_psth_all_passive');
save(save_fn,'tan_psth_all','tan_info_all')
fprintf('Saved %s\n',save_fn);


%% TAN PSTHs (task)

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

% Set MUA length (microns)
mua_length = 200;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

tan_psth_all = cell(length(animals),1);
tan_info_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    disp(animal);

    use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    tan_psth_day = cell(length(recordings),1);
    tan_info_day = cell(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Find striatum boundaries (skip if no marked striatum)
        AP_longstriatum_find_striatum_depth
        if any(isnan(striatum_depth))
            continue
        end

        % Discretize spikes by depth
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(2);
        depth_group_edges(end) = striatum_depth(2);
        n_depths = length(depth_group_edges)-1;

        % Initialize TAN PSTH by depth
        tan_info_day{curr_recording} = cell(n_depths,3);
        tan_psth_day{curr_recording} = cell(n_depths,1);

        % Classify tans
        AP_longstriatum_classify_striatal_units;        

        % Use TANs in "striatum" and depth-sort
        use_tans = find(striatum_celltypes.tan);

        % Align times (stim, move, reward)
        align_times = {stimOn_times(1:n_trials),stim_move_time(1:n_trials),reward_times_task};

        % If TANs, store PSTHs
        if any(use_tans)

            % Get TAN PSTH
            [use_spikes,use_spike_groups] = ismember(spike_templates,use_tans);
            tan_psth = ...
                ap.psth(spike_times_timelite(use_spikes), ...
                align_times,use_spike_groups(use_spikes), ...
                'smoothing',50);

            % Store TAN PSTH and info by depth group
            tan_depth_group = discretize(template_depths(use_tans),depth_group_edges);
            for curr_depth = unique(tan_depth_group)'
                % PSTH
                tan_psth_day{curr_recording}{curr_depth,1} = ...
                    tan_psth(tan_depth_group == curr_depth,:,:);

                % Info
                curr_depth_tans = use_tans(tan_depth_group == curr_depth);
                tan_info_day{curr_recording}{curr_depth,1} = repmat(animal,size(curr_depth_tans));
                tan_info_day{curr_recording}{curr_depth,2} = repmat(rec_day,size(curr_depth_tans));
                tan_info_day{curr_recording}{curr_depth,3} = curr_depth_tans;
            end            

        end

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});
        ap.print_progress_fraction(curr_recording,length(recordings));
    end

    % Store animal data
    tan_psth_all{curr_animal} = tan_psth_day;
    tan_info_all{curr_animal} = tan_info_day;

end

% Save
save_fn = fullfile(data_path,'tan_psth_all_task');
save(save_fn,'tan_psth_all','tan_info_all')
fprintf('Saved %s\n',save_fn);

%% Cortical maps

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

% Set MUA length (microns)
mua_length = 200;

ctx_map_all = cell(length(animals),1);
ctx_expl_var_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

%     use_workflow = 'lcr_passive';
    use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    ctx_map = cell(length(recordings),1);
    ctx_expl_var = cell(length(recordings),1);
    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts.ephys = true;
        load_parts.widefield = true;
        load_parts.widefield_align = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % If not aligned, skip
        if ~load_parts.widefield_align
            continue
        end

        % Find striatum boundaries (skip if no marked striatum)
        AP_longstriatum_find_striatum_depth
        if any(isnan(striatum_depth))
            continue
        end

        % Discretize spikes by depth
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(2);
        depth_group_edges(end) = striatum_depth(2);
        depth_group = discretize(spike_depths,depth_group_edges);

        % Get MUA binned by widefield
        sample_rate = (1/mean(diff(wf_t)));
        skip_seconds = 60;
        time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

        binned_spikes = zeros(max(depth_group),length(time_bins)-1);
        for curr_depth = 1:max(depth_group)
            curr_spike_times = spike_times_timelite(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end

        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;

        use_svs = 1:200;
        kernel_t = [0,0];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        lambda = 20;
        zs = [false,false];
        cvfold = 5;
        return_constant = false;
        use_constant = true;

        fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

        % Regress cortex to spikes
        [k,predicted_spikes,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample, ...
            binned_spikes_std,kernel_frames, ...
            lambda,zs,cvfold,return_constant,use_constant);

        ctx_map{curr_recording} = permute(k,[1,3,2]);
        ctx_expl_var{curr_recording} = explained_var.total;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

        ap.print_progress_fraction(curr_recording,length(recordings));

    end

    ctx_map_all{curr_animal} = ctx_map;
    ctx_expl_var_all{curr_animal} = ctx_expl_var;

    fprintf('Done %s\n',animal);

end

% Save
save_fn = fullfile(data_path,'ctx_maps_task');
save(save_fn,'ctx_map_all','ctx_expl_var_all');
fprintf('Saved %s\n',save_fn);


%% Cortical regression

clearvars

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set MUA length (microns)
mua_length = 200;

day_mua_all = cell(length(animals),1);
day_mua_predicted_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    disp(animal);

    use_workflow = 'lcr_passive';
    %     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_mua = cell(length(recordings),1);
    day_mua_predicted = cell(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts.ephys = true;
        load_parts.widefield = true;
        load_parts.widefield_align = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % If not aligned, skip
        if ~load_parts.widefield_align
            continue
        end

        % Find striatum boundaries (skip if no marked striatum)
        AP_longstriatum_find_striatum_depth
        if any(isnan(striatum_depth))
            continue
        end

        % Discretize spikes by depth
        depth_group_edges = striatum_depth(1):mua_length:striatum_depth(2);
        spike_depth_group = discretize(spike_depths,depth_group_edges);
        n_depths = max(spike_depth_group);

        % Get MUA binned by widefield
        sample_rate = (1/mean(diff(wf_t)));
        skip_seconds = 60;
        time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

        binned_spikes = zeros(n_depths,length(time_bins)-1);
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timelite(spike_depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end

        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;

        use_svs = 1:100;
        kernel_t = [-0.1,0.1];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        lambda = 5;
        zs = [false,false];
        cvfold = 5;
        return_constant = true;
        use_constant = true;

        fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

        % Regress cortex to spikes
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample, ...
            binned_spikes_std,kernel_frames, ...
            lambda,zs,cvfold,return_constant,use_constant);

        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});

        % Stim times (quiescent trials only)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        if isfield(trial_events.values,'TrialStimX')
            stim_x = vertcat(trial_events.values.TrialStimX);
            % sometimes not all trials shown?
            stim_x = stim_x(1:length(stimOn_times));
            align_times = cellfun(@(x) stimOn_times(stim_x == x & quiescent_trials), ...
                num2cell(unique(stim_x)),'uni',false);
        else
            % Stim times (task)
            % (skip first trial - somes very short iti?)
            % (also skip stim_to_move negative - bad quiescence))
            use_trials = true(n_trials,1);
            use_trials(1) = false;
            use_trials(stim_to_move < 0) = false;

            align_times = {stimOn_times(use_trials),stim_move_time(use_trials),reward_times_task};
        end

        % Align measured and predicted activity to align times
        aligned_spikes_mean = nan(n_depths,length(t_centers),length(align_times));
        aligned_predicted_spikes_mean = nan(n_depths,length(t_centers),length(align_times));
        for curr_align = 1:length(align_times)
            peri_event_t = align_times{curr_align} + t_centers;

            aligned_spikes = reshape(interp1(time_bin_centers,binned_spikes', ...
                peri_event_t), ...
                size(peri_event_t,1),length(t_centers),[]);

            aligned_predicted_spikes = reshape(interp1(time_bin_centers, ...
                predicted_spikes',peri_event_t), ...
                size(peri_event_t,1),length(t_centers),[]);

            aligned_spikes_mean(:,:,curr_align) = ...
                permute(nanmean(aligned_spikes,1),[3,2,1]);

            aligned_predicted_spikes_mean(:,:,curr_align) = ...
                permute(nanmean(aligned_predicted_spikes,1),[3,2,1]);

        end

        % Store mean aligned ROI trace
        day_mua{curr_recording} = aligned_spikes_mean;
        day_mua_predicted{curr_recording} = aligned_predicted_spikes_mean;

        % Prep for next loop
        clearvars('-except',preload_vars{:});
        ap.print_progress_fraction(curr_recording,length(recordings));
    end

    day_mua_all{curr_animal} = day_mua;
    day_mua_predicted_all{curr_animal} = day_mua_predicted;

end

% Save
save_fn = fullfile(data_path,'ctx_str_prediction');
save(save_fn,'day_mua_all','day_mua_predicted_all');
fprintf('Saved %s\n',save_fn);

%% ~~~~~~~~~~~ BATCH (trial data)

%% Passive

% Initialize
clearvars
trial_data_all = struct;

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

% Loop through animals
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    % Loop through recordings
    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts.ephys = true;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % Grab trial data
        AP_longstriatum_grab_trials;

         % Store trial data in master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_recording,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end

        % Prep for next loop
        clearvars('-except',preload_vars{:});
        ap.print_progress_fraction(curr_recording,length(recordings));

    end
    fprintf('Done %s\n',animal);
end

% Save
save_fn = fullfile(data_path,'trial_data_passive');
save(save_fn,'trial_data_all');
fprintf('Saved %s\n',save_fn);


%% Task

% Initialize
clearvars
trial_data_all = struct;

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

% Loop through animals
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    % Loop through recordings
    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts.ephys = true;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % Grab trial data
        AP_longstriatum_grab_trials;

         % Store trial data in master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_recording,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end

        % Prep for next loop
        clearvars('-except',preload_vars{:});
        ap.print_progress_fraction(curr_recording,length(recordings));

    end
    fprintf('Done %s\n',animal);
end

% Save
save_fn = fullfile(data_path,'trial_data_task');
save(save_fn,'trial_data_all');
fprintf('Saved %s\n',save_fn);


%% |--> ~~~~~~~~~~~ BATCH ANALYSIS

%% Plot all recording locations

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

animal_col = [brewermap(8,'Dark2');brewermap(8,'Set2')];
ccf_draw = ap.ccf_draw;
ccf_draw.draw_name('Caudoputamen');

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    probe_color = animal_col(curr_animal,:);
    ccf_draw.draw_probes_nte(animal,probe_color);
end

%% Plot striatum boundaries

% Get animals
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(data_path,'animals'));

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    figure;
    h = tiledlayout(1,length(recordings));
    title(h,animal);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Find striatum boundaries
        AP_longstriatum_find_striatum_depth

        % Plot units and striatum boundary
        ha = nexttile(h);
        ap.plot_unit_depthrate(spike_templates,template_depths,probe_areas,ha);
        title(ha,rec_day);
        yline([striatum_depth(1),striatum_depth(2)],'r','linewidth',3);
        drawnow;
    end
end




%% MUA/responsive/TANs grouped by k-means

%%%%%% Load just normal-task mice

% Load data
am_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(am_data_path,'bhv.mat'));
load(fullfile(am_data_path,'mua_passive.mat'));
% load(fullfile(am_data_path,'mua_passive_notan.mat'));

% load(fullfile(am_data_path,'mua_task.mat'));
% load(fullfile(am_data_path,'mua_task_k.mat'));

% load(fullfile(am_data_path,'ctx_maps_passive.mat'));
load(fullfile(am_data_path,'ctx_maps_task.mat'));

load(fullfile(am_data_path,'stim_cells.mat'));
load(fullfile(am_data_path,'tan_psth_all_task.mat'));

% %%%%%% Load/concatenate normal-task and fixed quiescence mice
% 
% am_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
% a1 = load(fullfile(am_data_path,'bhv.mat'));
% b1 = load(fullfile(am_data_path,'mua_passive.mat'));
% c1 = load(fullfile(am_data_path,'ctx_maps_task.mat'));
% d1 = load(fullfile(am_data_path,'stim_cells.mat'));
% e1 = load(fullfile(am_data_path,'tan_psth_all.mat'));
% 
% am_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data\fixed_quiescence';
% a2 = load(fullfile(am_data_path,'bhv.mat'));
% b2 = load(fullfile(am_data_path,'mua_passive.mat'));
% c2 = load(fullfile(am_data_path,'ctx_maps_task.mat'));
% d2 = load(fullfile(am_data_path,'stim_cells.mat'));
% e2 = load(fullfile(am_data_path,'tan_psth_all.mat'));
% 
% bhv = struct('learned_day',[{a1.bhv.learned_day},{a2.bhv.learned_day}]);
% day_mua_all = vertcat(b1.day_mua_all,b2.day_mua_all);
% ctx_expl_var_all = vertcat(c1.ctx_expl_var_all,c2.ctx_expl_var_all);
% ctx_map_all = vertcat(c1.ctx_map_all,c2.ctx_map_all);
% stim_responsive_cells = vertcat(d1.stim_responsive_cells,d2.stim_responsive_cells);
% tan_psth_all = vertcat(e1.tan_psth_all,e2.tan_psth_all);
% 
% %%%%%%%


% Convert cortex maps V to px
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);
ctx_map_px = cellfun(@(x) cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),x,'uni',false),ctx_map_all,'uni',false);

% Concatenate data
ctx_map_cat = cell2mat(permute(vertcat(ctx_map_px{:}),[2,3,1]));
ctx_expl_var_cat = cell2mat(vertcat(ctx_expl_var_all{:}));

str_mua_cat = cell2mat(vertcat(day_mua_all{:}));

stim_responsive_cells_cat = cell2mat(vertcat(stim_responsive_cells{:}));

% Get grouping indicies for each depth
ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day}','uni',false);

animal_idx = cell2mat(cellfun(@(animal,data) repmat(animal,size(horzcat(data{:}),2),1), ...
    num2cell(1:length(ctx_map_all))',ctx_map_all,'uni',false));

day_idx = cell2mat(cellfun(@(ld,data) repmat(ld,size(data,2),1), ...
    num2cell(cell2mat(cellfun(@(x) (1:length(x))',ctx_map_all,'uni',false))), ...
    vertcat(ctx_map_all{:}),'uni',false));

ld_idx = cell2mat(cellfun(@(ld,data) repmat(ld,size(data,2),1), ...
    num2cell(vertcat(ld{:})),vertcat(ctx_map_all{:}),'uni',false));

depth_idx = cell2mat(cellfun(@(data) (1:size(data,2))', ...
    vertcat(ctx_map_all{:}),'uni',false));

% K-means cortex maps
% (only use maps that explain x% variance)
expl_var_cutoff = 0.00; 
figure;
subplot(2,1,1);
plot(ctx_expl_var_cat,'.k');
yline(expl_var_cutoff,'r','Use maps cutoff');
subplot(2,1,2);
histogram(ctx_expl_var_cat,linspace(-0.2,1,100));
xline(expl_var_cutoff,'r','Use maps cutoff');

use_maps = ctx_expl_var_cat > expl_var_cutoff;

n_k = 3;
kidx = nan(size(ctx_map_cat,3),1);
[kidx(use_maps),kmeans_map] = kmeans(reshape(ctx_map_cat(:,:,use_maps), ...
    [],sum(use_maps))',n_k,'distance','correlation','replicates',5);
kmeans_map = reshape(kmeans_map',[size(ctx_map_cat,[1,2]),n_k]);

figure;tiledlayout(1,n_k,'tilespacing','none','tileindexing','columnmajor');
for i = 1:n_k
    nexttile;
    imagesc(kmeans_map(:,:,i));
    axis image off;
    clim(max(abs(clim)).*[-1,1]);
    colormap(ap.colormap('PWG',[],1.5));
    ap.wf_draw('ccf','k');
end

%%%%%

% get rid of maps that aren't at least x correlation to template 

k_corr_thresh = 0.5;

figure;tiledlayout(n_k,1);
for curr_kidx = 1:n_k
    curr_maps = find(use_maps & kidx == curr_kidx);
    curr_k_corr = corr(reshape(kmeans_map(:,:,curr_kidx),[],1), ...
        reshape(ctx_map_cat(:,:,curr_maps), ...
        prod(size(ctx_map_cat,[1,2])),[]));

    kidx(curr_maps(curr_k_corr < k_corr_thresh)) = NaN;

    nexttile;
    plot(curr_k_corr,'.k');
    yline(k_corr_thresh,'r','Kernel correlation threshold');
end

%%%%%%

% Plot MUA by k-means cluster
use_align = 3;
use_t = [500,700];
plot_ld = -3:3;
figure;h = tiledlayout(2,n_k,'tileindexing','columnmajor');
for curr_kidx = 1:n_k

    curr_use_ctx = kidx == curr_kidx;
    [grp,~,grp_idx] = unique([animal_idx(curr_use_ctx),ld_idx(curr_use_ctx)],'rows');
    
    use_grp = ismember(grp(:,2),plot_ld);
    ld_unique = unique(grp(use_grp,2));

%     ctx_map_depthgrp = ap.groupfun(@mean,ctx_map_cat(:,:,curr_use_ctx),[],[],grp_idx);
%     ctx_map_ldgrp = ap.groupfun(@mean,ctx_map_depthgrp,[],[],grp(:,2));
%     figure; h2 = tiledlayout('flow','TileSpacing','none');
%     c = max(abs(ctx_map_ldgrp),[],'all').*[-1,1].*0.8;
%     for i = ld_unique'
%         nexttile(h2);
%         imagesc(ctx_map_ldgrp(:,:,ld_unique==i));
%         axis image off;
%         clim(c);
%         ap.wf_draw('ccf','k');
%         title(i);
%     end
%     colormap(ap.colormap('PWG',[],1.5));

    str_mua_depthgrp = ap.groupfun(@sum,str_mua_cat(curr_use_ctx,:,use_align),grp_idx,[]);

    softnorm = 10;
    str_mua_depthgrp_norm = (str_mua_depthgrp-nanmean(str_mua_depthgrp(:,1:500),2))./ ...
        (nanmean(str_mua_depthgrp(:,1:500),2) + softnorm);
    str_mua_depthgrp_norm_tmax = max(str_mua_depthgrp_norm(:,use_t(1):use_t(2)),[],2);
    str_mua_ldgrp_norm = ap.groupfun(@mean,str_mua_depthgrp_norm(use_grp,:),grp(use_grp,2),[]);

    nexttile(h); hold on; colororder(gca,ap.colormap('BKR',length(ld_unique)));
    plot(str_mua_ldgrp_norm','linewidth',2);

    nexttile(h); hold on; colororder(gca,[0.5,0.5,0.5]);
    arrayfun(@(x) plot(grp(use_grp & grp(:,1)==x,2), ...
        str_mua_depthgrp_norm_tmax(use_grp & grp(:,1)==x)),unique(grp(use_grp,1)));
    xlabel('Learned day');
    ylabel('MUA')

    [m,e] = grpstats(str_mua_depthgrp_norm_tmax(use_grp),grp(use_grp,2),{'mean','sem'});
    errorbar(ld_unique,m,e,'k','linewidth',2);
    drawnow;

    % Stats (shuffle stat LD-1 vs LD<-1)
    prelearn_days = grp(:,2) <= -1;
    stat_use_act = str_mua_depthgrp_norm_tmax(prelearn_days);
    stat_use_grp = grp(prelearn_days,:);

    stat_diff = mean(stat_use_act(stat_use_grp(:,2) == -1)) - mean(stat_use_act(stat_use_grp(:,2) < -1));

    n_shuff = 1000;
    stat_diff_shuff = nan(n_shuff,1);
    for i = 1:n_shuff
        curr_grp_shuff = AP_shake(stat_use_grp(:,2),[],stat_use_grp(:,1));
        stat_diff_shuff(i) = mean(stat_use_act(curr_grp_shuff == -1)) - mean(stat_use_act(curr_grp_shuff < -1));
    end
    stat_rank = tiedrank(vertcat(stat_diff,stat_diff_shuff));
    stat_p = 1-(stat_rank(1)./(n_shuff+2));

    fprintf('MUA stats - Kidx %d: p = %.3f\n',curr_kidx,stat_p);

end
linkaxes(h.Children(1:2:end));
linkaxes(h.Children(2:2:end));


% Plot responsive cell fraction by k-means cluster
figure;h = tiledlayout(1,n_k);
for curr_kidx = 1:n_k

    curr_use_ctx = kidx == curr_kidx;
    [grp,~,grp_idx] = unique([animal_idx(curr_use_ctx),ld_idx(curr_use_ctx)],'rows');

    stim_cells_combine = ap.groupfun(@sum,stim_responsive_cells_cat(curr_use_ctx,:),grp_idx,[]);
    stim_cells_frac  = stim_cells_combine(:,1)./stim_cells_combine(:,2);

    nexttile(h); hold on;
    arrayfun(@(x) plot(grp(grp(:,1)==x,2), stim_cells_frac(grp(:,1)==x), ...
        'color',[0.5,0.5,0.5]),unique(grp(:,1)));
    xlabel('Learned day');
    ylabel('Frac stim-responsive')

    [m,e] = grpstats(stim_cells_frac,grp(:,2),{'mean','sem'});
    errorbar(unique(grp(:,2)),m,e,'k','linewidth',2);
    drawnow;

    % Stats (shuffle stat LD-1 vs LD<-1)
    prelearn_days = grp(:,2) <= -1;
    stat_use_act = stim_cells_frac(prelearn_days);
    stat_use_grp = grp(prelearn_days,:);

    stat_diff = mean(stat_use_act(stat_use_grp(:,2) == -1)) - mean(stat_use_act(stat_use_grp(:,2) < -1));

    n_shuff = 1000;
    stat_diff_shuff = nan(n_shuff,1);
    for i = 1:n_shuff
        curr_grp_shuff = AP_shake(stat_use_grp(:,2),[],stat_use_grp(:,1));
        stat_diff_shuff(i) = mean(stat_use_act(curr_grp_shuff == -1)) - mean(stat_use_act(curr_grp_shuff < -1));
    end
    stat_rank = tiedrank(vertcat(stat_diff,stat_diff_shuff));
    stat_p = 1-(stat_rank(1)./(n_shuff+2));

    fprintf('MUA stats - Kidx %d: p = %.3f\n',curr_kidx,stat_p);

end


% Plot task kernel by k-means cluster
plot_kernel = 1;
mua_task_k_flat = cellfun(@transpose,vertcat(mua_task_k_all{:}),'uni',false);
mua_task_k = vertcat(mua_task_k_flat{:});
mua_task_k_cat = cell2mat(mua_task_k(:,plot_kernel));

use_t = [8,14];

figure;h = tiledlayout(2,n_k,'tileindexing','columnmajor');
for curr_kidx = 1:n_k

    curr_use_ctx = kidx == curr_kidx;
    [grp,~,grp_idx] = unique([animal_idx(curr_use_ctx),ld_idx(curr_use_ctx)],'rows');
    ld_unique = unique(grp(:,2));

    str_mua_depthgrp = ap.groupfun(@sum,str_mua_cat(curr_use_ctx,:,1),grp_idx,[]);
    str_mua_depthgrp_baseline = nanmean(str_mua_depthgrp(:,1:500),2);
    
    str_k_depthgrp = ap.groupfun(@sum,mua_task_k_cat(curr_use_ctx,:),grp_idx,[])/(1/35); % convert to rate - currently in spikes

    str_mua_depthgrp_norm = str_k_depthgrp./str_mua_depthgrp_baseline;
    str_mua_depthgrp_norm_tmax = max(str_mua_depthgrp_norm(:,use_t(1):use_t(2)),[],2);

    str_mua_ldgrp_norm = ap.groupfun(@mean,str_mua_depthgrp_norm,grp(:,2),[]);
    nexttile(h); hold on; colororder(gca,ap.colormap('BKR',7));
    plot(str_mua_ldgrp_norm(ismember(ld_unique,-3:3),:)','linewidth',2);

    nexttile(h); hold on; colororder(gca,lines);
    arrayfun(@(x) plot(grp(grp(:,1)==x,2), ...
        str_mua_depthgrp_norm_tmax(grp(:,1)==x)),unique(grp(:,1)));
    xlabel('Learned day');
    ylabel('MUA')

    [m,e] = grpstats(str_mua_depthgrp_norm_tmax,grp(:,2),{'mean','sem'});
    errorbar(ld_unique,m,e,'k');
    drawnow;

end
linkaxes(h.Children(2:2:end));


% Plot behavior by LD
ld_unique = unique(vertcat(ld{:}));

rxn_idx_cat = cat(1,bhv.rxn_idx);
rxn_med_cat = cat(1,bhv.stim_to_move_med);
stim_move_ratio_cat = cat(1,bhv.stim_move_frac_ratio);
frac_move_stimalign_cat = cat(1,bhv.frac_move_stimalign);

rxn_med_ldgrp = ap.groupfun(@nanmean,rxn_med_cat,vertcat(ld{:}),[]);
stim_move_ratio_ldgrp = ap.groupfun(@nanmean,stim_move_ratio_cat,vertcat(ld{:}),[]);
frac_move_stimalign_ldgrp = ap.groupfun(@nanmean,frac_move_stimalign_cat,vertcat(ld{:}),[]);

figure; tiledlayout(1,2);

nexttile;
yyaxis left;
plot(ld_unique,rxn_med_ldgrp,'-o','linewidth',2);
ylabel('Reaction time');
yyaxis right;
plot(ld_unique,stim_move_ratio_ldgrp,'-o','linewidth',2);
ylabel('Stim move ratio');
xline(0,'r');

nexttile;
plot(frac_move_stimalign_ldgrp(ismember(ld_unique,-3:3),:)','linewidth',2);
colororder(gca,ap.colormap('BKR',7));



% Behavior vs striatal activity
curr_kidx = 3;
curr_use_ctx = kidx == curr_kidx;
[grp,~,grp_idx] = unique([animal_idx(curr_use_ctx),ld_idx(curr_use_ctx)],'rows');
str_mua_depthgrp = ap.groupfun(@sum,str_mua_cat(curr_use_ctx,:,use_align),grp_idx,[]);
str_mua_depthgrp_norm = (str_mua_depthgrp-nanmean(str_mua_depthgrp(:,1:500),2))./ ...
    nanmean(str_mua_depthgrp(:,1:500),2);
str_mua_depthgrp_norm_tmax = max(str_mua_depthgrp_norm(:,use_t(1):use_t(2)),[],2);

[ld_unique,~,ld_grp] = unique(grp(:,2));
size_grid = [max(grp(:,1)),max(ld_grp)];
idx = sub2ind(size_grid,grp(:,1),ld_grp);
act_grid = nan(size_grid);
act_grid(idx) = str_mua_depthgrp_norm_tmax;


animal_idx = cell2mat(cellfun(@(x,ld) ...
    repmat(x,length(ld),1),num2cell(1:length(ld))',ld,'uni',false));
bhv_grp = [animal_idx,vertcat(ld{:})];

[ld_unique,~,ld_grp] = unique(bhv_grp(:,2));
size_grid = [max(bhv_grp(:,1)),max(ld_grp)];
idx = sub2ind(size_grid,bhv_grp(:,1),ld_grp);
bhv_grid = nan(size_grid);
bhv_grid(idx) = rxn_med_cat;

% (this is super messy)
bhv_grid_use = bhv_grid(:,ismember(ld_unique,unique(grp(:,2))));

act_bhv_corr = corr(act_grid,bhv_grid_use,'rows','pairwise');
figure;imagesc(unique(grp(:,2)),unique(grp(:,2)),act_bhv_corr);
set(gca,'YDir','normal');
clim([-1,1]);
colormap(ap.colormap('BWR',[],1));
axis image;
xlim([-3.5,2.5]);ylim(xlim);
xlabel('LD: Behavior');
ylabel('LD: Striatal activity');


ld_x = unique(grp(:,2));
figure; hold on
plot(reshape(bhv_grid_use,[],1),reshape(act_grid,[],1),'.k');
plot(reshape(bhv_grid_use(:,ld_x<-1),[],1),reshape(act_grid(:,ld_x<-1),[],1),'or');
plot(reshape(bhv_grid_use(:,ismember(ld_x,-1)),[],1),reshape(act_grid(:,ismember(ld_x,-1)),[],1),'om');
plot(reshape(bhv_grid_use(:,ld_x>=0),[],1),reshape(act_grid(:,ld_x>=0),[],1),'ob');

figure; hold on;
ld_x_discrete = discretize(ld_x,[-Inf,-2,-1,0,Inf],'IncludedEdge','right');
ld_x_discrete_col = [0,0,0;0,0,1;1,0,0;0,0.7,0];

scatter(reshape(bhv_grid_use,[],1),reshape(act_grid,[],1),10, ...
    ld_x_discrete_col(reshape(repmat(ld_x_discrete',size(act_grid,1),1),[],1),:));
gscatter(nanmean(bhv_grid_use,1),nanmean(act_grid,1),ld_x_discrete,ld_x_discrete_col,'o',10,'filled')
errorbar(nanmean(bhv_grid_use,1),nanmean(act_grid,1), ...
    AP_sem(act_grid,1),AP_sem(act_grid,1), ...
    AP_sem(bhv_grid_use,1),AP_sem(bhv_grid_use,1),'k');
legend({'','LD<-1','LD=-1','LD==0','LD>1'},'location','best')
xlabel('Behavior');
ylabel('Striatal MUA');

figure; hold on
ld_x_diff = ld_x(2:end);
a = diff(bhv_grid_use,[],2);
b = diff(act_grid,[],2);
plot(reshape(a,[],1),reshape(b,[],1),'.k');
plot(reshape(a(:,ismember(ld_x_diff,-3:-2)),[],1),reshape(b(:,ismember(ld_x_diff,-3:-2)),[],1),'or');
plot(reshape(a(:,ismember(ld_x_diff,-1)),[],1),reshape(b(:,ismember(ld_x_diff,-1)),[],1),'om');
plot(reshape(a(:,ismember(ld_x_diff,0:2)),[],1),reshape(b(:,ismember(ld_x_diff,0:2)),[],1),'ob');
xlabel('\DeltaBehavior');
ylabel('\DeltaActivity');




% (testing: make template map from visually-responsive regions)

% % (OLD: weighted sum of maps by stim cell fraction)
% x = stim_responsive_cells_cat(:,1)./stim_responsive_cells_cat(:,2);
% x(isnan(x)) = 0;
% r = ctx_map_cat.*permute(x,[2,3,1]);
% r2 = sqrt(sum(r.^2,3));
% 
% figure;imagesc(r2); 
% axis image off; colormap(ap.colormap('WK',[],2));
% clim(max(abs(clim)).*[-1,1]);
% ap.wf_draw('ccf','k');
% colormap(ap.colormap('PWG',[],1.5));
% 
% a = corr(r2(:),reshape(ctx_map_cat.^2,[],size(ctx_map_cat,3)));
% a_thresh = 0.45;
% figure;plot(a,'.k')
% yline(a_thresh,'r','thresh');
% 
% kidx = nan(size(ctx_map_cat,3),1);
% kidx(a>a_thresh) = 1;
% n_k = 1;
% 
% figure;imagesc(nanmean(ctx_map_cat(:,:,kidx==1),3));axis image
% clim(max(abs(clim)).*[-1,1]);
% ap.wf_draw('ccf','k');
% colormap(ap.colormap('PWG',[],1.5));

% NEW: threshold for fraction of stim cells, nonvis maps are maps from same
% day under threshold, vis maps are more correlated to vis than nonvis
% (above included more data and looked better after analysis for some
% reason)

vis_rec_idx = x >= 0.1;

use_animaldays = unique([animal_idx(vis_maps),day_idx(vis_rec_idx)],'rows');
nonvis_rec_idx = ismember([animal_idx,day_idx],use_animaldays,'rows') & ~vis_rec_idx;

vis_map = nanmean(ctx_map_cat(:,:,vis_rec_idx),3);
nonvis_map = nanmean(ctx_map_cat(:,:,nonvis_rec_idx),3);

vis_map_corr = corr(vis_map(:),reshape(ctx_map_cat,[],size(ctx_map_cat,3)));
nonvis_map_corr = corr(nonvis_map(:),reshape(ctx_map_cat,[],size(ctx_map_cat,3)));

vis_map_idx = vis_map_corr - nonvis_map_corr > 0;

% (to use kidx code)
n_k = 2;
kidx = nan(size(ctx_map_cat,3),1);
kidx(vis_map_idx) = 1;
kidx(~vis_map_idx) = 2;

figure;tiledlayout(2,2);

nexttile;
imagesc(vis_map);
axis image off
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));
title('Vis template');

nexttile;
imagesc(nonvis_map);
axis image off
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));
title('Non-vis template');

nexttile;
imagesc(nanmean(ctx_map_cat(:,:,vis_map_idx),3));
axis image off
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));
title('Vis map mean');

nexttile;
imagesc(nanmean(ctx_map_cat(:,:,~vis_map_idx),3));
axis image off
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));
title('Non-vis map mean');



%% TAN analysis

% Load data
am_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(am_data_path,'bhv.mat'));
load(fullfile(am_data_path,'ctx_maps_task.mat'));

load(fullfile(am_data_path,'stim_cells.mat'));

% Convert cortex maps V to px
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);
ctx_map_px = cellfun(@(x) cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),x,'uni',false),ctx_map_all,'uni',false);

% Concatenate data
ctx_map_cat = cell2mat(permute(vertcat(ctx_map_px{:}),[2,3,1]));
ctx_expl_var_cat = cell2mat(vertcat(ctx_expl_var_all{:}));
stim_responsive_cells_cat = cell2mat(vertcat(stim_responsive_cells{:}));

% Get indicies for animal / learned day / depth
ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day}','uni',false);

animal_ld_idx_cell = cellfun(@(animal,ld,data) cellfun(@(ld,data) ...
    [repmat(animal,size(data,2),1),repmat(ld,size(data,2),1),(1:size(data,2))'], ...
    num2cell(ld),data,'uni',false), ...
    num2cell(1:length(ctx_map_all))',ld,ctx_map_all,'uni',false);
animal_ld_idx = cell2mat(vertcat(animal_ld_idx_cell{:}));

% Get depths corresponding to visually responsive cells
x = stim_responsive_cells_cat(:,1)./stim_responsive_cells_cat(:,2);
x(isnan(x)) = 0;
r = ctx_map_cat.*permute(x,[2,3,1]);
r2 = sqrt(sum(r.^2,3));

figure;imagesc(r2); 
axis image off; colormap(ap.colormap('WK',[],2));
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));

a = corr(r2(:),reshape(ctx_map_cat.^2,[],size(ctx_map_cat,3)));
a_thresh = 0.45;
figure;plot(a,'.k')
yline(a_thresh,'r','thresh');

kidx = nan(size(ctx_map_cat,3),1);
kidx(a>a_thresh) = 1;
n_k = 1;

figure;imagesc(nanmean(ctx_map_cat(:,:,kidx==1),3));axis image
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));


% Load TANs from task and passive
load(fullfile(am_data_path,'tan_psth_all_passive.mat'));

tan_psth_passive_cat = vertcat(tan_psth_all{:});
tan_psth_passive_cat = vertcat(tan_psth_passive_cat{:});

load(fullfile(am_data_path,'tan_psth_all_task.mat'));

tan_psth_task_cat = vertcat(tan_psth_all{:});
tan_psth_task_cat = vertcat(tan_psth_task_cat{:});

% Combine TANs by LD section
% ld_split = [animal_ld_idx(:,2) < -1, ...
%     animal_ld_idx(:,2) == -1, ... 
%     animal_ld_idx(:,2) >= 0];
ld_split = animal_ld_idx(:,2) == (-3:2);

figure; colormap(ap.colormap('BWR',[],1.5));
task_tl = tiledlayout(3,size(ld_split,2),'TileIndexing','ColumnMajor');

figure; colormap(ap.colormap('BWR',[],1.5));
passive_tl = tiledlayout(3,size(ld_split,2),'TileIndexing','ColumnMajor');

for curr_ld = 1:size(ld_split,2)
    curr_task_psth = vertcat(tan_psth_task_cat{kidx == 1 & ld_split(:,curr_ld)});
    curr_task_psth = curr_task_psth - mean(curr_task_psth(:,1:500,1),[2,3]);

    curr_passive_psth = vertcat(tan_psth_passive_cat{kidx == 1 & ld_split(:,curr_ld)});
    curr_passive_psth = curr_passive_psth - mean(curr_passive_psth(:,1:500,1),[2,3]);

    [~,sort_idx] = sort(mean(curr_task_psth(:,550:700,1),2));

    for curr_align = 1:size(curr_task_psth,3)
        nexttile(task_tl);
        imagesc(curr_task_psth(sort_idx,:,curr_align));
        clim([-10,10]);

        nexttile(passive_tl);
        imagesc(curr_passive_psth(sort_idx,:,curr_align));
        clim([-10,10]);
    end

end
title(task_tl,'Task');
title(passive_tl,'Passive');



%% Widefield responses

% Load data
am_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(am_data_path,'bhv.mat'));
load(fullfile(am_data_path,'wf_passive.mat'));
% load(fullfile(am_data_path,'wf_task.mat'));

% Load U master
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);

% Concatenate data
V_cat = cat(4,day_V_all{:});

% Get indicies for animal / learned day
ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day}','uni',false);
ld_unique = unique(vertcat(ld{:}));

V_animal_day_idx = cell2mat(cellfun(@(x,animal,ld) ...
    [repmat(animal,size(x,4),1),ld], ...
    day_V_all,num2cell(1:length(day_V_all))',ld,'uni',false));

% Average widefield by learned day
V_ldavg = permute(ap.groupfun(@nanmean,V_cat,[],[],[],V_animal_day_idx(:,2)),[1,2,4,3]);
px_ldavg = plab.wf.svd2px(U_master,V_ldavg);

% Plot widefield by learned day
plot_ld = -2:2;
plot_align = 3;
ap.imscroll(px_ldavg(:,:,:,ismember(ld_unique,plot_ld),plot_align));
axis image;
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1));

% |--> (draw ROI from above and plot heatmap, lines, and max fluorescence)
figure;
h = tiledlayout(1,3);

nexttile;imagesc([],plot_ld,roi.trace);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG',[],1.5));

ha = nexttile;
plot(roi.trace','linewidth',2);
colororder(ha,ap.colormap('BKR',size(roi.trace,1)));

use_frames = 15:30;
roi_max = squeeze(max(ap.wf_roi(U_master,permute(V_cat(:,use_frames,3,:),[1,2,4,3]),[],[],roi.mask),[],2));

use_grps = ismember(V_animal_day_idx(:,2),plot_ld);
roi_max_avg = ap.groupfun(@nanmedian,roi_max,V_animal_day_idx(:,2),[]);
roi_max_std = ap.groupfun(@AP_sem,roi_max,V_animal_day_idx(:,2),[]);

nexttile; hold on;
arrayfun(@(x) plot(V_animal_day_idx(V_animal_day_idx(:,1)==x,2), ...
    roi_max(V_animal_day_idx(:,1)==x)),1:max(V_animal_day_idx(:,1)));

plot_grp_idx = ismember(ld_unique,plot_ld);
errorbar(ld_unique(plot_grp_idx),roi_max_avg(plot_grp_idx), ...
    roi_max_std(plot_grp_idx),'k','linewidth',2);
xlabel('Learned day');
ylabel('Max fluorescence');

title(h,'ROI');








%% Cortex prediction
% (unused at the moment - cortex predicts striatum, but that would be
% expected if the synapse was strengthened?)

am_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(am_data_path,'bhv.mat'));
load(fullfile(am_data_path,'ctx_maps_passive.mat'));
load(fullfile(am_data_path,'ctx_str_prediction.mat'));


% Load master U, convert V to px
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);
ctx_map_px = cellfun(@(x) cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),x,'uni',false),ctx_map_all,'uni',false);

use_animals = ~cellfun(@isempty,ctx_map_all);

ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day}','uni',false);

% (learned day)
animal_ld_idx_cell = cellfun(@(animal,ld,data) cellfun(@(ld,data) ...
    [repmat(animal,size(data,2),1),repmat(ld,size(data,2),1),(1:size(data,2))'], ...
    num2cell(ld),data,'uni',false), ...
    num2cell(find(use_animals)),ld(use_animals),ctx_map_all(use_animals),'uni',false);
animal_ld_idx = cell2mat(vertcat(animal_ld_idx_cell{:}));

% (cardinal day)
animal_day_idx_cell = cellfun(@(animal,ld,data) cellfun(@(ld,data) ...
    [repmat(animal,size(data,2),1),repmat(ld,size(data,2),1),(1:size(data,2))'], ...
    num2cell(1:length(ld))',data,'uni',false), ...
    num2cell(find(use_animals)),ld(use_animals),ctx_map_all(use_animals),'uni',false);
animal_day_idx = cell2mat(vertcat(animal_day_idx_cell{:}));


str_mua_cat = cell2mat(vertcat(day_mua_all{use_animals}));
softnorm = 1;
r_norm = (str_mua_cat-nanmean(str_mua_cat(:,200:400),[2,3]))./(nanmean(str_mua_cat(:,200:400),[2,3])+softnorm);
use_t = [550,750];
r2 = nanmean(r_norm(:,use_t(1):use_t(2)),2);

% Concatenate maps and explained variance
ctx_map_cat = cell2mat(permute(vertcat(ctx_map_px{:}),[2,3,1]));
ctx_expl_var_cat = cell2mat(vertcat(ctx_expl_var_all{:}));



% Concatenate V's
V_cat = cat(4,day_V_all{:});

V_animal_day_idx = cell2mat(cellfun(@(x,animal,ld) ...
    [repmat(animal,size(x,4),1),ld], ...
    day_V_all(use_animals),num2cell(find(use_animals)),ld(use_animals),'uni',false));

ap.imscroll(ctx_map_cat);
axis image off;
clim(max(abs(clim)).*[-1,1].*0.5);
colormap(ap.colormap('BWR',[],1.5));
ap.wf_draw('ccf','k');


% (draw ROI over mVis)
vis_map_weights = roi.trace';

% (draw ROI over mpfc)
mpfc_map_weights = roi.trace';




% Get vis cortex maps, group striatal responses, plot by LD
ctx_thresh = 0.005; %0.005
use_expl_var = ctx_expl_var_cat > 0.05;
% curr_use_ctx = vis_map_weights >= ctx_thresh & mpfc_map_weights < ctx_thresh & use_expl_var; % vis-only
% curr_use_ctx = vis_map_weights < ctx_thresh & mpfc_map_weights >= ctx_thresh & use_expl_var; % mpfc-only
% curr_use_ctx = vis_map_weights >= ctx_thresh & mpfc_map_weights >= ctx_thresh & use_expl_var; % vis+mpfc
% curr_use_ctx = vis_map_weights >= ctx_thresh & use_expl_var; % vis (all)
curr_use_ctx = mpfc_map_weights >= ctx_thresh & use_expl_var; % mpfc (all)

figure;plot(vis_map_weights,mpfc_map_weights,'.k')
xline(ctx_thresh,'r');yline(ctx_thresh,'r');
xlabel('vis weights');
ylabel('mpfc weights');

[grp,~,grp_idx] = unique(animal_ld_idx(curr_use_ctx,1:2),'rows');

use_align = 3;
x1 = ap.groupfun(@sum,str_mua_cat(curr_use_ctx,:,use_align),grp_idx,[]);
x1_norm = (x1-nanmean(x1(:,1:500),2))./nanmean(x1(:,1:500),2);

x11t = ap.groupfun(@nanmean,x1_norm,grp(:,2),[]);
x11 = max(x1_norm(:,use_t(1):use_t(2)),[],2);

figure;imagesc([],unique(grp(:,2)),x11t)
ylabel('LD');
title('MUA');

x2 = ap.groupfun(@nanmean,ctx_map_cat(:,:,curr_use_ctx),[],[],grp_idx);
x22 = ap.groupfun(@nanmean,x2,[],[],grp(:,2));
figure; tiledlayout('flow','TileSpacing','none');
c = max(abs(x22),[],'all').*[-1,1].*0.8;
ld_unique = unique(grp(:,2));
for i = -3:3
    nexttile;
    imagesc(x22(:,:,ld_unique==i));
    axis image off;
    clim(c);
    ap.wf_draw('ccf','k');
    title(i);
end
colormap(ap.colormap('BWR',[],1.5));

ap.imscroll(x22,ld_unique);
axis image off;
clim(c);
ap.wf_draw('ccf','k');
colormap(ap.colormap('BWR',[],1.5));


figure; hold on;
set(gca,'colororder',jet(length(unique(grp(:,1)))));
arrayfun(@(x) plot(grp(grp(:,1)==x,2), x11(grp(:,1)==x)),unique(grp(:,1)));
xlabel('Learned day');
ylabel('MUA')

[m,e] = grpstats(x11,grp(:,2),{'mean','sem'});
errorbar(unique(grp(:,2)),m,e,'k','linewidth',2);


% Average widefield from animals/days used above
use_align = 3;
% use_v = ismember(V_animal_day_idx,grp,'rows');
use_v = true(size(V_animal_day_idx,1),1);

V_avg = permute(ap.groupfun(@nanmean,V_cat(:,:,use_align,use_v),[],[],[],V_animal_day_idx(use_v,2)),[1,2,4,3]);
px_avg = plab.wf.svd2px(U_master,V_avg);
ld_unique = unique(V_animal_day_idx(use_v,2));
ap.imscroll(px_avg(:,:,:,ismember(ld_unique,-3:3)));
axis image;
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));

use_frames = 15:28;
px_tavg = squeeze(max(px_avg(:,:,use_frames,:),[],3));
ap.imscroll(px_tavg,ld_unique);
axis image;
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));


use_frames = 20:23;
px_tmax_cat = squeeze(mean(plab.wf.svd2px(U_master,V_cat(:,use_frames,use_align,:)),3));

ap.imscroll(px_tmax_cat);
axis image;
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));

% (draw_roi);
use_roi_trace = roi.trace(use_v);
use_V_animal_day_idx = V_animal_day_idx(use_v,:);

figure; hold on;
arrayfun(@(animal) ...
    plot(use_V_animal_day_idx(use_V_animal_day_idx(:,1) == animal,2), ...
    use_roi_trace(use_V_animal_day_idx(:,1) == animal),'color',[0.5,0.5,0.5]), ...
    1:max(use_V_animal_day_idx(:,1)));

m = ap.groupfun(@nanmean,use_roi_trace,[],use_V_animal_day_idx(:,2));
e = ap.groupfun(@AP_sem,use_roi_trace,[],use_V_animal_day_idx(:,2));
errorbar(ld_unique,m,e,'k','linewidth',2);


% K-means cortex maps
% (only use maps that explain x% variance)
use_maps = ctx_expl_var_cat > 0.05;

n_k = 4;
kidx = nan(size(ctx_map_cat,3),1);
[kidx(use_maps),kmeans_map] = kmeans(reshape(ctx_map_cat(:,:,use_maps),[],sum(use_maps))',n_k,'distance','correlation','replicates',5);
kmeans_map = reshape(kmeans_map',[size(ctx_map_cat,[1,2]),n_k]);

figure;tiledlayout(1,n_k,'tilespacing','none','tileindexing','columnmajor');
for i = 1:n_k
    nexttile;
    imagesc(kmeans_map(:,:,i));
    axis image off;
    clim(max(abs(clim)).*[-1,1]);
    colormap(ap.colormap('PWG',[],1.5));
    ap.wf_draw('ccf','k');
end

% Plot MUA by k-means cluster
use_align = 3;
figure;h = tiledlayout(2,n_k,'tileindexing','columnmajor');
for curr_kidx = 1:n_k

    curr_use_ctx = kidx == curr_kidx;
    [grp,~,grp_idx] = unique(animal_ld_idx(curr_use_ctx,1:2),'rows');

    x2 = ap.groupfun(@mean,ctx_map_cat(:,:,curr_use_ctx),[],[],grp_idx);
    x22 = ap.groupfun(@mean,x2,[],[],grp(:,2));
    figure; h2 = tiledlayout('flow','TileSpacing','none');
    c = max(abs(x22),[],'all').*[-1,1].*0.8;
    ld_unique = unique(grp(:,2));
    for i = ld_unique'
        nexttile(h2);
        imagesc(x22(:,:,ld_unique==i));
        axis image off;
        clim(c);
        ap.wf_draw('ccf','k');
        title(i);
    end
    colormap(ap.colormap('PWG',[],1.5));

    x1 = ap.groupfun(@sum,str_mua_cat(curr_use_ctx,:,use_align),grp_idx,[]);
    x1_norm = (x1-nanmean(x1(:,1:500),2))./nanmean(x1(:,1:500),2);
    x11 = max(x1_norm(:,use_t(1):use_t(2)),[],2);

    x1_norm_grp = ap.groupfun(@mean,x1_norm,grp(:,2),[]);
    nexttile(h); colororder(ap.colormap('BKR',7));
    plot(x1_norm_grp(ismember(ld_unique,-3:3),:)','linewidth',2);

    nexttile(h); hold on;
    arrayfun(@(x) plot(grp(grp(:,1)==x,2), x11(grp(:,1)==x)),unique(grp(:,1)));
    xlabel('Learned day');
    ylabel('MUA')

    [m,e] = grpstats(x11,grp(:,2),{'mean','sem'});
    errorbar(ld_unique,m,e,'k');
    drawnow;

    %%% TESTING (for stats)
    ld_unique = unique(grp(:,2));
    act_grid = nan(max(grp(:,1)),length(ld_unique));
    for curr_entry = 1:size(x11,1)
        act_grid(grp(curr_entry,1),ld_unique==grp(curr_entry,2)) = ...
            x11(curr_entry);
    end
    p = signrank(nanmean(act_grid(:,ld_unique <= -2),2),act_grid(:,ld_unique == -1));
    fprintf('Kidx %d: p = %.3f\n',curr_kidx,p);
    %%%%%

end
linkaxes(h.Children(2:2:end));


% Plot responsive cell fraction by k-means cluster
figure;h = tiledlayout(1,n_k);

stim_responsive_cells_cat = cell2mat(vertcat(stim_responsive_cells{use_animals}));
for curr_kidx = 1:n_k

    curr_use_ctx = kidx == curr_kidx;
    [grp,~,grp_idx] = unique(animal_ld_idx(curr_use_ctx,1:2),'rows');

    stim_cells_combine = ap.groupfun(@sum,stim_responsive_cells_cat(curr_use_ctx,:),grp_idx,[]);
    stim_cells_frac  = stim_cells_combine(:,1)./stim_cells_combine(:,2);

    nexttile(h); hold on;
    arrayfun(@(x) plot(grp(grp(:,1)==x,2), stim_cells_frac(grp(:,1)==x)),unique(grp(:,1)));
    xlabel('Learned day');
    ylabel('Frac stim-responsive')

    [m,e] = grpstats(stim_cells_frac,grp(:,2),{'mean','sem'});
    errorbar(unique(grp(:,2)),m,e,'k');
    drawnow;

    %%% TESTING (for stats)
    ld_unique = unique(grp(:,2));
    act_grid = nan(max(grp(:,1)),length(ld_unique));
    for curr_entry = 1:length(stim_cells_frac)
        act_grid(grp(curr_entry,1),ld_unique==grp(curr_entry,2)) = ...
            stim_cells_frac(curr_entry);
    end
    p = signrank(nanmean(act_grid(:,ld_unique <= -2),2),act_grid(:,ld_unique == -1));
    fprintf('Kidx %d: p = %.3f\n',curr_kidx,p);
    %%%%%

end

% Cortical prediction (k-kmeans)
day_mua_cat = cell2mat(cellfun(@cell2mat,day_mua_all,'uni',false));
day_mua_predicted_cat = cell2mat(cellfun(@cell2mat,day_mua_predicted_all,'uni',false));

use_align = 3;
figure;h = tiledlayout(2,n_k,'TileIndexing','columnmajor');
for curr_kidx = 1:n_k

    curr_use_ctx = kidx == curr_kidx;
    [grp,~,grp_idx] = unique(animal_ld_idx(curr_use_ctx,1:2),'rows');
    % Measured
    x1 = ap.groupfun(@sum,day_mua_cat(curr_use_ctx,:,use_align),grp_idx,[]);
    x1_norm = (x1-nanmean(x1(:,1:500),2))./nanmean(x1(:,1:500),2);
    x11 = max(x1_norm(:,use_t(1):use_t(2)),[],2);

    nexttile(h); hold on;
    arrayfun(@(x) plot(grp(grp(:,1)==x,2), x11(grp(:,1)==x)),unique(grp(:,1)));
    xlabel('Learned day');
    ylabel('MUA')
    title('Striatum');

    [m,e] = grpstats(x11,grp(:,2),{'mean','sem'});
    errorbar(unique(grp(:,2)),m,e,'k');
    drawnow;

    % Cortex-predicted
    x1 = ap.groupfun(@sum,day_mua_predicted_cat(curr_use_ctx,:,use_align),grp_idx,[]);
    x1_norm = (x1-nanmean(x1(:,1:500),2))./nanmean(x1(:,1:500),2);
    x11 = max(x1_norm(:,use_t(1):use_t(2)),[],2);

    nexttile(h); hold on;
    arrayfun(@(x) plot(grp(grp(:,1)==x,2), x11(grp(:,1)==x)),unique(grp(:,1)));
    xlabel('Learned day');
    ylabel('MUA')
    title('Cortex-predicted');

    [m,e] = grpstats(x11,grp(:,2),{'mean','sem'});
    errorbar(unique(grp(:,2)),m,e,'k');
    drawnow;

end

% Plot single animal over days
plot_animal = 2;
figure;
h = tiledlayout(1,length(day_mua_all{plot_animal}));
title(h,animals{plot_animal});
for curr_day = 1:length(day_mua_all{plot_animal})
    curr_activity = day_mua_all{plot_animal}{curr_day};

    softnorm = 1;
    curr_activity_norm = (curr_activity-nanmean(curr_activity(:,200:400,:),2))./ ...
        (nanmean(curr_activity(:,200:400,:),2)+softnorm);

    nexttile; hold on;
    spacing = 3;
    col = lines(size(curr_activity,3));
    for curr_align = 1:size(curr_activity_norm,3)
        AP_stackplot(curr_activity_norm(:,:,curr_align)',[],spacing,[],col(curr_align,:));
    end
    title(ld{plot_animal}(curr_day));

end
linkaxes(h.Children);


%% Trial data (striatum)

% Load data
am_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(am_data_path,'bhv.mat'));
load(fullfile(am_data_path,'ctx_maps_task.mat'));
load(fullfile(am_data_path,'trial_data_passive.mat'));
% load(fullfile(am_data_path,'trial_data_task.mat'));

% (time not saved: re-create)
t = -0.5:1/50:1;

% Convert cortex maps V to px
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);
ctx_map_px = cellfun(@(x) cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),x,'uni',false),ctx_map_all,'uni',false);

% Concatenate data
ctx_map_cat = cell2mat(permute(vertcat(ctx_map_px{:}),[2,3,1]));
ctx_expl_var_cat = cell2mat(vertcat(ctx_expl_var_all{:}));

% Get indicies for animal / learned day / depth
ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day}','uni',false);

animal_ld_idx_cell = cellfun(@(animal,ld,data) cellfun(@(ld,data) ...
    [repmat(animal,size(data,2),1),repmat(ld,size(data,2),1),(1:size(data,2))'], ...
    num2cell(ld),data,'uni',false), ...
    num2cell(1:length(ctx_map_all))',ld,ctx_map_all,'uni',false);
animal_ld_idx = cell2mat(vertcat(animal_ld_idx_cell{:}));

% Concatenate trial data
striatum_mua_cat = cell2mat(cellfun(@(x)  ...
    permute(reshape(permute(x,[2,1,3]),size(x,2),[]),[2,1,3]), ...
    vertcat(trial_data_all.striatum_mua{:}),'uni',false));

% Group area by cortical map k-means
n_k = 3;
kidx = nan(size(ctx_map_cat,3),1);
use_maps = squeeze(any(ctx_map_cat,[1,2]));
[kidx(use_maps),kmeans_map] = ...
    kmeans(reshape(ctx_map_cat(:,:,use_maps),[],sum(use_maps))', ...
    n_k,'distance','correlation','replicates',5);
kmeans_map = reshape(kmeans_map',[size(ctx_map_cat,[1,2]),n_k]);

figure;tiledlayout(n_k,1,'tilespacing','none');
colormap(ap.colormap('PWG',[],1.5));
for i = 1:n_k
    nexttile;
    imagesc(kmeans_map(:,:,i));
    axis image off;
    clim(max(abs(clim)).*[-1,1]);
    ap.wf_draw('ccf','k');
end

kidx_day = mat2cell(kidx,cellfun(@(x) size(x,3),vertcat(ctx_map_px{:})),1);

%%%%% TEMP: fixing number mismatch (end days with no striatum, cell size
%%%%% wasn't initialized first)
for curr_animal = 1:length(trial_data_all.striatum_mua)
    if length(trial_data_all.striatum_mua{curr_animal}) ~= ...
            length(bhv(curr_animal).learned_day)

        trial_data_all.striatum_mua{curr_animal}{ ...
            length(bhv(curr_animal).learned_day)} = [];

    end
end
%%%%%%%%%%

% Make indicies per trial
day_animal = cell2mat(cellfun(@(x,y) repmat(x,length(y),1), ...
    num2cell(1:length(bhv)),{bhv.learned_day},'uni',false)');
trial_animal = cell2mat(cellfun(@(x,y) repmat(x,prod(size(y,[1,3])),1), ...
    num2cell(day_animal),vertcat(trial_data_all.striatum_mua{:}),'uni',false));

trial_ld = cell2mat(cellfun(@(x,y) repmat(x,prod(size(y,[1,3])),1), ...
    num2cell(vertcat(ld{:})),vertcat(trial_data_all.striatum_mua{:}),'uni',false));

trial_idx = cell2mat(cellfun(@(x) reshape(repmat((1:size(x,1))',size(x,3),1),[],1), ...
    vertcat(trial_data_all.striatum_mua{:}),'uni',false));

trial_frac_idx = cell2mat(cellfun(@(x) reshape(repmat((1:size(x,1))'./size(x,1),size(x,3),1),[],1), ...
    vertcat(trial_data_all.striatum_mua{:}),'uni',false));

% (passive)
if isfield(trial_data_all,'stim')
    trial_stim = cell2mat(cellfun(@(x,y) repmat(x(1:size(y,1)),size(y,3),1), ...
        vertcat(trial_data_all.stim{:}), ...
        vertcat(trial_data_all.striatum_mua{:}),'uni',false));
end
if isfield(trial_data_all,'wheel_move')
    trial_wheel = cell2mat(cellfun(@(x,y) repmat(x(1:size(y,1),:),size(y,3),1), ...
        vertcat(trial_data_all.wheel_move{:}), ...
        vertcat(trial_data_all.striatum_mua{:}),'uni',false));
end

% (task)
if isfield(trial_data_all,'stim_to_move')
    trial_rxn = cell2mat(cellfun(@(x,y) repmat(x(1:size(y,1)),size(y,3),1), ...
        vertcat(trial_data_all.stim_to_move{:}), ...
        vertcat(trial_data_all.striatum_mua{:}),'uni',false));
end
if isfield(trial_data_all,'outcome')
    trial_outcome = cell2mat(cellfun(@(x,y) repmat(x(1:size(y,1)),size(y,3),1), ...
        vertcat(trial_data_all.outcome{:}), ...
        vertcat(trial_data_all.striatum_mua{:}),'uni',false));
end

trial_depth = cell2mat(cellfun(@(x) reshape(repmat(1:size(x,3),size(x,1),1),[],1), ...
    vertcat(trial_data_all.striatum_mua{:}),'uni',false));

trial_kidx = cell2mat(cellfun(@(x,y) reshape(repmat(x',size(y,1),1),[],1), ...
    kidx_day,vertcat(trial_data_all.striatum_mua{:}),'uni',false));

% (PLOT PASSIVE)
plot_ld = [-3:3];
figure; h = tiledlayout(n_k,length(plot_ld),'TileSpacing','tight');
colormap(ap.colormap('WK',[],1.5));

for curr_kidx = 1:n_k

    use_trials = trial_stim == 90 & ~any(trial_wheel(:,t>0 & t<=0.3),2) & trial_kidx == curr_kidx;
    [grp,~,grp_idx] = unique([trial_animal(use_trials),trial_ld(use_trials),trial_idx(use_trials)],'rows');

    curr_striatum_mua_sum = ap.groupfun(@sum,striatum_mua_cat(use_trials,:),grp_idx,[]);

    softnorm = 20;
    curr_striatum_mua_norm = (curr_striatum_mua_sum-nanmean(curr_striatum_mua_sum(:,t<0),2))./ ...
        (nanmean(curr_striatum_mua_sum(:,t<0),2) + softnorm);

    for curr_ld = plot_ld
        nexttile(h);
        imagesc(imgaussfilt(curr_striatum_mua_norm(grp(:,2) == curr_ld,:),[3,1]));
        colormap(ap.colormap('WK',[],1.5));
        clim([0,2]);
        title(sprintf('LD %d',curr_ld));
    end

    % (Plot average across animals)
    [trial_grp,~,trial_grp_idx] = unique(grp(:,1:2),'rows');
    curr_striatum_trialavg = ap.groupfun(@mean,curr_striatum_mua_norm,trial_grp_idx,[]);

    [ld_grp,~,ld_grp_idx] = unique(trial_grp(:,2));
    curr_striatum_ldavg = ap.groupfun(@mean,curr_striatum_trialavg,ld_grp_idx,[]);

    figure;
    plot(t,curr_striatum_ldavg(ismember(ld_grp,plot_ld),:)','linewidth',2);
    colororder(gca,ap.colormap('BKR',length(plot_ld)));

end

% (PLOT TASK)
plot_ld = [-3:3];

figure; h1 = tiledlayout(n_k,length(plot_ld),'TileSpacing','tight');
colormap(ap.colormap('WK',[],1.5));

figure; h2 = tiledlayout(n_k,1,'TileSpacing','tight');

for curr_kidx = 1:n_k

    use_trials = find(trial_kidx == curr_kidx);
    [grp,grp_first_idx,grp_idx] = unique([trial_animal(use_trials),trial_ld(use_trials),trial_idx(use_trials)],'rows');

    curr_striatum_mua_sum = ap.groupfun(@sum,striatum_mua_cat(use_trials,:),grp_idx,[]);

    softnorm = 20;
    curr_striatum_mua_norm = (curr_striatum_mua_sum-nanmean(curr_striatum_mua_sum(:,t<0),2))./ ...
        (nanmean(curr_striatum_mua_sum(:,t<0),2) + softnorm);

    curr_rxn = trial_rxn(use_trials(grp_first_idx));

    for curr_ld = plot_ld
        plot_trials = find(grp(:,2) == curr_ld);
        [~,sort_idx] = sort(curr_rxn(plot_trials));

        nexttile(h1);
        imagesc(t,[],imgaussfilt(curr_striatum_mua_norm(plot_trials(sort_idx),:),[5,1]));
        clim([0,2]); hold on;
        xline(0,'r','linewidth',2);
        plot(curr_rxn(plot_trials(sort_idx)),1:length(plot_trials),'b','linewidth',2)
        title(sprintf('LD %d',curr_ld));
    end

    % (Plot average across animals)
    [trial_grp,~,trial_grp_idx] = unique(grp(:,1:2),'rows');
    curr_striatum_trialavg = ap.groupfun(@mean,curr_striatum_mua_norm,trial_grp_idx,[]);

    [ld_grp,~,ld_grp_idx] = unique(trial_grp(:,2));
    curr_striatum_ldavg = ap.groupfun(@mean,curr_striatum_trialavg,ld_grp_idx,[]);

    nexttile(h2);
    plot(t,curr_striatum_ldavg(ismember(ld_grp,plot_ld),:)','linewidth',2);
    colororder(gca,ap.colormap('BKR',length(plot_ld)));

end


%% Trial data (widefield)

% Load data
am_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(am_data_path,'bhv.mat'));
load(fullfile(am_data_path,'ctx_maps_task.mat'));
load(fullfile(am_data_path,'trial_data_passive.mat'));
% load(fullfile(am_data_path,'trial_data_task.mat'));

% (time not saved: re-create)
t = -0.5:1/50:1;

% Convert cortex maps V to px
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);
ctx_map_px = cellfun(@(x) cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),x,'uni',false),ctx_map_all,'uni',false);

% Concatenate data
ctx_map_cat = cell2mat(permute(vertcat(ctx_map_px{:}),[2,3,1]));
ctx_expl_var_cat = cell2mat(vertcat(ctx_expl_var_all{:}));

% Get indicies for animal / learned day / depth
ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day}','uni',false);

animal_ld_idx_cell = cellfun(@(animal,ld,data) cellfun(@(ld,data) ...
    [repmat(animal,size(data,2),1),repmat(ld,size(data,2),1),(1:size(data,2))'], ...
    num2cell(ld),data,'uni',false), ...
    num2cell(1:length(ctx_map_all))',ld,ctx_map_all,'uni',false);
animal_ld_idx = cell2mat(vertcat(animal_ld_idx_cell{:}));

% Concatenate trial data (and subtract baseline)
wf_cat = cell2mat(vertcat(trial_data_all.widefield{:}));
wf_cat = wf_cat - nanmean(wf_cat(:,t<0,:),2);

n_vs = size(wf_cat,3);

% Make indicies per trial
day_animal = cell2mat(cellfun(@(x,y) repmat(x,length(y),1), ...
    num2cell(1:length(bhv)),{bhv.learned_day},'uni',false)');
trial_animal = cell2mat(cellfun(@(x,y) repmat(x,size(y,1),1), ...
    num2cell(day_animal),vertcat(trial_data_all.widefield{:}),'uni',false));

trial_ld = cell2mat(cellfun(@(x,y) repmat(x,size(y,1),1), ...
    num2cell(vertcat(ld{:})),vertcat(trial_data_all.widefield{:}),'uni',false));

trial_idx = cell2mat(cellfun(@(x) repmat((1:size(x,1))',1), ...
    vertcat(trial_data_all.widefield{:}),'uni',false));

% (passive)
if isfield(trial_data_all,'stim')
    trial_stim = cell2mat(vertcat(trial_data_all.stim{:}));
    trial_wheel = cell2mat(vertcat(trial_data_all.wheel_move{:}));
end

% (task)
if isfield(trial_data_all,'stim_to_move')
    trial_rxn = cell2mat(vertcat(trial_data_all.stim_to_move{:}));
end
if isfield(trial_data_all,'outcome')
    trial_outcome = cell2mat(vertcat(trial_data_all.outcome{:}));
end

% Plot average movie for drawing ROI
avg_px = plab.wf.svd2px(U_master(:,:,1:n_vs), ...
    permute(nanmean(wf_cat,1),[3,2,1]));
ap.imscroll(avg_px,t);
axis image; 
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG',[],1.5));
ap.wf_draw('ccf','k');

% (draw roi)
roi_trace = permute(ap.wf_roi(U_master(:,:,1:n_vs), ...
    permute(wf_cat,[3,2,1]),[],[],roi.mask),[3,2,1]);

% (PLOT TASK)
plot_ld = [-3:3];
use_trials = find(trial_rxn > 0.3);find(true(size(trial_animal,1),1));

figure; h = tiledlayout(1,length(plot_ld),'TileSpacing','tight');
colormap(ap.colormap('PWG',[],1.5));
for curr_ld = plot_ld
    plot_trials = use_trials(trial_ld(use_trials) == curr_ld);
    [~,sort_idx] = sort(trial_rxn(plot_trials));

    nexttile(h);
    imagesc(t,[],imgaussfilt(roi_trace(plot_trials(sort_idx),:),[5,1]));
    clim([-1,1].*0.015); hold on;
    xline(0,'r','linewidth',2);
    plot(trial_rxn(plot_trials(sort_idx)),1:length(plot_trials),'b','linewidth',2)
    title(sprintf('LD %d',curr_ld));
end

% (Plot average across animals)
[grp,grp_first_idx,grp_idx] = unique([trial_animal(use_trials),trial_ld(use_trials),trial_idx(use_trials)],'rows');

[trial_grp,~,trial_grp_idx] = unique(grp(:,1:2),'rows');
roi_trialavg = ap.groupfun(@mean,roi_trace(use_trials,:),trial_grp_idx,[]);

[ld_grp,~,ld_grp_idx] = unique(trial_grp(:,2));
roi_ldavg = ap.groupfun(@mean,roi_trialavg,ld_grp_idx,[]);

figure;
plot(t,roi_ldavg(ismember(ld_grp,plot_ld),:)','linewidth',2);
colororder(gca,ap.colormap('BKR',length(plot_ld)));


% (plot average movie of specific trials)
use_trials = trial_ld < -1;
% use_trials = ~any(trial_wheel(:,t>0 & t<=0.5),2) & trial_stim == 90 & trial_ld == 0;
avg_px = plab.wf.svd2px(U_master(:,:,1:n_vs), ...
    permute(nanmean(wf_cat(use_trials,:,:),1),[3,2,1]));
ap.imscroll(avg_px,t);
axis image; 
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG',[],1.5));
ap.wf_draw('ccf','k');


% (PLOT PASSIVE)
plot_ld = -3:3;

% use_trials = ~any(trial_wheel(:,t>0 & t<=0.3),2) & ...
%     ismember(trial_ld,plot_ld) & trial_stim == 90;
use_trials = ismember(trial_ld,plot_ld);

figure; h = tiledlayout(1,length(plot_ld),'TileSpacing','tight');
colormap(ap.colormap('PWG',[],1.5));

for curr_ld = plot_ld
    plot_trials = use_trials & trial_ld == curr_ld;

    nexttile(h);
    imagesc(t,[],imgaussfilt(roi_trace(plot_trials,:),[1,1]));
    clim([-1,1].*0.015); hold on;
    xline(0,'r','linewidth',2);
    title(sprintf('LD %d',curr_ld));
end

% (Plot average across animals)
[trial_grp,~,trial_grp_idx] = unique([trial_animal(use_trials),trial_ld(use_trials)],'rows');
roi_trialavg = ap.groupfun(@mean,roi_trace(use_trials,:),trial_grp_idx,[]);

[ld_grp,~,ld_grp_idx] = unique(trial_grp(:,2));
roi_ldavg = ap.groupfun(@mean,roi_trialavg,ld_grp_idx,[]);

figure;
plot(t,roi_ldavg(ismember(ld_grp,plot_ld),:)','linewidth',2);
colororder(gca,ap.colormap('BKR',length(plot_ld)));

% %%%%%%%%%%%
% % (plot ROI across days, sorted)
% figure; h = tiledlayout(1,length(plot_ld),'TileSpacing','tight');
% colormap(ap.colormap('PWG',[],1.5));
% 
% for curr_ld = plot_ld
%     plot_trials = find(use_trials & trial_ld == curr_ld & trial_rxn > 0.2);
%     [~,x] = sort(trial_rxn(plot_trials));
% %     [~,x] = sort(nanmean(roi_trace(plot_trials,t>0&t<0.3),2),'descend');
%     nexttile(h);
%     imagesc(t,[],imgaussfilt(roi_trace(plot_trials(x),:),[1,1]));
%     clim([-1,1].*0.015); hold on;
%     xline(0,'r','linewidth',2);
%     title(sprintf('LD %d',curr_ld));
% end


%% Trial data (striatum + widefield)

% Load data
am_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(am_data_path,'bhv.mat'));
load(fullfile(am_data_path,'ctx_maps_task.mat'));
% load(fullfile(am_data_path,'trial_data_passive.mat'));
load(fullfile(am_data_path,'trial_data_task.mat'));

% (time not saved: re-create)
t = -0.5:1/50:1;

% Convert cortex maps V to px
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);
ctx_map_px = cellfun(@(x) cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),x,'uni',false),ctx_map_all,'uni',false);

% Concatenate data
ctx_map_cat = cell2mat(permute(vertcat(ctx_map_px{:}),[2,3,1]));
ctx_expl_var_cat = cell2mat(vertcat(ctx_expl_var_all{:}));

% Get indicies for animal / learned day / depth
ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day}','uni',false);

animal_ld_idx_cell = cellfun(@(animal,ld,data) cellfun(@(ld,data) ...
    [repmat(animal,size(data,2),1),repmat(ld,size(data,2),1),(1:size(data,2))'], ...
    num2cell(ld),data,'uni',false), ...
    num2cell(1:length(ctx_map_all))',ld,ctx_map_all,'uni',false);
animal_ld_idx = cell2mat(vertcat(animal_ld_idx_cell{:}));

% Make indicies per trial
day_animal = cell2mat(cellfun(@(x,y) repmat(x,length(y),1), ...
    num2cell(1:length(bhv)),{bhv.learned_day},'uni',false)');
trial_animal = cell2mat(cellfun(@(x,y) repmat(x,size(y,1),1), ...
    num2cell(day_animal),vertcat(trial_data_all.widefield{:}),'uni',false));

trial_ld = cell2mat(cellfun(@(x,y) repmat(x,size(y,1),1), ...
    num2cell(vertcat(ld{:})),vertcat(trial_data_all.widefield{:}),'uni',false));

trial_idx = cell2mat(cellfun(@(x) repmat((1:size(x,1))',1), ...
    vertcat(trial_data_all.widefield{:}),'uni',false));

trial_frac_idx = cell2mat(cellfun(@(x) repmat((1:size(x,1))'./size(x,1),1), ...
    vertcat(trial_data_all.widefield{:}),'uni',false));

% (passive)
if isfield(trial_data_all,'stim')
    trial_stim = cell2mat(vertcat(trial_data_all.stim{:}));
end
if isfield(trial_data_all,'wheel_move')
    trial_wheel = cell2mat(vertcat(trial_data_all.wheel_move{:}));
end
if isfield(trial_data_all,'wheel_velocity')
    trial_wheel_vel = cell2mat(vertcat(trial_data_all.wheel_velocity{:}));
end
% (task)
if isfield(trial_data_all,'stim_to_move')
    trial_rxn = cell2mat(vertcat(trial_data_all.stim_to_move{:}));
end
if isfield(trial_data_all,'outcome')
    trial_outcome = cell2mat(vertcat(trial_data_all.outcome{:}));
end

% Widefield: concatenate trial data (and subtract baseline)
wf_cat = cell2mat(vertcat(trial_data_all.widefield{:}));
wf_cat = wf_cat - nanmean(wf_cat(:,t<0,:),2);

n_vs = size(wf_cat,3);

% Striatum: concatenate trial data, combine/normalize fraom kidx (above)

%%%%% TEMP: fixing number mismatch (end days with no striatum, cell size
%%%%% wasn't initialized first)
trial_data_all.striatum_mua{2}{4} = [];
trial_data_all.striatum_mua{8}{8} = [];
%%%%%%%%%%

use_str_depth = kidx == 1;
use_str_depth_day = mat2cell(use_str_depth,cellfun(@(x) size(x,3),vertcat(ctx_map_px{:})),1);

softnorm = 20;
str_cat = cellfun(@(x,use_depth) (sum(x(:,:,use_depth),3) - nanmean(sum(x(:,t<0,use_depth),3),2)) ./ ...
    (nanmean(sum(x(:,t<0,use_depth),3),2) + softnorm), ...
    vertcat(trial_data_all.striatum_mua{:}),use_str_depth_day,'uni',false, ...
    'ErrorHandler',@(varargin) []);

% (create NaNs for days w/o striatum: this is shit but it's quick for now)
x = vertcat(trial_data_all.widefield{:});
str_cat(cellfun(@isempty,str_cat)) = cellfun(@(x) nan(size(x,1),size(x,2)),x(cellfun(@isempty,str_cat)),'uni',false);
str_cat = cell2mat(str_cat);
clear x

% (also dirty for now: days w/o visual striatum are all zeros, NaN-out)
str_cat(all(str_cat==0,2),:) = NaN;


%%%%% Split by one modality, look at other

% Plot average movie for drawing ROI
avg_px = plab.wf.svd2px(U_master(:,:,1:n_vs), ...
    permute(nanmean(wf_cat,1),[3,2,1]));
ap.imscroll(avg_px,t);
axis image; 
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG',[],1.5));
ap.wf_draw('ccf','k');

% (after drawing ROI, grab fluorescence)
roi_trace = permute(ap.wf_roi(U_master(:,:,1:n_vs), ...
    permute(wf_cat,[3,2,1]),[],[],roi.mask),[3,2,1]);

% (split data within animal)
use_t = t >= 0.05 & t <= 0.15;

% use_trials_all = trial_stim == 90 &  ~any(trial_wheel(:,t>0 & t<=0.3),2) & trial_ld == -1;
use_trials_all = trial_ld <= -1;

n_split = 2;
trial_split = nan(size(trial_ld));
for curr_animal = unique(trial_animal)'
    curr_animal_trials = trial_animal == curr_animal;
    curr_trials = use_trials_all & curr_animal_trials;

    trial_split(curr_trials) = discretize(tiedrank( ...
        nanmean(str_cat(curr_trials,use_t),2))./ ...
        sum(curr_trials),linspace(0,1,n_split+1));
end

use_trials_split = ~isnan(trial_split);

str_grp_mean = ap.groupfun(@nanmean,str_cat(use_trials_split,:), ...
    trial_split(use_trials_split),[]);

curr_wf_split = ap.groupfun(@nanmean,wf_cat(use_trials_split,:,:), ...
    trial_split(use_trials_split),[],[]);
wf_roi_grp_mean = permute(ap.wf_roi(U_master(:,:,1:n_vs), ...
    permute(curr_wf_split,[3,2,1]),[],[],roi.mask),[3,2,1]);

figure; tiledlayout(1,2); 
nexttile; hold on; set(gca,'colororder',ap.colormap('KR',n_split));
plot(t,str_grp_mean','linewidth',2);
nexttile; hold on; set(gca,'colororder',ap.colormap('KR',n_split));
plot(t,wf_roi_grp_mean','linewidth',2);
 
curr_wf_split_px = plab.wf.svd2px(U_master(:,:,1:n_vs),permute(curr_wf_split,[3,2,1]));
ap.imscroll(curr_wf_split_px,t);
axis image; 
clim(max(abs(clim)).*[-1,1]*0.8);
colormap(ap.colormap('PWG',[],1.5));
ap.wf_draw('ccf','k');

wf_max_use_t = t > 0 & t < 0.2;
ap.imscroll(squeeze(max(curr_wf_split_px(:,:,wf_max_use_t,:),[],3)));
axis image; 
clim(max(abs(clim)).*[-1,1]*0.8);
colormap(ap.colormap('PWG',[],1.5));
ap.wf_draw('ccf','k');


% avg str by animal 
roi_trace = permute(ap.wf_roi(U_master(:,:,1:n_vs), ...
    permute(wf_cat,[3,2,1]),[],[],roi.mask),[3,2,1]);

[grp,~,tr_idx] = unique([trial_animal(use_trials_split),trial_split(use_trials_split)],'rows');

figure; tiledlayout(1,3);
nexttile; hold on; set(gca,'colororder',ap.colormap('KR',n_split));
x1 = ap.groupfun(@nanmean,str_cat(use_trials_split,:),tr_idx,[]);
x2 = ap.groupfun(@nanmean,x1,grp(:,2),[]);
plot(t,x2','linewidth',2);
title('Striatum');

nexttile; hold on; set(gca,'colororder',ap.colormap('KR',n_split));
x1 = ap.groupfun(@nanmean,roi_trace(use_trials_split,:),tr_idx,[]);
x2 = ap.groupfun(@nanmean,x1,grp(:,2),[]);
plot(t,x2','linewidth',2);
title('Cortex');

nexttile; hold on; set(gca,'colororder',ap.colormap('KR',n_split));
x1 = ap.groupfun(@nanmean,trial_wheel_vel(use_trials_split,:),tr_idx,[]);
x2 = ap.groupfun(@nanmean,x1,grp(:,2),[]);
plot(t,x2','linewidth',2);
title('Wheel');




% Regress cortex to spikes (using trials)
for grp = 1:4

    switch grp
        case 1
            use_trials = trial_ld < -1;
        case 2
            use_trials = trial_ld == -1;
        case 3
            use_trials = trial_ld == 0;
        case 4
            use_trials = trial_ld > 0;
    end

    kernel_frames = -10:10;
    lambda = 100;
    zs = [false,false];
    cvfold = 1;
    return_constant = false;
    use_constant = true;

    use_t = t > 0 & t < 0.5;

    discontinuities = false(sum(use_trials),sum(use_t));
    discontinuities(:,end) = true;
    discontinuities = reshape(discontinuities',[],1)';

    [k,predicted_spikes,explained_var] = ...
        ap.regresskernel( ...
        reshape(permute(wf_cat(use_trials,use_t,:),[2,1,3]),[],size(wf_cat,3))', ...
        reshape(str_cat(use_trials,use_t)',[],1)', ...
        kernel_frames, ...
        lambda,zs,cvfold,return_constant,use_constant,discontinuities);

    k_px = plab.wf.svd2px(U_master(:,:,1:length(k)),k);

    % figure;imagesc(k_px);
    ap.imscroll(k_px);
    axis image
    clim(max(abs(clim)).*[-1,1]*0.8);
    colormap(ap.colormap('PWG',[],1.5));
    ap.wf_draw('ccf','k');
    drawnow; 

end


% splitting by rxn, then striatum
figure; h = tiledlayout(2,2,'TileIndexing','ColumnMajor');       
for rxn = 1:2

    switch rxn
        case 1
            use_rxn = trial_rxn > 0.0 & trial_rxn <= 0.2;
        case 2
            use_rxn = trial_rxn > 0.2 & trial_rxn <= 0.4;
    end

    use_trials_all = use_rxn;
%     use_trials_all = ~any(trial_wheel(:,t>0 & t<=0.3),2) & trial_stim == 90;

    use_t = t >= 0.05 & t <= 0.15;

    n_split = 3;
    trial_split = nan(size(trial_ld));
    for curr_animal = unique(trial_animal)'
        curr_animal_trials = trial_animal == curr_animal;
        curr_trials = use_trials_all & curr_animal_trials;

        trial_split(curr_trials) = discretize(tiedrank( ...
            mean(str_cat(curr_trials,use_t),2))./ ...
            sum(curr_trials),linspace(0.1,0.9,n_split+1));
    end

    use_trials_split = ~isnan(trial_split) & trial_ld < 0;

    str_grp_mean = ap.groupfun(@nanmean,str_cat(use_trials_split,:), ...
        trial_split(use_trials_split),[]);

    curr_wf_split = ap.groupfun(@nanmean,wf_cat(use_trials_split,:,:), ...
        trial_split(use_trials_split),[],[]);
    wf_roi_grp_mean = squeeze(ap.wf_roi(U_master(:,:,1:n_vs),permute(curr_wf_split,[3,2,1]),[],[],roi.mask));

    nexttile(h); hold on; set(gca,'colororder',ap.colormap('KR',n_split));
    plot(t,str_grp_mean','linewidth',2);
    nexttile(h); hold on; set(gca,'colororder',ap.colormap('KR',n_split));
    plot(t,wf_roi_grp_mean,'linewidth',2);

    curr_wf_split_px = plab.wf.svd2px(U_master(:,:,1:n_vs),permute(curr_wf_split,[3,2,1]));
    ap.imscroll(curr_wf_split_px,t);
    axis image;
    clim(max(abs(clim)).*[-1,1]*0.8);
    colormap(ap.colormap('PWG',[],1.5));
    ap.wf_draw('ccf','k');

    figure;imagesc(reshape(mean(curr_wf_split_px(:,:,use_t,:),3), ...
        size(curr_wf_split_px,1),[]));
    axis image off;
    clim(max(abs(clim)).*[-1,1]*0.8);
    colormap(ap.colormap('PWG',[],1.5));

    drawnow;

end
linkaxes(h.Children(1:2:end));
linkaxes(h.Children(2:2:end));



% Plot trials split by fraction within day
n_trial_frac = 4;
str_cat_ld_trialsplit = nan(7,length(t),n_trial_frac);
figure; tiledlayout(1,n_trial_frac);
for curr_trial_frac = 1:n_trial_frac
    use_trials = ismember(trial_ld,-3:3) & ...
        trial_frac_idx > (curr_trial_frac-1)/n_trial_frac & ...
        trial_frac_idx <= curr_trial_frac/n_trial_frac;
    [grp,~,tr_idx] = unique([trial_animal(use_trials),trial_ld(use_trials)],'rows');

    str_cat_trialgrp = ap.groupfun(@nanmean,str_cat(use_trials,:),tr_idx,[]);
    str_cat_ld = ap.groupfun(@nanmean,str_cat_trialgrp,grp(:,2),[]);
    
    nexttile;
    plot(t,str_cat_ld');
    colororder(gca,ap.colormap('BKR',7));
    title(sprintf('Trial %d/%d',curr_trial_frac,n_trial_frac));

    str_cat_ld_trialsplit(:,:,curr_trial_frac) = str_cat_ld;
end

use_t = t > 0.05 & t < 0.15;
str_cat_ld_trialsplit_tmean = permute(nanmean(str_cat_ld_trialsplit(:,use_t,:),2),[1,3,2]);
figure;plot(unique(grp(:,2)),str_cat_ld_trialsplit_tmean);


%% Utility: count days with given kidx

kidx_set = reshape(unique(kidx(~isnan(kidx))),1,[]);
kidx_day = mat2cell(kidx,cellfun(@(x) size(x,2), vertcat(ctx_map_all{:})));

kidx_day_setmatch = cell2mat(cellfun(@(kidx) any(kidx == kidx_set,1), ...
    kidx_day,'uni',false));

[kidx_ld_n,ld_grp] = grpstats(kidx_day_setmatch,vertcat(ld{:}),{'sum','gname'});
ld_grp = cellfun(@str2num,ld_grp);

figure;
tiledlayout(1,2);

nexttile;
imagesc(AP_padcatcell(kidx_day));
xlabel('Day');
ylabel('Depth');
title('Kidx');
clim([-0.5,max(kidx_set)]);

nexttile;
plot(ld_grp,kidx_ld_n)
xlabel('Learned day');
ylabel('Kidx count');























