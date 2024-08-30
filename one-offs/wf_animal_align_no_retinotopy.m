% Create master L/R image (from AM data)
% (note - this contains the data to be aligned, alignment was overwritten)

am_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';
load(fullfile(am_data_path,'wf_passive.mat'));

% Load master U, convert V to px
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);

% Get total average passive widefield
V_cat = cat(4,day_V_all{:});
wf_passive_avg = plab.wf.svd2px(U_master,squeeze(nanmean(V_cat,4)));

% Get max L - max R
use_frames = 15:25;
wf_passive_max = squeeze(max(wf_passive_avg(:,:,use_frames,:),[],3));
wf_lr_master = wf_passive_max(:,:,1) - wf_passive_max(:,:,3);

clearvars -except wf_lr_master


%% Create L/R image and animal align

animal = 'AM015';

recordings = plab.find_recordings(animal,[],'lcr_passive');
wf_recordings = recordings(cellfun(@any,{recordings.widefield}));

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.03;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

wf_lr_day_mean = cell(size(wf_recordings));

for curr_recording = 1:length(wf_recordings)

    % Grab pre-load vars
    preload_vars = who;

    % Load data
    rec_day = recordings(curr_recording).day;
    rec_time = recordings(curr_recording).recording{end};

    load_parts = struct;
    load_parts.widefield = true;
    ap.load_recording;

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
    aligned_V_mean = nan(size(wf_U,3),length(t_centers),length(align_times));
    for curr_align = 1:length(align_times)
        peri_event_t = align_times{curr_align} + t_centers;

        aligned_V = interp1(wf_t,wf_V',peri_event_t);

        aligned_V_baselinesub = aligned_V - ...
            mean(aligned_V(:,t_centers < 0,:),2);

        aligned_V_mean(:,:,curr_align) = ...
            permute(nanmean(aligned_V_baselinesub,1),[3,2,1]);
    end
    aligned_px_mean = plab.wf.svd2px(wf_U,aligned_V_mean);
    wf_day_mean{curr_recording} = aligned_px_mean;

    ap.print_progress_fraction(curr_recording,length(wf_recordings));
end

use_frames = 15:25;
% (use all non-empty days)
wf_mean = nanmean(cat(5,wf_day_mean{cellfun(@(x) ~isempty(x),wf_day_mean)}),5);
% (use subset of days)
% wf_mean = nanmean(cat(5,wf_day_mean{5:7}),5);

wf_max = squeeze(max(wf_mean(:,:,use_frames,:),[],3));
wf_lr_unaligned = wf_max(:,:,1) - wf_max(:,:,3);

% Align animal to master
plab.wf.wf_align(wf_lr_unaligned,animal,[],'new_animal',wf_lr_master);

% NOTE:
% If the alignemnt doesn't look good, trying cutting out bad days
% (e.g. AM015 looked better using days 3-6)








