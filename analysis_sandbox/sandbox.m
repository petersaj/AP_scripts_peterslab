%% SANDBOX
%
% Temporary code

%% Get move times in task (put into load)

% Align wheel movement to stim onset
stimOn_times = photodiode_times(1:2:end);

[wheel_velocity,wheel_move] = AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);
surround_sample_rate = timelite.daq_info(1).rate;
wheel_window = [-0.5,2];
wheel_window_t = wheel_window(1):1/surround_sample_rate:wheel_window(2);
wheel_window_t_peri_event = stimOn_times + wheel_window_t;
stim_aligned_wheel = interp1(timelite.timestamps, ...
    wheel_velocity,wheel_window_t_peri_event);

% (set a threshold in speed and time for wheel movement)
thresh_displacement = 0.025;
time_over_thresh = 0.05; % ms over velocity threshold to count
samples_over_thresh = time_over_thresh.*surround_sample_rate;
wheel_over_thresh_fullconv = convn( ...
    abs(stim_aligned_wheel) > thresh_displacement, ...
    ones(1,samples_over_thresh)) >= samples_over_thresh;
wheel_over_thresh = wheel_over_thresh_fullconv(:,end-size(stim_aligned_wheel,2)+1:end);

[move_trial,wheel_move_sample] = max(wheel_over_thresh,[],2);
wheel_move_time = arrayfun(@(x) wheel_window_t_peri_event(x,wheel_move_sample(x)),1:size(wheel_window_t_peri_event,1))';
wheel_move_time(~move_trial) = NaN;

stim_to_move = wheel_move_time - stimOn_times;


%% Ephys raster

% (passive)
align_times_all = photodiode_times(1:2:end);
align_category_all = vertcat(trial_events.values.TrialStimX);
% (get only quiescent trials)
[wheel_velocity,wheel_move] = AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);
framerate = 30;
wheel_window = [0,0.5];
wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
wheel_window_t_peri_event = align_times_all + wheel_window_t;
event_aligned_move = interp1(timelite.timestamps, ...
    +wheel_move,wheel_window_t_peri_event,'previous');
quiescent_trials = ~any(event_aligned_move,2);
align_times = align_times_all(quiescent_trials);
align_category = align_category_all(quiescent_trials);

AP_cellraster(align_times,align_category);

% % (task)
% stimOn_times = photodiode_times(1:2:end);
% wheel_starts = timelite.timestamps(diff([0;wheel_move]) == 1);
% 
% % (get wheel starts when no stim on screen: not sure this works yet)
% wheel_starts_iti = ...
%     wheel_starts(interp1(timelite.timestamps, ...
%     +photodiode_thresh,wheel_starts,'previous') == 0);
% 
% [~,rxn_sort_idx] = sort(stim_to_move);
% AP_cellraster({stimOn_times,wheel_move_time,wheel_starts_iti,reward_times}, ...
%     {rxn_sort_idx,rxn_sort_idx,1:length(wheel_starts_iti),1:length(reward_times)});


%% Ephys sparse noise? 

stim_screen = discretize(noise_locations,[0,128,129,255],-1:1);
% Get stim times vector (x,y)
nY = size(stim_screen,1);
nX = size(stim_screen,2);
stim_times_grid = cell(nY,nX);
for x = 1:nX
    for y = 1:nY      
        align_stims = (stim_screen(y,x,1:end-1) == 0) & ...
            (stim_screen(y,x,2:end) == 0);
        align_times = stim_times(find(align_stims)+1);
        stim_times_grid{y,x} = align_times;
    end
end

% Vectorize stim times by y,x position
[stim_x,stim_y] = meshgrid(1:nX,1:nY);
stim_positions = cellfun(@(x,y,times) repmat([y,x],length(times),1), ...
    num2cell(stim_x),num2cell(stim_y),stim_times_grid,'uni',false);

% Pick unit
% use_spikes = spike_times_timeline(spike_templates == 276);
use_spikes = spike_times_timeline(spike_templates == 229);
% use_spikes = spike_times_timeline(ismember(spike_templates, ...
%     find(template_depths > 1700 & template_depths < 2300)));

% Get stim times vector (x,y)
stim_aligned_avg = cell(nY,nX);
raster_window = [-0.05,0.2];
smooth_size = 2;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
for x = 1:nX
    for y = 1:nY      
        
        psth_bin_size = 0.01;
        [psth,bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
            use_spikes,stim_times_grid{y,x}, ...
            raster_window, psth_bin_size);
        
        psth_smooth = conv2(psth,smWin,'same');
        stim_aligned_avg{y,x} = psth_smooth;
        
    end
end
stim_aligned_avg_cat = cell2mat(cellfun(@(x) permute(x,[1,3,2]),stim_aligned_avg,'uni',false));
AP_imscroll(stim_aligned_avg_cat,bins);
axis image;




%% Align widefield to event

% (passive)
align_times_all = photodiode_times(1:2:end);
align_category_all = vertcat(trial_events.values.TrialStimX);
% (get only quiescent trials)
[wheel_velocity,wheel_move] = AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);
framerate = 30;
wheel_window = [0,0.5];
wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
wheel_window_t_peri_event = align_times_all + wheel_window_t;
event_aligned_move = interp1(timelite.timestamps, ...
    +wheel_move,wheel_window_t_peri_event,'previous');
quiescent_trials = ~any(event_aligned_move,2);
align_times = align_times_all(quiescent_trials);
align_category = align_category_all(quiescent_trials);

% % (task stim)
% align_times = photodiode_times(1:2:end);
% align_category = ones(size(align_times));

% % (task rewards)
% align_times = reward_times;
% align_category = ones(size(align_times));

% % (sparse noise)
% px_x = 23;
% px_y = 4;
% align_times = ...
%     stim_times(find( ...
%     (noise_locations(px_y,px_x,1:end-1) == 128 & ...
%     noise_locations(px_y,px_x,2:end) == 255) | ...
%     (noise_locations(px_y,px_x,1:end-1) == 128 & ...
%     noise_locations(px_y,px_x,2:end) == 0))+1);
% align_category = ones(size(align_times));


surround_window = [-0.5,1];

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
clim(prctile(abs(aligned_px_avg(:)),100).*[-1,1]);
axis image;


%% Align mousecam to event

use_cam = mousecam_fn;
use_t = mousecam_times;

% % Get wheel movements during stim, only use quiescent trials
% framerate = 30;
% wheel_window = [0,0.5];
% wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
% wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
% event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
%     wheel_velocity,wheel_window_t_peri_event);
% wheel_thresh = 0;
% quiescent_trials = ~any(abs(event_aligned_wheel) > 0,2);
%
% use_stim = 3;
% use_align = stimOn_times(stimIDs == use_stim & quiescent_trials);
%
align_times = photodiode_times(1:2:end);
align_category = vertcat(trial_events.values.TrialStimX);

use_align = align_times(align_category == 90);

surround_frames = 30;

% Initialize video reader, get average and average difference
vr = VideoReader(use_cam);
cam_im1 = read(vr,1);

cam_align_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    surround_frames*2+1);
cam_align_diff_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    surround_frames*2);

frame_t_offset = nan(size(use_align));
for curr_align = 1:length(use_align)

    % Find closest camera frame to timepoint
    [frame_t_offset(curr_align),curr_frame] = ...
        min(abs(use_align(curr_align) - use_t));

    % Pull surrounding frames
    curr_surround_frames = curr_frame + [-surround_frames,surround_frames];
    curr_clip = double(squeeze(read(vr,curr_surround_frames)));
    curr_clip_diff = abs(diff(curr_clip,[],3));

    cam_align_avg = cam_align_avg + curr_clip./length(use_align);
    cam_align_diff_avg = cam_align_diff_avg + curr_clip_diff./length(use_align);

    AP_print_progress_fraction(curr_align,length(use_align));
end

surround_t = [-surround_frames:surround_frames]./vr.FrameRate;
AP_imscroll(cam_align_avg,surround_t)
axis image;

surround_t = [-surround_frames:surround_frames]./vr.FrameRate;
AP_imscroll(cam_align_diff_avg,surround_t(2:end))
axis image;

% Plot difference within window
use_t = [0,0.2];
use_t_idx = surround_t >= use_t(1) & surround_t <= use_t(2);
figure;
imagesc(nanmean(cam_align_diff_avg(:,:,use_t_idx(2:end)),3));
axis image off;

%% Align wheel to event

align_times = photodiode_times(1:2:end);

surround_time = [-10,10];
surround_sample_rate = 100;
surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);
pull_times = align_times + surround_time_points;

[wheel_velocity,wheel_move] = ...
    AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);

event_aligned_wheel_vel = interp1(timelite.timestamps, ...
    wheel_velocity,pull_times);
event_aligned_wheel_move = interp1(timelite.timestamps, ...
    +wheel_move,pull_times,'previous');


figure;
subplot(2,2,1);
imagesc(surround_time_points,[],event_aligned_wheel_vel)
xline(0,'color','r');
clim(max(abs(clim)).*[-1,1])
colormap(gca,AP_colormap('BWR'));

subplot(2,2,3);
plot(surround_time_points,nanmean(event_aligned_wheel_vel,1));
xline(0,'color','r');

subplot(2,2,2);
imagesc(surround_time_points,[],event_aligned_wheel_move)
xline(0,'color','r');
ylabel('Velocity');
xlabel('Time from event');

subplot(2,2,4);
plot(surround_time_points,nanmean(event_aligned_wheel_move,1));
xline(0,'color','r');
ylabel('Move prob.');
xlabel('Time from event');


[~,sort_idx] = sort([trial_events.values.TrialQuiescence]);
figure;
imagesc(surround_time_points,[],event_aligned_wheel_move(sort_idx,:));
xline(0,'color','r','linewidth',2);
hold on;
plot(-[trial_events.values(sort_idx).TrialQuiescence],1:length(trial_events.values),'b','linewidth',2);


%% Testing MCMS API

% MCMS API documentation is here:
% https://oxford-uat.colonymanagement.org/api/swagger-ui/index.html

% Get authentication token

basicUrl = 'https://oxford-uat.colonymanagement.org/api';
authenticateEndpoint = [basicUrl '/authenticate'];

usr = 'ap7';
psw = 'Bluecookiejar';

headers = struct;
headers.Accept = '*/*';
headers.username = usr;
headers.password = psw;

header_cell = [fieldnames(headers),struct2cell(headers)];

options = weboptions( ...
    'MediaType','application/json', ...
    'ContentType','json', ...
    'RequestMethod','post', ...
    'HeaderFields',header_cell);
mcms_token = webread(authenticateEndpoint,options);

% Get procedure list

proceduresEndpoint = [basicUrl '/procedures'];
headers = struct;
headers.Accept = 'application/json';
headers.Authorization = ['Bearer ' mcms_token.token];

header_cell = [fieldnames(headers),struct2cell(headers)];

options = weboptions( ...
    'MediaType','application/json', ...
    'ContentType','json', ...
    'HeaderFields',header_cell);

data = webread(proceduresEndpoint,options);

% Get weights

curr_animal = '02150140';

endpoint = [basicUrl '/animalweights/animal/' curr_animal];
headers = struct;
headers.Accept = 'application/json';
headers.Authorization = ['Bearer ' mcms_token.token];

header_cell = [fieldnames(headers),struct2cell(headers)];

options = weboptions( ...
    'MediaType','application/json', ...
    'ContentType','json', ...
    'HeaderFields',header_cell);

data = webread(endpoint,options);


data_timestamps = datetime({data.sampleDate},'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSZ','TimeZone','local');

[~,sort_idx] = sort(data_timestamps);
[data(sort_idx).weightValue]


% Get name

curr_animal = 'TOAA1.1a';
endpoint = [basicUrl '/animals/name/' curr_animal];
headers = struct;
headers.Accept = 'application/json';
headers.Authorization = ['Bearer ' mcms_token.token];

header_cell = [fieldnames(headers),struct2cell(headers)];

options = weboptions( ...
    'MediaType','application/json', ...
    'ContentType','json', ...
    'HeaderFields',header_cell);

data = webread(endpoint,options);


% Get project licenses

endpoint = [basicUrl '/projectlicenses'];
headers = struct;
headers.Accept = 'application/json';
headers.Authorization = ['Bearer ' mcms_token.token];

header_cell = [fieldnames(headers),struct2cell(headers)];

options = weboptions( ...
    'MediaType','application/json', ...
    'ContentType','json', ...
    'HeaderFields',header_cell);

data = webread(endpoint,options);

%% Sparse noise retinotopy
% Best SNR comes from using raw blue signal

surround_window = [0.3,0.5]; % 6s = [0.3,0.5], deconv = [0.05,0.15]
framerate = 1./nanmean(diff(wf_times));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
response_n = nan(n_y_squares,n_x_squares);
response_grid = cell(n_y_squares,n_x_squares);
for px_y = 1:n_y_squares
    for px_x = 1:n_x_squares

        align_times = ...
            stim_times(find( ...
            (noise_locations(px_y,px_x,1:end-1) == 128 & ...
            noise_locations(px_y,px_x,2:end) == 255) | ...
            (noise_locations(px_y,px_x,1:end-1) == 128 & ...
            noise_locations(px_y,px_x,2:end) == 0))+1);

        response_n(px_y,px_x) = length(align_times);

        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < wf_times(2) | ...
            align_times + surround_time(2) > wf_times(end)) = [];

        % Get stim-aligned responses, 2 choices:

        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]);

        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = align_times + surround_time;
        frame_edges = [wf_times;wf_times(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);

        % Get stim-aligned baseline (at stim onset)
        align_baseline_times = align_times;
        align_frames_baseline = discretize(align_baseline_times,frame_edges);

        % Don't use NaN frames (delete, dirty)
        nan_stim = any(isnan(align_frames),2) | isnan(align_frames_baseline);
        align_frames(nan_stim,:) = [];
        align_frames_baseline(nan_stim,:) = [];

        % Define the peri-stim V's as subtracting first frame (baseline)
        peri_stim_v = ...
            reshape(wf_V_raw{1}(:,align_frames)',size(align_frames,1),size(align_frames,2),[]) - ...
            nanmean(reshape(wf_V_raw{1}(:,align_frames_baseline)',size(align_frames_baseline,1),size(align_frames_baseline,2),[]),2);

        mean_peri_stim_v = permute(mean(peri_stim_v,2),[3,1,2]);

        % Store V's
        response_grid{px_y,px_x} = mean_peri_stim_v;

    end
end

% Get position preference for every pixel
U_downsample_factor = 1; %2 if max method
screen_resize_scale = 1; %3 if max method
filter_sigma = (screen_resize_scale*2);

% Downsample U
[Uy,Ux,nSV] = size(wf_U);
Ud = imresize(wf_U,1/U_downsample_factor,'bilinear');

% Convert V responses to pixel responses
use_svs = 1:2000; % de-noises, otherwise size(U,3)
n_boot = 10;

response_mean_bootstrap = cellfun(@(x) bootstrp(n_boot,@mean,x')',response_grid,'uni',false);



%%%%%%% Get retinotopy (for each bootstrap)
use_method = 'com'; % max or com
vfs_boot = nan(size(Ud,1),size(Ud,2),n_boot);
for curr_boot = 1:n_boot

    response_mean = cell2mat(cellfun(@(x) x(:,curr_boot),response_mean_bootstrap(:),'uni',false)');
    stim_im_px = reshape(permute(plab.wf.svd2px(Ud(:,:,use_svs), ...
        response_mean(use_svs,:)),[3,1,2]),n_y_squares,n_x_squares,[]);
    gauss_filt = fspecial('gaussian',[n_y_squares,n_x_squares],filter_sigma);
    stim_im_smoothed = imfilter(imresize(stim_im_px,screen_resize_scale,'bilinear'),gauss_filt);

    switch use_method
        case 'max'
            % Upsample each pixel's response map and find maximum
            [~,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
            [m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
            m_yr = reshape(m_y,size(Ud,1),size(Ud,2));
            m_xr = reshape(m_x,size(Ud,1),size(Ud,2));

        case 'com'
            % Conversely, do COM on original^2
            [xx,yy] = meshgrid(1:size(stim_im_smoothed,2),1:size(stim_im_smoothed,1));
            m_xr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,xx),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));
            m_yr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,yy),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));
    end

    % Calculate and plot sign map (dot product between horz & vert gradient)

    % 1) get gradient direction
    [~,Vdir] = imgradient(imgaussfilt(m_yr,1));
    [~,Hdir] = imgradient(imgaussfilt(m_xr,1));

    % 3) get sin(difference in direction) if retinotopic, H/V should be
    % orthogonal, so the closer the orthogonal the better (and get sign)
    angle_diff = sind(Vdir-Hdir);
    angle_diff(isnan(angle_diff)) = 0;

    vfs_boot(:,:,curr_boot) = angle_diff;
end

vfs_boot_mean = imgaussfilt(nanmean(vfs_boot,3),2);

figure;
imagesc(vfs_boot_mean)
axis image off;
colormap(AP_colormap('RWB'));
title(animal)


%% Load DCIMG frames

img_fn = 'P:\Data\AP005\2023-05-11\widefield\AP005_2023-05-1100001.dcimg';

% Open file for reading
dcimg_fid = plab.wf.dcimgmex('open',img_fn);

% Load frames
load_frames = find(wf_times >= 428,1)+1:2:find(wf_times <= 438,1,'last');

im_width = plab.wf.dcimgmex( 'getparam', dcimg_fid, 'IMAGE_WIDTH' );
im_height = plab.wf.dcimgmex( 'getparam', dcimg_fid, 'IMAGE_HEIGHT' );

im = zeros(im_height,im_width,length(load_frames),'uint16');

for curr_frame_idx = 1:length(load_frames)
    im(:,:,curr_frame_idx) = ...
        plab.wf.dcimgmex('readframe',dcimg_fid,load_frames(curr_frame_idx))';
end

% Close file
plab.wf.dcimgmex('close', dcimg_fid)

% Show image
AP_imscroll(im);axis image


%% TESTING BATCH PASSIVE

animal = 'AP004';
use_workflow = 'lcr_passive';
recordings = ap.find_recordings(animal,use_workflow);

% temp: cut out bad days
% recordings(ismember({recordings.day},{'2023-05-04','2023-05-05'})) = [];
recordings(ismember({recordings.day},{'2023-05-08'})) = [];

wf_px = cell(size(recordings));

for curr_recording = 1:length(recordings)
    
    % Grab pre-load vars
    preload_vars = who;

    % Load data
    rec_day = recordings(curr_recording).day;
    rec_time = recordings(curr_recording).protocol{end};
    if ~recordings(curr_recording).widefield(end)
        continue
    end
    ap.load_experiment;

    % Align to stim and store
    % (passive)
    align_times_all = photodiode_times(1:2:end);
    align_category_all = vertcat(trial_events.values.TrialStimX);
    % (get only quiescent trials)
    [wheel_velocity,wheel_move] = AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);
    framerate = 30;
    wheel_window = [0,0.5];
    wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
    wheel_window_t_peri_event = align_times_all + wheel_window_t;
    event_aligned_move = interp1(timelite.timestamps, ...
        +wheel_move,wheel_window_t_peri_event,'previous');
    quiescent_trials = ~any(event_aligned_move,2);
    align_times = align_times_all(quiescent_trials);
    align_category = align_category_all(quiescent_trials);

    surround_window = [-0.5,1];
    surround_samplerate = 35;
    t = surround_window(1):1/surround_samplerate:surround_window(2);
    peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

    aligned_v = reshape(interp1(wf_times,wf_V',peri_event_t,'previous'), ...
        length(align_times),length(t),[]);

    align_id = findgroups(align_category);
    aligned_v_avg = permute(splitapply(@nanmean,aligned_v,align_id),[3,2,1]);
    aligned_v_avg_baselined = aligned_v_avg - nanmean(aligned_v_avg(:,t < 0,:),2);

    aligned_px_avg = plab.wf.svd2px(wf_U,aligned_v_avg_baselined);

    wf_px{curr_recording} = aligned_px_avg;

    AP_print_progress_fraction(curr_recording,length(recordings));

    % Clear vars except pre-load for next loop
    clearvars('-except',preload_vars{:});

end

surround_window = [-0.5,1];
surround_samplerate = 35;
t = surround_window(1):1/surround_samplerate:surround_window(2);
a = cellfun(@(x) max(x(:,:,t > 0 & t < 0.2,3),[],3), ...
    wf_px(cellfun(@(x) ~isempty(x),wf_px)),'uni',false);
c = [0,max(cellfun(@(x) max(x(:)),a))];

figure('Name',animal');
tiledlayout('flow')
for i = 1:length(a)
nexttile;imagesc(a{i});clim(c);axis image off;
end


%% TESTING BATCH BEHAVIOR

animal = 'AP005';
use_workflow = {'stim_wheel_right_stage1','stim_wheel_right_stage2'};
recordings = ap.find_recordings(animal,use_workflow);

surround_time = [-5,5];
surround_sample_rate = 100;
surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

n_trials_water = nan(length(recordings),2);
rxn_med = nan(length(recordings),1);
frac_move_stimalign = nan(length(recordings),length(surround_time_points));

for curr_recording = 1:length(recordings)
    
    % Grab pre-load vars
    preload_vars = who;

    % Load data
    rec_day = recordings(curr_recording).day;
    rec_time = recordings(curr_recording).protocol{end};
    load_parts.widefield = false;
    ap.load_experiment;

    % Get total trials/water
    n_trials_water(curr_recording,:) = [length(trial_events.timestamps), ...
        sum(([trial_events.values.Outcome] == 1)*6)];

    % Get median stim-outcome time
    n_trials = length([trial_events.timestamps.Outcome]);
    rxn_med(curr_recording) = median(seconds([trial_events.timestamps(1:n_trials).Outcome] - ...
    cellfun(@(x) x(1),{trial_events.timestamps(1:n_trials).StimOn})));

    % Align wheel movement to stim onset
    align_times = photodiode_times(1:2:end);
    pull_times = align_times + surround_time_points;

    [wheel_velocity,wheel_move] = ...
        AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);

    event_aligned_wheel_vel = interp1(timelite.timestamps, ...
        wheel_velocity,pull_times);
    event_aligned_wheel_move = interp1(timelite.timestamps, ...
        +wheel_move,pull_times,'previous');

    frac_move_stimalign(curr_recording,:) = nanmean(event_aligned_wheel_move,1);

    AP_print_progress_fraction(curr_recording,length(recordings));

    % Clear vars except pre-load for next loop
    clearvars('-except',preload_vars{:});

end

relative_day = days(datetime({recordings.day}) - datetime({recordings(1).day}))+1;
nonrecorded_day = setdiff(1:length(recordings),relative_day);

figure('Name',animal);
tiledlayout('flow');

nexttile;
yyaxis left; plot(relative_day,n_trials_water(:,1));
ylabel('# trials');
yyaxis right; plot(relative_day,n_trials_water(:,2));
ylabel('Reward amount');
xlabel('Day');
if any(nonrecorded_day)
    xline(nonrecorded_day,'--k');
end

nexttile;
yyaxis left
plot(relative_day,rxn_med)
set(gca,'YScale','log');
ylabel('Med. rxn');
xlabel('Day');
if any(nonrecorded_day)
    xline(nonrecorded_day,'--k');
end

yyaxis right
prestim_max = max(frac_move_stimalign(:,surround_time_points < 0),[],2);
poststim_max = max(frac_move_stimalign(:,surround_time_points > 0),[],2);
plot(relative_day,(poststim_max-prestim_max)./(poststim_max+prestim_max));
yline(0);
ylabel('pre/post move idx');
xlabel('Day');

nexttile;
imagesc(surround_time_points,[],frac_move_stimalign);
clim([0,1]);
colormap(gca,AP_colormap('WK'));
set(gca,'YTick',1:length(recordings),'YTickLabel',string(datetime({recordings.day},'format','MM-dd')));
xlabel('Time from stim');

nexttile; hold on
set(gca,'ColorOrder',copper(length(recordings)));
plot(surround_time_points,frac_move_stimalign');
xline(0,'color','k');
ylabel('Fraction moving');
xlabel('Time from stim');







