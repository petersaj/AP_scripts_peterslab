%% SANDBOX
%
% Temporary code

%% Align widefield to event

% (passive)
align_times = photodiode_times(1:2:end);
align_category = vertcat(trial_events.values.TrialStimX);

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

% % (get only quiescent trials)
% [wheel_velocity,wheel_move] = AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);
% framerate = 30;
% wheel_window = [0,0.5];
% wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
% wheel_window_t_peri_event = align_times + wheel_window_t;
% event_aligned_move = interp1(timelite.timestamps, ...
%     +wheel_move,wheel_window_t_peri_event,'previous');
% quiescent_trials = ~any(event_aligned_move,2);
% 
% align_times = align_times(quiescent_trials);
% align_category = align_category(quiescent_trials);


surround_window = [-0.5,1];

surround_samplerate = 35;
t = surround_window(1):1/surround_samplerate:surround_window(2);
peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

use_U = wf_U;
use_V = AP_deconv_wf(V_neuro_hemocorr,[],wf_framerate); %wf_V;
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
align_category = vertcat(trial_events.values.TrialStimX);

surround_time = [-0.5,2];
surround_sample_rate = 100;
surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);
pull_times = align_times + surround_time_points;

[wheel_velocity,wheel_move] = ...
    AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);

stim_aligned_wheel = interp1(timelite.timestamps, ...
    wheel_velocity,pull_times);




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


surround_window = [0.3,0.5]; % 6s = [0.3,0.5]
framerate = 1./nanmedian(diff(wf_times));
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


%% Testing hemo correction (ROI)

% IN PROGRESS: 
% trying moving mean, looks the same as moving median

[blue,roi] = AP_svd_roi(wf_U_raw{1},wf_V_raw{1},wf_avg);
violet = AP_svd_roi(wf_U_raw{2},wf_V_raw{2},[],[],roi);

blue_hc = AP_svd_roi(wf_U,wf_Vneuro_hemocorrected,[],[],roi);

skip_frames = 10;
use_frames = skip_frames:length(wf_t_all{1})-skip_frames;
% use_frames = 1:length(wf_t_all{1});
blue = blue(use_frames);
violet = violet(use_frames);
blue_hc = blue_hc(use_frames);

% Interpolate violet onto blue times (measured in alternating)
violet_bt = interp1(wf_t_all{2}(use_frames),violet,wf_t_all{1}(use_frames),'linear','extrap')';

% Get violet scale: filter in heartbeat range, regress
heartbeat_freq = [5,15];
[b,a] = butter(2,heartbeat_freq/(wf_framerate/2));
fblue = filtfilt(b,a,double(blue));
fviolet = filtfilt(b,a,double(violet_bt));

violet_scale = fviolet'\fblue';

% Baseline subtract: moving mean 
movmed_n = wf_framerate*60*2;
blue_movmed = movmean(blue,movmed_n);
violet_movmed = movmean(violet_bt,movmed_n);

blue_baselined = blue - blue_movmed;
violet_baselined = violet_bt - violet_movmed;

% Hemocorrect: subtract scaled baselined violet from baselined blue
blue_hemocorr = blue_baselined - violet_baselined*violet_scale;

% Plots
figure; tiledlayout('flow');

% Plot heartbeat-filtered traces
nexttile; hold on;
plot(fblue,'b');
plot(fviolet,'color',[0.7,0,0.7]);

% Plot traces
nexttile; hold on;
plot(blue,'b');
plot(violet_bt,'color',[0.7,0,0.7]);
plot(blue_hemocorr,'color',[0.8,0,0]);

% Plot spectrum
Fs = wf_framerate;
L = length(use_trace);
NFFT = 2^nextpow2(L);

[blue_spectrum,F] = pwelch(double(blue)',[],[],NFFT,Fs);
[violet_spectrum,F] = pwelch(double(violet_bt)',[],[],NFFT,Fs);
[blue_hemocorr_spectrum,F] = pwelch(double(blue_hemocorr)',[],[],NFFT,Fs);
[blue_hc_spectrum,F] = pwelch(double(blue_hc)',[],[],NFFT,Fs);

nexttile; hold on;
plot(F,log10(smooth(blue_spectrum,50)),'b');
plot(F,log10(smooth(violet_spectrum,50)),'color',[0.7,0,0.7]);
plot(F,log10(smooth(blue_hemocorr_spectrum,50)),'color',[0.8,0,0]);
plot(F,log10(smooth(blue_hc_spectrum,50)),'k');

xlabel('Frequency');
ylabel('Log Power');
legend({'Blue','Violet','Blue hemocorr'})


%% Testing hemo correction (whole brain)

% (I think this works fine - it should be almost the exact same as the
% other one but for some reason looks much better at getting rid of
% heartbeat, anyway the visual responses look the same)

px_downsample = 3;

% Downsample U
[Uy,Ux,nSV] = size(wf_U);
Ud_blue = imresize(wf_U_raw{1},1/px_downsample,'bilinear');
Ud_violet= imresize(wf_U_raw{2},1/px_downsample,'bilinear');

% Get all pixel traces (flat) from downsampled U
px_blue = reshape(plab.wf.svd2px(Ud_blue,wf_V_raw{1}),[],size(wf_V_raw{1},2))';
px_violet = reshape(plab.wf.svd2px(Ud_violet,wf_V_raw{2}),[],size(wf_V_raw{2},2))';

% Interpolate violet into blue timepoints (extrapolate if last unpaired)
px_violet_bt = interp1(wf_t_all{2},px_violet,wf_t_all{1},'linear','extrap');

% Filter both colors at heartbeat frequency, subtract mean
heartbeat_freq = [5,15];
[b,a] = butter(2,heartbeat_freq/(wf_framerate/2));
px_blue_heartbeat = filter(b,a,px_blue);
px_violet_bt_heartbeat = filter(b,a,px_violet_bt);

% Get scaling of violet to blue from heartbeat (don't use trace edges)
skip_frames_scale = 500;
use_frames_scale = skip_frames_scale:size(px_blue_heartbeat,1)-skip_frames_scale;
% (scaling = cov(blue-mean,viol-mean)/var(viol-mean) )
violet_scale_downsamp = sum(detrend(px_blue_heartbeat(use_frames_scale,:),'constant').*...
    detrend(px_violet_bt_heartbeat(use_frames_scale,:),'constant'))./ ...
    sum(detrend(px_violet_bt_heartbeat(use_frames_scale,:),'constant').^2);

% Get transform matrix to convert pixel scaling into V-space
T = pinv(reshape(Ud_blue,[],size(Ud_blue,3)))* ...
    diag(violet_scale_downsamp)* ...
    reshape(Ud_blue,[],size(Ud_blue,3));

% Hemo-correct blue
% (subtract scaled violet moving mean)
wf_V_violet_bt = interp1(wf_t_all{2},wf_V_raw{2}',wf_t_all{1},'linear','extrap')';

movmed_n = wf_framerate*60*2;
violet_movmean = movmean(wf_V_violet_bt,movmed_n,2);

Vhemocorr = wf_V_raw{1} - transpose((wf_V_violet_bt - violet_movmean)'*T');

Vhemocorr_df = plab.wf.svd_dff(wf_U_raw{1},Vhemocorr,wf_avg_all{1});

% Check results
AP_expscroll(wf_U_raw{1},Vhemocorr,wf_t_all{1},mousecam_fn,mousecam_times)

[a,roi] = AP_svd_roi(wf_U_raw{1},wf_V_raw{1},wf_avg);
b = AP_svd_roi(wf_U_raw{1},Vhemocorr,[],[],roi);
blue_hc = AP_svd_roi(wf_U,wf_Vneuro_hemocorrected,[],[],roi);

figure; hold on;
plot(a);plot(b);


% Plot spectrum
Fs = wf_framerate;
L = length(a);
NFFT = 2^nextpow2(L);

[blue_spectrum,F] = pwelch(double(a)',[],[],NFFT,Fs);
[blue_hemocorr_spectrum,F] = pwelch(double(b)',[],[],NFFT,Fs);
[blue_hc_spectrum,F] = pwelch(double(blue_hc)',[],[],NFFT,Fs);

nexttile; hold on;
plot(F,log10(smooth(blue_spectrum,50)),'b');
plot(F,log10(smooth(blue_hemocorr_spectrum,50)),'color',[0.8,0,0]);
plot(F,log10(smooth(blue_hc_spectrum,50)),'k');

xlabel('Frequency');
ylabel('Log Power');
legend({'Blue','Blue hemocorr','Blue old hemocorr'})


V_neuro_hemocorr = plab.wf.hemo_correct( ...
    wf_U_raw{1},wf_V_raw{1},wf_t_all{1}, ...
    wf_U_raw{2},wf_V_raw{2},wf_t_all{2});

















