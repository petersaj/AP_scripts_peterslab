%% SANDBOX
% 
% Temporary code

%% Align widefield to event

align_times = photodiode_times(1:2:end);
align_category = vertcat(trial_events.values.TrialStimX);

% align_times = reward_times;
% align_category = ones(size(align_times));

% (get only quiescent trials)
[wheel_velocity,wheel_move] = AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);
framerate = 30;
wheel_window = [0,0.5];
wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
wheel_window_t_peri_event = align_times + wheel_window_t;
event_aligned_move = interp1(timelite.timestamps, ...
    +wheel_move,pull_times,'previous');
quiescent_trials = ~any(event_aligned_move,2);

align_times = align_times(quiescent_trials);
align_category = align_category(quiescent_trials);


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
















