
% animal = 'AP006';
% 
% % use_workflow = 'lcr_passive';
% use_workflow = {'stim_wheel_right_stage1','stim_wheel_right_stage2'};
% % use_workflow = 'sparse_noise';
% 
% recordings = ap.find_workflow(animal,use_workflow);
% 
% % use_rec = 5;
% 
% rec_day = '2023-05-23';
% use_rec = strcmp(rec_day,{recordings.day});
% 
% rec_day = recordings(use_rec).day;
% % rec_time = recordings(use_rec).protocol{1}; % (if multiple, use first)
% rec_time = recordings(use_rec).protocol{end}; % (if multiple, use last)

%% Define parts to load

%% Define what to load

% If nothing specified, load everything (but not LFP)
if ~exist('load_parts','var')
    load_parts.mousecam = true;
    load_parts.widefield = true;
    load_parts.ephys = true;
else
    % If only some things specified, don't load others
    if ~isfield(load_parts,'mousecam')
        load_parts.mousecam = false;
    end
    if ~isfield(load_parts,'widefield')
        load_parts.widefield = false;
    end
    if ~isfield(load_parts,'ephys')
        load_parts.ephys = false;
    end
end


%% Load timelite and associated inputs

% Set level for TTL threshold
ttl_thresh = 2;

% Load timelite
timelite_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'timelite.mat');
timelite = load(timelite_fn);

% % Plot data
% figure;
% stackedplot(timelite.timestamps,timelite.data,'DisplayLabels',{timelite.daq_info.channel_name})

% Flipper times
flipper_idx = strcmp({timelite.daq_info.channel_name}, 'flipper');
flipper_thresh = timelite.data(:,flipper_idx) >= ttl_thresh;
flipper_times = timelite.timestamps(find(diff(flipper_thresh) ~= 0) + 1);

% Mousecam times (note: all exposures, including those not recorded)
mousecam_idx = strcmp({timelite.daq_info.channel_name}, 'mouse_camera');
mousecam_thresh = timelite.data(:,mousecam_idx) >= ttl_thresh;
mousecam_expose_times = timelite.timestamps(find(diff(mousecam_thresh) == 1) + 1);

% Widefield times
widefield_idx = strcmp({timelite.daq_info.channel_name}, 'widefield_camera');
widefield_thresh = timelite.data(:,widefield_idx) >= ttl_thresh;
widefield_expose_times = timelite.timestamps(find(diff(widefield_thresh) == 1) + 1);

% Wheel position and velocity
timelite_wheel_idx = strcmp({timelite.daq_info.channel_name}, 'wheel');
wheel_position = timelite.data(:,timelite_wheel_idx);
[wheel_velocity,wheel_move] = AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);

% Screen on times
screen_idx = strcmp({timelite.daq_info.channel_name}, 'stim_screen');
screen_on = timelite.data(:,screen_idx) > ttl_thresh;

% Photodiode flips (interpolate from previous across screen flicker)
photodiode_thresh_level = 1; % low: has a relatively slow rise time
photodiode_idx = strcmp({timelite.daq_info.channel_name}, 'photodiode');
photodiode_thresh_screen_on = medfilt1(timelite.data(screen_on,photodiode_idx),3);
photodiode_thresh = interp1(timelite.timestamps(screen_on),photodiode_thresh_screen_on, ...
    timelite.timestamps,'previous','extrap') > photodiode_thresh_level;
photodiode_times = timelite.timestamps(find(diff(photodiode_thresh) ~= 0) + 1);

% Reward times (if on past certain time, valve signal flips rapidly to
% avoid burnout - take the reward onset as flipping up and staying high for
% some length of samples)
reward_idx = strcmp({timelite.daq_info.channel_name}, 'reward_valve');
reward_thresh = timelite.data(:,reward_idx) >= ttl_thresh;
reward_on_pattern = [0,ones(1,3)]; % flip up, be consecutively high
reward_times = timelite.timestamps(strfind(reward_thresh',reward_on_pattern));

%% Load Bonsai

% Bonsai events (put in all my protocols)
try
    bonsai_file = 'bonsai_events.csv';
    bonsai_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'bonsai',bonsai_file);
    trial_events = AP_load_bonsai(bonsai_fn);
catch me
end

% Sparse noise: get noise locations and times
try
    bonsai_file = 'NoiseLocations.bin';
    bonsai_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'bonsai',bonsai_file);
    fid = fopen(bonsai_fn);

    n_x_squares = trial_events.parameters.ScreenExtentX./trial_events.parameters.StimSize;
    n_y_squares = trial_events.parameters.ScreenExtentY./trial_events.parameters.StimSize;

    noise_locations = reshape(fread(fid),n_y_squares,n_x_squares,[]);
    fclose(fid);

    % Get stim times from photodiode (extrapolate: sparse noise photodiode
    % flips every N stim to give a more robust signal)
    photodiode_stim_idx = 1:trial_events.parameters.NthPhotodiodeFlip:size(noise_locations,3);
    % (check that the number of photodiode flips is expected)
    if length(photodiode_stim_idx) ~= length(photodiode_times)
        error('Sparse noise: mismatch photodiode times and stim number')
    end

    stim_times = interp1(photodiode_stim_idx,photodiode_times, ...
        1:size(noise_locations,3),'linear','extrap')';

catch me
end


%%%%% testing: bonsai times into timelite times
% bonsai_reward_datetime = [trial_events.timestamps([trial_events.values.Outcome] == 1).Outcome];
% bonsai_relative_t = bonsai_reward_datetime(1);
% 
% bonsai_reward_t = seconds(bonsai_reward_datetime - bonsai_relative_t);
% 
% % event_datetime = cellfun(@(x) x(1),{trial_events.timestamps.StimOn}');
% event_datetime = vertcat(trial_events.timestamps.QuiescenceStart);
% % event_datetime = vertcat(trial_events.timestamps.QuiescenceReset);
% 
% event_t_bonsai = seconds(event_datetime - bonsai_relative_t);
% 
% event_t_tl = interp1(bonsai_reward_t,reward_times,event_t_bonsai,'linear','extrap');



%% Load mousecam

if load_parts.mousecam

mousecam_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'mousecam','mousecam.mj2');
mousecam_header_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'mousecam','mousecam_header.bin');

% Read mousecam header and get flipper times
mousecam_flipper_pin = 2;
mousecam_header = plab.mousecam.read_mousecam_header(mousecam_header_fn, mousecam_flipper_pin);
mousecam_flipper_times = mousecam_header.timestamps(find(diff(mousecam_header.flipper) ~= 0) + 1);

% Check that timelite and mousecam have equal flipper flips
if length(flipper_times) ~= length(mousecam_flipper_times)
    warning('Flipper times not matched in timelite and mousecam');
    min_flipper_n = min(length(flipper_times),length(mousecam_flipper_times));
    %%% temporary?
    flipper_times = flipper_times(1:min_flipper_n);
    mousecam_flipper_times = mousecam_flipper_times(1:min_flipper_n);
end

% Get frame time after flips in timeline and mousecam
mousecam_postflips_idx_tl = arrayfun(@(x) ...
    find(mousecam_expose_times > flipper_times(x),1), ...
    1:length(flipper_times))';

mousecam_postflips_idx_cam = find(diff(mousecam_header.flipper) ~= 0) + 1;
% %%% temporary?
% mousecam_postflips_idx_cam = mousecam_postflips_idx_cam(1:min_flipper_n);

% For sync: only use frames where flip happened in window before frame
% started (if flip happens during/close to exposure - camera pin state can
% be ambiguous)
mousecam_use_flips = ...
    (mousecam_expose_times(mousecam_postflips_idx_tl) - flipper_times) > 0.005 & ...
    (mousecam_expose_times(mousecam_postflips_idx_tl) - flipper_times) < 0.02;

use_flipframes = setdiff(1:length(flipper_times), ...
    find(diff(mousecam_postflips_idx_tl) ~= diff(mousecam_postflips_idx_cam)) + [0,1]);

% Get offset between frame index in timelite and mousecam
mousecam_idx_offset = unique( ...
    mousecam_postflips_idx_tl(mousecam_use_flips) - ...
    mousecam_postflips_idx_cam(mousecam_use_flips));

% If there's more than one offset value, something's misaligned
if length(mousecam_idx_offset) ~= 1
    error('Mousecam frames misaligned: >1 offset value')
end

% Get the corresponding timelite frame times for each mousecam frame
mousecam_tl_idx = (1:length(mousecam_header.timestamps)) + mousecam_idx_offset;
mousecam_tl_captured = mousecam_tl_idx > 0 & mousecam_tl_idx <= length(mousecam_expose_times);
mousecam_times = mousecam_expose_times(mousecam_tl_idx(mousecam_tl_captured));

end

%% Load widefield

if load_parts.widefield && ...
        exist(plab.locations.make_server_filename(animal,rec_day,[],'widefield'),'dir')

% Load widefield data for all colors
widefield_colors = {'blue','violet'};
[wf_avg_all,wf_U_raw,wf_V_raw,wf_t_all] = deal(cell(length(widefield_colors),1));
for curr_wf = 1:length(widefield_colors)
    mean_image_fn = plab.locations.make_server_filename(animal,rec_day,[], ...
        'widefield',sprintf('meanImage_%s.npy',widefield_colors{curr_wf}));
    svdU_fn = plab.locations.make_server_filename(animal,rec_day,[], ...
        'widefield',sprintf('svdSpatialComponents_%s.npy',widefield_colors{curr_wf}));
    svdV_fn = plab.locations.make_server_filename(animal,rec_day,rec_time, ...
        'widefield',sprintf('svdTemporalComponents_%s.npy',widefield_colors{curr_wf}));

    wf_avg_all{curr_wf} = readNPY(mean_image_fn);
    wf_U_raw{curr_wf} = readNPY(svdU_fn);
    wf_V_raw{curr_wf} = readNPY(svdV_fn);

    % Timestamps: assume colors go in order (dictated by Arduino)
    wf_t_all{curr_wf} = widefield_expose_times(curr_wf:length(widefield_colors):end);
end

% Correct hemodynamics
V_neuro_hemocorr = plab.wf.hemo_correct( ...
    wf_U_raw{1},wf_V_raw{1},wf_t_all{1}, ...
    wf_U_raw{2},wf_V_raw{2},wf_t_all{2});

% Get DF/F
wf_Vdf = plab.wf.svd_dff(wf_U_raw{1},V_neuro_hemocorr,wf_avg_all{1});

% Deconvolve
wf_framerate = mean(1./diff(wf_t_all{1}));
wf_Vdf_deconv = AP_deconv_wf(wf_Vdf,[],wf_framerate);

% Set final processed widefield variables
wf_U = wf_U_raw{1};
wf_V = wf_Vdf_deconv;
wf_times = wf_t_all{1};
wf_avg = wf_avg_all{1};

end 

%% Experiment scroller


% AP_expscroll(wf_U_raw{1},AP_deconv_wf(wf_V_raw{1},[],wf_framerate),wf_t_all{1},mousecam_fn,mousecam_times)

% AP_expscroll(wf_U_raw{1},wf_V_raw{1},wf_t_all{1},mousecam_fn,mousecam_times)
% AP_expscroll(wf_U_raw{2},wf_V_raw{2},wf_t_all{2},mousecam_fn,mousecam_times)

% AP_expscroll(wf_U,wf_Vdf,wf_times,mousecam_fn,mousecam_times)

% AP_expscroll(wf_U,wf_V,wf_times,mousecam_fn,mousecam_times)













