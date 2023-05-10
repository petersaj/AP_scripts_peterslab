
animal = 'AP005';
rec_day = '2023-05-07';
rec_time = '1541';

%% Load timelite and associated inputs

% Set level for TTL threshold
ttl_thresh = 2;

% Load timelite
timelite_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'timelite.mat');
timelite = load(timelite_fn);

% % Plot data
% figure;
% stackedplot(timelite.timestamps,timelite.data,'DisplayLabels',{timelite.daq_info.channel_name})
% drawnow;

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

try
    bonsai_file = 'bonsai_events.csv';
    bonsai_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'bonsai',bonsai_file);
    trial_events = AP_load_bonsai(bonsai_fn);
catch me
end

%% Load mousecam

mousecam_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'mousecam','mousecam.mj2');
mousecam_header_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'mousecam','mousecam_header.bin');

% Read mousecam header and get flipper times
mousecam_flipper_pin = 2;
mousecam_header = plab.mousecam.read_mousecam_header(mousecam_header_fn, mousecam_flipper_pin);
mousecam_flipper_times = mousecam_header.timestamps(find(diff(mousecam_header.flipper) ~= 0) + 1);

% Check that timelite and mousecam have equal flipper flips
if length(flipper_times) ~= length(mousecam_flipper_times)
    error('Flipper times not matched in timelite and mousecam');
end

% Get frame time after flips in timeline and mousecam
mousecam_postflips_idx_tl = arrayfun(@(x) ...
    find(mousecam_expose_times > flipper_times(x),1), ...
    1:length(flipper_times))';

mousecam_postflips_idx_cam = find(diff(mousecam_header.flipper) ~= 0) + 1;

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
mousecam_times = mousecam_expose_times(mousecam_tl_idx);


%% Load widefield

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
    % Assume colors go in order: dictated by Arduino
    wf_t_all{curr_wf} = widefield_expose_times(curr_wf:length(widefield_colors):end);
end

wf_framerate = mean(1./diff(wf_t_all{1}));

%%%%%%%%%%%%%%%%%%%%  hemodynamic correction - put this in one function?

% Interpolate violet V values to blue timestamps
% (to estimate values at simultaneous timestamps)
wf_Vhemo_tshift = interp1(wf_t_all{2},wf_V_raw{2}',wf_t_all{1},'linear','extrap')';

% Change violet V basis set from violet U into blue U
wf_Vhemo_Uneuro = plab.wf.change_U(wf_U_raw{1},wf_V_raw{1},wf_U_raw{2});

% Get hemo tform matrix to map Vhemo onto Vneuro
hemo_freq = [5,15]; % get correction levels from heartbeat freq
hemo_local_px = 3; % pixel downsampling for local correction
skip_seconds = 20; % the beginning and end can be wonky

skip_frames = 1 + round(skip_seconds*wf_framerate);
[~,hemo_tform] = plab.wf.hemo_correct_local(wf_U_raw{1}, ...
    wf_V_raw{1}(:,skip_frames:end-skip_frames-1), ...
    wf_Vhemo_Uneuro(:,skip_frames:end-skip_frames-1), ...
    wf_framerate,hemo_freq,hemo_local_px);

% Hemo correct by subtracting scaled Vhemo from Vneuro
wf_Vneuro_hemocorrected = (wf_V_raw{1}' - (wf_Vhemo_Uneuro - nanmean(wf_Vhemo_Uneuro,2))' * hemo_tform)';

%%%%%% FILTER 

% High-pass filter (causal filter: moves forwards in time)
highpassCutoff = 0.01; % Hz
[b100s, a100s] = butter(2, highpassCutoff/(wf_framerate/2), 'high');
wf_Vneuro_hemocorrected_filtered = filter(b100s,a100s,wf_Vneuro_hemocorrected,[],2);

% Get DF/F
wf_Vdf = plab.wf.svd_dff(wf_U_raw{1}, wf_Vneuro_hemocorrected_filtered, wf_avg_all{1});

% Set final processed widefield variables
wf_U = wf_U_raw{1};
wf_V = wf_Vdf;
wf_times = wf_t_all{1};
wf_avg = wf_avg_all{1};


%%%% TEMPORARY: MOVE THIS FUNCTION 
wf_V = AP_deconv_wf(wf_Vdf,[],wf_framerate);


%% Experiment scroller (check alignment)

% AP_expscroll(wf_U_raw{2},wf_V_raw{2},wf_t_all{2},mousecam_fn,mousecam_times)


AP_expscroll(wf_U,wf_V,wf_times,mousecam_fn,mousecam_times)













