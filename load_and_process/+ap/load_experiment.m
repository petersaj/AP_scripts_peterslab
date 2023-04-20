
animal = 'aptest1';
rec_day = '2023-04-19';
rec_time = '1634';

%% Load timelite and associated inputs

ttl_thresh = 2;

% Load timelite
timelite_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'timelite.mat');
timelite = load(timelite_fn);

% Plot data
figure;
stackedplot(timelite.timestamps,timelite.data,'DisplayLabels',{timelite.daq_info.channel_name})

% Flipper times
flipper_idx = strcmp({timelite.daq_info.channel_name}, 'flipper');
flipper_thresh = timelite.data(:,flipper_idx) >= ttl_thresh;
flipper_times = timelite.timestamps(find(diff(flipper_thresh) ~= 0) + 1);

% Mousecam times
mousecam_idx = strcmp({timelite.daq_info.channel_name}, 'face_camera');
mousecam_thresh = timelite.data(:,mousecam_idx) >= ttl_thresh;
mousecam_times = timelite.timestamps(find(diff(mousecam_thresh) == 1) + 1);

% Widefield times
widefield_idx = strcmp({timelite.daq_info.channel_name}, 'widefield_camera');
widefield_thresh = timelite.data(:,widefield_idx) >= ttl_thresh;
widefield_times = timelite.timestamps(find(diff(widefield_thresh) == 1) + 1);

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

bonsai_file = 'test.csv';
bonsai_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,bonsai_file);
trial_events = AP_load_bonsai(bonsai_fn);


%% Load mousecam

mousecam_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'mousecam.mj2');
mousecam_header_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'mousecam_header.bin');

% Read mousecam header and get flipper times
mousecam_flipper_pin = 2;
mousecam_header = read_mousecam_header(mousecam_header_fn, mousecam_flipper_pin);
mousecam_flipper_times = mousecam_header.timestamps(find(diff(mousecam_header.flipper) ~= 0) + 1);

% Check that timelite and mousecam have equal flipper flips
if length(flipper_times) ~= length(mousecam_flipper_times)
    error('Flipper times not matched in timelite and mousecam');
end

%%%%% WORKING HERE: syncing mousecam is a lot harder than expected

% Convert timestamps from mousecam to timelite: 
% - use flips that didn't happen during an exposure

mousecam_times_prepostflips = cell2mat(arrayfun(@(x) ...
    mousecam_times([ ...
    find(mousecam_times <= flipper_times(x),1,'last'), ...
    find(mousecam_times >= flipper_times(x),1,'first')]), ...
    1:length(flipper_times),'uni',false))';

mousecam_use_flips = ...
    min(abs(flipper_times - mousecam_times_prepostflips),[],2) > 0.01;


% - get timelite timestamps for mousecam after each flip
mousecam_postflips_idx_tl = arrayfun(@(x) ...
    find(mousecam_times > flipper_times(x),1), ...
    1:length(flipper_times))';

mousecam_postflips_idx_cam = find(diff(mousecam_header.flipper) ~= 0) + 1;

%%%% I think this just about works:
b = round(interp1(mousecam_postflips_idx_tl(mousecam_use_flips), ...
    mousecam_postflips_idx_cam(mousecam_use_flips), ...
    1:length(mousecam_header.timestamps),'linear','extrap'));



%%%%%% things that didn't work 

b = interp1(mousecam_postflips_idx_tl, ...
    mousecam_postflips_idx_cam, ...
    1:length(mousecam_header.timestamps));

% - get mousecam timestamps for all frame with flips
mousecam_flip_frametimes = mousecam_header.timestamps(find(diff(mousecam_header.flipper) ~= 0) + 1);
% - use postflip frames as sync points to convert timestamps
mousecam_timestamps_tl = interp1( ...
    mousecam_flip_frametimes(mousecam_use_flips), ...
    mousecam_times_postflips(mousecam_use_flips), ...
    mousecam_header.timestamps,'linear');



mousecam_timestamps_tl = interp1( ...
    mousecam_flip_frametimes, ...
    flipper_times, ...
    mousecam_header.timestamps,'linear');


% only use flips that are spaced similarly between timelite and mousecam
% (won't work: we want to use ones that might have spacing, because that
% means flip happened before frame)
flip_duration_leeway = 0.01; % ms leeway difference
mousecam_use_flips = setdiff(1:length(flipper_times), ...
    find((diff(flipper_times) - diff(mousecam_flip_frametimes)) ...
    > flip_duration_leeway) + [0,1]);

figure;plot((diff(flipper_times(mousecam_use_flips)) - ...
    diff(mousecam_flip_frametimes(mousecam_use_flips))),'.')


% Use flip periods containing same number of frames 
% (sometimes don't if flip happened mid-frame)
% WON'T WORK: if happened mid-frame, can get right number but not sure
% which one it corresponds to
n_flipframes_tl = histcounts(mousecam_times,flipper_times)';
n_flipframes_mousecam = diff(find(diff(mousecam_header.flipper) ~= 0) + 1);



figure; hold on
plot(timelite.timestamps,timelite.data(:,flipper_idx))
plot(timelite.timestamps,timelite.data(:,mousecam_idx))


















































































































