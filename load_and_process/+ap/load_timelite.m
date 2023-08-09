% Load timelite

if verbose; disp('Loading Timelite...'); end

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
% (if stuck high at the end, bad frame: remove)
if widefield_thresh(end) && ~isempty(widefield_expose_times)
    widefield_expose_times(end) = [];
end

% Wheel position and velocity
timelite_wheel_idx = strcmp({timelite.daq_info.channel_name}, 'wheel');
wheel_position = timelite.data(:,timelite_wheel_idx);
[wheel_velocity,wheel_move] = ap.parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);

% Screen on times
screen_idx = strcmp({timelite.daq_info.channel_name}, 'stim_screen');
screen_on = timelite.data(:,screen_idx) > ttl_thresh;

% Photodiode flips
photodiode_idx = strcmp({timelite.daq_info.channel_name}, 'photodiode');
if all(screen_on)
    % (if no screen flicker, use as is)
    photodiode_trace = timelite.data(screen_on,photodiode_idx);
else
    % (if screen flicker: median filter and interpolate across flicker)
    photodiode_trace = ...
        interp1(timelite.timestamps(screen_on), ...
        medfilt1(timelite.data(screen_on,photodiode_idx),3), ...
        timelite.timestamps,'previous','extrap');
end
% (discretize into black/NaN/white, interpolate across NaN/intermediate)
% (e.g. simulate instantaneous flips and ignore intermediate values)
photodiode_bw_thresh = [0.2,2.8]; % [black,white]
photodiode_bw = nan(size(photodiode_trace));
photodiode_bw(photodiode_trace < photodiode_bw_thresh(1)) = 0;
photodiode_bw(photodiode_trace > photodiode_bw_thresh(2)) = 1;
photodiode_bw_interp = interp1(find(~isnan(photodiode_bw)), ...
    photodiode_bw(~isnan(photodiode_bw)), ...
    1:length(photodiode_bw),'next','extrap')';

photodiode_flip_idx = find(diff(photodiode_bw_interp) ~= 0 & ...
    ~isnan(photodiode_bw_interp(2:end))) + 1;
photodiode_times = timelite.timestamps(photodiode_flip_idx);
photodiode_values = photodiode_bw_interp(photodiode_flip_idx);

% Reward times (if on past certain time, valve signal flips rapidly to
% avoid burnout - take the reward onset as flipping up and staying high for
% some length of samples)
reward_idx = strcmp({timelite.daq_info.channel_name}, 'reward_valve');
reward_thresh = timelite.data(:,reward_idx) >= ttl_thresh;
reward_on_pattern = [0,ones(1,3)]; % flip up, be consecutively high
reward_times = timelite.timestamps(strfind(reward_thresh',reward_on_pattern));

