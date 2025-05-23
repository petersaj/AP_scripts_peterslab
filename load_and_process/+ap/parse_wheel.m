function [wheel_velocity,wheel_move] = parse_wheel(wheel_position,sample_rate)
% [wheel_velocity,wheel_move,wheel_velocity_split] = parse_wheel(wheel_position,sample_rate)
% 
% Get velocity and parse movements from wheel position trace 
%
% Inputs: 
% wheel_position: assumes rotary encoder +/- 1 steps, evenly sampled
% sample_rate: sample rate for wheel position
%
% Outputs: 
% wheel_velocity: velocity of the wheel (as clicks/s)
% wheel_move: binary vector of times with movement or quiescence

% Set timing for smoothing wheel velocity
wheel_smooth_t = 0.05; % seconds

% Set number of samples and actual smoothing time
% (ensure odd number of samples to center the filter)
wheel_smooth_samples = round(wheel_smooth_t*sample_rate) + ...
    mod(round(wheel_smooth_t*sample_rate)+1,2);
wheel_smooth_samples_t = wheel_smooth_samples/sample_rate;

% Ensure wheel_position is column vector
wheel_position = reshape(wheel_position,[],1);

% Turn position into single clicks (position diff with leading 0)
% (round to ensure integer, but should be since on DAQ counter)
wheel_clicks = round([0;diff(wheel_position)]); % rotary encoder clicks

% Clean wheel clicks by removing cases where: 
% - single-timepoint wheel clicks
% - clicks that sum to zero within window
wheel_clicks_use = ...
    movsum(wheel_clicks ~= 0,wheel_smooth_samples) > 1 & ...
    movsum(wheel_clicks,wheel_smooth_samples) ~= 0;

wheel_clicks_clean = wheel_clicks.*wheel_clicks_use;

% Get velocity (from gaussian-averaged window, in clicks/s)
wheel_velocity = smoothdata(wheel_clicks_clean, ...
    'gaussian',wheel_smooth_samples)*sample_rate;

% Threshold wheel for movement, get start/stops
wheel_velocity_thresh = abs(wheel_velocity) > 0;
% (if constantly moving or not moving: set constant and return)
if length(unique(wheel_velocity_thresh)) ~= 2
    wheel_move = wheel_velocity_thresh;
    return
end

velocity_starts_all = find(diff([false;wheel_velocity_thresh;false]) == 1);
velocity_stops_all = find(diff([wheel_velocity_thresh;false]) == -1);

% Combine movements with small gaps between
combine_move_t = 0.3; % in s (empirical/arbitrary)
combine_move_samples = round(combine_move_t*sample_rate);
combine_move = find((velocity_starts_all(2:end) - velocity_stops_all(1:end-1)) > combine_move_samples);
velocity_starts_trim = velocity_starts_all([1;combine_move+1]);
velocity_stops_trim = velocity_stops_all([combine_move;end]);

% Refine wheel start/stops to first/last wheel click within movement
% (should already be from starts because window is causal)
wheel_starts = cellfun(@(curr_start,curr_stop) ...
    (curr_start-1)+find(wheel_clicks(curr_start:curr_stop)~= 0,1,'first'), ...
    num2cell(velocity_starts_trim),num2cell(velocity_stops_trim));

wheel_stops = cellfun(@(curr_start,curr_stop) ...
    (curr_start-1)+find(wheel_clicks(curr_start:curr_stop)~= 0,1,'last'), ...
    num2cell(velocity_starts_trim),num2cell(velocity_stops_trim));

% Make logical vector of movement (1) / quiescence (0)
wheel_move = logical(interp1([wheel_starts;wheel_stops;0], ...
    [ones(size(wheel_starts));zeros(size(wheel_stops));0], ...
    transpose(1:length(wheel_position)),'previous','extrap'));













