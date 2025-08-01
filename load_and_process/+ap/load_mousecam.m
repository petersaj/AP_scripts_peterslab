% Load mousecam
%
% Get timelite times for each corresponding mousecam frame
%
% NOTE: possibly not robust if there is a flip mid-exposure on the first
% recorded frame, which would cause all frames to be pushed forward or back
% 1 frames.
%
% (updated after AP022 2024-05-20 can see reflection of stim in headplate
% holder, which gives timing validation)

if verbose; disp('Loading Mousecam...'); end

mousecam_fn = plab.locations.filename('server',animal,rec_day,rec_time,'mousecam','mousecam.mj2');
mousecam_header_fn = plab.locations.filename('server',animal,rec_day,rec_time,'mousecam','mousecam_header.bin');

if ~exist(mousecam_fn,'file')
    warning('%s %s: no mousecam found',animal,rec_day);
    return
end

% Read mousecam header and get flipper times
mousecam_flipper_pin = 2;
mousecam_header = plab.mousecam.read_mousecam_header(mousecam_header_fn, mousecam_flipper_pin);
mousecam_flipper_times = mousecam_header.timestamps(find(diff(mousecam_header.flipper) ~= 0) + 1);

% Count mousecam exposures per flip (mousecam and timelite)
mousecam_exposures_per_flip_mousecam = diff(find(diff(mousecam_header.flipper) ~= 0));
mousecam_exposures_per_flip_timelite = histcounts(mousecam_exposeOn_times,flipper_times)';

% Find flip epoch delay between mousecam and timelite
% (best match between number of exposures in flip epochs)
mousecam_flip_epoch_delay = ...
    finddelay(mousecam_exposures_per_flip_mousecam,mousecam_exposures_per_flip_timelite);

% Get flip epochs with matching number of frames
mousecam_num_epochs = min(length(mousecam_exposures_per_flip_mousecam), ...
    length(mousecam_exposures_per_flip_mousecam)-mousecam_flip_epoch_delay);
mousecam_match_epochs = find(mousecam_exposures_per_flip_mousecam(1:mousecam_num_epochs) == ...
    mousecam_exposures_per_flip_timelite(mousecam_flip_epoch_delay+1: ...
    mousecam_flip_epoch_delay+mousecam_num_epochs));

% Get first frame index for each matching epoch in timelite/mousecam
mousecam_match_epoch_startframe_mousecam = ...
    interp1(mousecam_header.timestamps,1:length(mousecam_header.timestamps), ...
    mousecam_flipper_times(mousecam_match_epochs),'next');
mousecam_match_epoch_startframe_timelite = ...
    interp1(mousecam_exposeOn_times,1:length(mousecam_exposeOn_times), ...
    flipper_times(mousecam_flip_epoch_delay+mousecam_match_epochs),'next');

% Define the first recorded frame as the most common index difference
% (not always the same because sometimes flip happened mid-exposure, so
% using the most common difference gives this leeway)
mousecam_first_frame_idx = mode( ...
    mousecam_match_epoch_startframe_timelite - ...
    mousecam_match_epoch_startframe_mousecam)+1;

% Get mousecam frame times on timelite clock
% (extrapolate pre/post-timelite times)
mousecam_frame_timelite_idx = mousecam_first_frame_idx: ...
    mousecam_first_frame_idx+length(mousecam_header.timestamps)-1;
mousecam_times = interp1(1:length(mousecam_exposeOn_times),mousecam_exposeOn_times, ...
    mousecam_frame_timelite_idx,'linear','extrap');








