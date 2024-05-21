% Load mousecam
%
% Get timelite times for each corresponding mousecam frame

if verbose; disp('Loading Mousecam...'); end

mousecam_fn = plab.locations.filename('server',animal,rec_day,rec_time,'mousecam','mousecam.mj2');
mousecam_header_fn = plab.locations.filename('server',animal,rec_day,rec_time,'mousecam','mousecam_header.bin');

if load_parts.mousecam && exist(mousecam_fn,'file')

    % Read mousecam header and get flipper times
    mousecam_flipper_pin = 2;
    mousecam_header = plab.mousecam.read_mousecam_header(mousecam_header_fn, mousecam_flipper_pin);
    mousecam_flipper_times = mousecam_header.timestamps(find(diff(mousecam_header.flipper) ~= 0) + 1);

    % Count mousecam exposures per flip (mousecam and timelite)
    mousecam_exposures_per_flip_mousecam = diff(unique([0;length(mousecam_header.flipper); ...
        find(diff(mousecam_header.flipper) ~= 0)]));

    mousecam_exposures_per_flip_timelite = ...
        accumarray(discretize(mousecam_expose_times, ...
        [0;flipper_times;timelite.timestamps(end)]),1);

    % Find flip epoch delay between mousecam and timelite
    % (best match between number of exposures in flip epochs)
    mousecam_flip_epoch_delay = ...
        finddelay(mousecam_exposures_per_flip_mousecam,mousecam_exposures_per_flip_timelite);

    % Get indicies, then expose times, for mousecam frames
    mousecam_first_frame_idx = ...
        (1+sum(mousecam_exposures_per_flip_timelite(1:mousecam_flip_epoch_delay+1)) - mousecam_exposures_per_flip_mousecam(1));

    mousecam_last_frame_idx = length(mousecam_expose_times) - ...
        (mousecam_exposures_per_flip_timelite(mousecam_flip_epoch_delay+length(mousecam_exposures_per_flip_mousecam)) - ...
        mousecam_exposures_per_flip_mousecam(end));

    mousecam_times = mousecam_expose_times(mousecam_first_frame_idx:mousecam_last_frame_idx);

    if length(mousecam_times) ~= length(mousecam_header.timestamps)
        error('Mousecam: frame time mismatch, incorrectly aligned to timelite')
    end

end









