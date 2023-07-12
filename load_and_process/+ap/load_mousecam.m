% Load mousecam

if verbose; disp('Loading Mousecam...'); end

mousecam_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'mousecam','mousecam.mj2');
mousecam_header_fn = plab.locations.make_server_filename(animal,rec_day,rec_time,'mousecam','mousecam_header.bin');

if load_parts.mousecam && exist(mousecam_fn,'file')

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
