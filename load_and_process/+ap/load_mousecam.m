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
        accumarray(discretize(mousecam_exposeOn_times, ...
        [0;flipper_times;timelite.timestamps(end)]),1);

    % Find flip epoch delay between mousecam and timelite
    % (best match between number of exposures in flip epochs)
    mousecam_flip_epoch_delay = ...
        finddelay(mousecam_exposures_per_flip_mousecam,mousecam_exposures_per_flip_timelite);

    % Get index for first recorded frame
    mousecam_first_frame_idx = ...
        (1+sum(mousecam_exposures_per_flip_timelite(1:mousecam_flip_epoch_delay+1)) - mousecam_exposures_per_flip_mousecam(1));

    % Get mousecam times from first recorded frame to mousecam N frames
    mousecam_times = mousecam_exposeOn_times(mousecam_first_frame_idx: ...
        mousecam_first_frame_idx+length(mousecam_header.timestamps)-1);

end

% ~~~~~~~~~~~~~~~
% Possibly more robust method, but would take more work: get "matching"
% frames between mousecam and timelite, interpolate to full frame index,
% ensure all integer values after interpolating. Doesn't currently work -
% likely because "matching" epochs include by chance many times where there
% was a flip mid-exposure, so these aren't actually "matching" frames:

% % Get flip epochs with matching number of frames
% mousecam_match_epochs = find(mousecam_exposures_per_flip_mousecam == ...
%     mousecam_exposures_per_flip_timelite(mousecam_flip_epoch_delay+1:end));
% 
% % Find corresponding frames in timelite/mousecam
% mousecam_flip_epochs_mousecam = discretize(mousecam_header.timestamps,[0;mousecam_flipper_times;mousecam_header.timestamps(end)]);
% mousecam_flip_epochs_timelite = discretize(mousecam_exposeOn_times,[0;flipper_times;timelite.timestamps(end)]);
% 
% mousecam_frame_match_idx_mousecam = find(ismember(mousecam_flip_epochs_mousecam,mousecam_match_epochs));
% mousecam_frame_match_idx_timelite = find(ismember(mousecam_flip_epochs_timelite,mousecam_match_epochs+mousecam_flip_epoch_delay));
% 
% % (remove corresponding mousecam exposures with mid-exposure flip)
% mousecam_midexpose_flip = find(arrayfun(@(x) ...
%     any(flipper_times >= mousecam_exposeOn_times(x) & ...
%     flipper_times <= mousecam_exposeOff_times(x)), ...
%     1:length(mousecam_exposeOn_times)));
% mousecam_match_use = ~ismember(mousecam_frame_match_idx_timelite,mousecam_midexpose_flip);
% 
% % Get timelite index for each mousecam frame
% mousecam_timelite_idx = interp1( ...
%     mousecam_frame_match_idx_mousecam(mousecam_match_use), ...
%     mousecam_frame_match_idx_timelite(mousecam_match_use), ...
%     1:length(mousecam_header.timestamps),'linear','extrap');












