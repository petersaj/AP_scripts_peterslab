% Load widefield data

if verbose; disp('Loading Widefield...'); end

%% Set locations

% Load widefield data for all colors
widefield_colors = {'blue','violet'};

wf_day_path = plab.locations.filename('server',animal,rec_day,[],'widefield');
wf_rec_path = plab.locations.filename('server',animal,rec_day,rec_time,'widefield');

%% Load data

[wf_avg_all,wf_U_raw,wf_V_raw] = deal(cell(length(widefield_colors),1));
for curr_wf = 1:length(widefield_colors)
    mean_image_fn = fullfile(wf_day_path,sprintf('meanImage_%s.npy',widefield_colors{curr_wf}));
    svdU_fn = fullfile(wf_day_path,sprintf('svdSpatialComponents_%s.npy',widefield_colors{curr_wf}));
    svdV_fn = fullfile(wf_rec_path,sprintf('svdTemporalComponents_%s.npy',widefield_colors{curr_wf}));

    wf_avg_all{curr_wf} = readNPY(mean_image_fn);
    wf_U_raw{curr_wf} = readNPY(svdU_fn);
    wf_V_raw{curr_wf} = readNPY(svdV_fn);
end

%% Get frame times by color
% (check for timelite/frame mismatches and dropped frames)

% Get widefield colors (assume alternating)
widefield_frame_colors = mod((1:length(widefield_expose_times))'-1,length(widefield_colors))+1;

wf_cam_tl_frame_diff = sum(cellfun(@(x) size(x,2),wf_V_raw)) - ...
    length(widefield_expose_times);

if wf_cam_tl_frame_diff == 0
    % Same frames as expose times: all good
    wf_use_frames = true(size(widefield_expose_times));

elseif wf_cam_tl_frame_diff > 0
    % More frames than timelite exposures: timelite probably cut off early 
    warning('Widefield: timelite missed [%d] exposures, cutting from end',wf_cam_tl_frame_diff)

    % (get expected frames by color, truncate V
    n_color_frames = accumarray(widefield_frame_colors,1);
    wf_V_raw = cellfun(@(v,n_frames) v(:,1:n_frames),wf_V_raw,num2cell(n_color_frames),'uni',false);
    wf_use_frames = true(size(widefield_expose_times));

elseif wf_cam_tl_frame_diff < 0
    % More timelite exposures than frames: widefield dropped frames
    warning('Widefield: dropped [%d] frame(s), removing dropped exposure times',-wf_cam_tl_frame_diff);

    widefield_metadata_fn = ...
        plab.locations.filename('server',animal,rec_day,[], ...
        'widefield',sprintf('widefield_%s_metadata.bin',rec_time));

    [wf_dropped_frames,wf_dropped_frame_idx] = ...
        plab.wf.find_dropped_frames(widefield_metadata_fn,false);

    if -wf_cam_tl_frame_diff == length(wf_dropped_frames)
        % If dropped frames were all identified, use non-dropped frames
        wf_use_frames = ~wf_dropped_frame_idx;

    elseif (-wf_cam_tl_frame_diff - length(wf_dropped_frames)) < 300
        if wf_cam_tl_frame_diff < -5 && isempty(wf_dropped_frames)
            % RARE CASE: if few missing frames and no dropped frames
            % detected, assume it was the last frame that was dropped (since
            % this drop is undetectable. Discovered in AP029 2024-12-09 1356)
            wf_use_frames = [true(length(widefield_expose_times)-1,1);false];
            warning('Widefield: dropped [%d] frames that could not be detected, assuming drops at the end',-wf_cam_tl_frame_diff);
        else
            % If few dropped frames and not all identified, error out
            error('Widefield: could not identify dropped frames')
        end

    elseif (-wf_cam_tl_frame_diff - length(wf_dropped_frames)) > 300
        % If large number of dropped frames, assume the remainder at end
        % (note: this is very messy - made to handle cases where the camera
        % broke at the end of the recording)
        wf_use_frames = false(size(widefield_expose_times));
        wf_use_frames(1:length(wf_dropped_frame_idx)) = ~wf_dropped_frame_idx;
        warning('Widefield: %d dropped frames could not be identified, assuming at end of recording', ...
            -wf_cam_tl_frame_diff - length(wf_dropped_frames));

    end

end

% Get timestamps for widefield frames by color
wf_t_all = arrayfun(@(x) ...
    widefield_expose_times(wf_use_frames & widefield_frame_colors == x), ...
    (1:length(widefield_colors))','uni',false);

% Get framerate from timestamps
wf_framerate = mean(1./diff(wf_t_all{1}));


%% Correct hemodynamics

[V_neuro_hemocorr,hemocorr_t] = plab.wf.hemo_correct( ...
    wf_U_raw{1},wf_V_raw{1},wf_t_all{1}, ...
    wf_U_raw{2},wf_V_raw{2},wf_t_all{2});


%% Normalize to DF/F

wf_Vdf = plab.wf.svd_dff(wf_U_raw{1},V_neuro_hemocorr,wf_avg_all{1});

%% Deconvolve

wf_Vdf_deconv = ap.wf_deconv(wf_Vdf,wf_framerate);

%% Set final processed widefield variables

wf_U = wf_U_raw{1};
wf_V = wf_Vdf_deconv;
wf_t = hemocorr_t;
wf_avg = wf_avg_all{1};

%% Align widefield (if load_parts.widefield_align not turned off)

% load_parts.widefield_align: apply day/animal alignments to widefield
% load_parts.widefield_master: convert widefield V's into master U's

if ~isfield(load_parts,'widefield_align') || ...
        load_parts.widefield_align
    try
        wf_avg = plab.wf.wf_align(wf_avg,animal,rec_day);
        wf_U = plab.wf.wf_align(wf_U,animal,rec_day);
        if verbose; disp('Aligned widefield U/avg...');end
        load_parts.widefield_align = true;
        
        if isfield(load_parts,'widefield_master') && load_parts.widefield_master
            % Convert basis set to U master
            [U_master,V_master] = plab.wf.u2master(wf_U,wf_V);
            % Set U/V to master
            wf_U = U_master;
            wf_V = V_master;
            % Clear temporary conversions
            clear U_master V_master
            if verbose; disp('Changed widefield basis to U master...');end
        end

    catch  me
        warning(me.identifier,'Widefield: %s',me.message);
        load_parts.widefield_align = false;
        load_parts.widefield_master = false;
    end
end










