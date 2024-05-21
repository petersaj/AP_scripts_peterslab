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

if sum(cellfun(@(x) size(x,2),wf_V_raw)) == length(widefield_expose_times)
    % Same frames as expose times: all good
    wf_use_frames = true(size(widefield_expose_times));

elseif sum(cellfun(@(x) size(x,2),wf_V_raw)) > length(widefield_expose_times)
    % More frames than timelite exposures: timelite probably cut off early
    warning('Widefield: timelite missed exposures, cutting from end')
    wf_V_raw = cellfun(@(v,t) v(:,1:length(t)),wf_V_raw,wf_t_all,'uni',false);
    wf_use_frames = true(size(widefield_expose_times));

elseif sum(cellfun(@(x) size(x,2),wf_V_raw)) < length(widefield_expose_times)
    % More timelite exposures than frames: widefield dropped frames
    warning('Widefield: dropped frames');

    widefield_metadata_fn = ...
        plab.locations.filename('server',animal,rec_day,[], ...
        'widefield',sprintf('widefield_%s_metadata.bin',rec_time));

    [wf_dropped_frames,wf_dropped_frame_idx] = ...
        plab.wf.find_dropped_frames(widefield_metadata_fn,false);

    if length(widefield_expose_times) - sum(cellfun(@(x) size(x,2),wf_V_raw)) ~= ...
            length(wf_dropped_frames)
        error('Widefield: could not identify dropped frames')
    end

    wf_use_frames = ~wf_dropped_frame_idx;

end


% Get widefield colors (assume alternating)
widefield_frame_colors = mod((1:length(widefield_expose_times))-1,length(widefield_colors))+1;

% Get timestamps for widefield frames by color
wf_t_all = arrayfun(@(x) ...
    widefield_expose_times(wf_use_frames & widefield_frame_colors == x), ...
    1:length(widefield_colors),'uni',false);

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
    end
end










