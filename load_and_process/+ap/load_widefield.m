% Load widefield data

if verbose; disp('Loading Widefield...'); end

%% Set locations

% Load widefield data for all colors
widefield_colors = {'blue','violet'};

wf_day_path = plab.locations.filename('server',animal,rec_day,[],'widefield');
wf_rec_path = plab.locations.filename('server',animal,rec_day,rec_time,'widefield');

% %%%%%%% TESTING 
% widefield_colors = {'color1','color2'};
% 
% wf_day_path = fullfile('D:\widefield_test','1color');
% wf_rec_path = fullfile('D:\widefield_test','1color',rec_time);

% wf_day_path = fullfile('D:\widefield_test','2color');
% wf_rec_path = fullfile('D:\widefield_test','2color',rec_time);

%% Load data


[wf_avg_all,wf_U_raw,wf_V_raw,wf_t_all] = deal(cell(length(widefield_colors),1));
for curr_wf = 1:length(widefield_colors)
    mean_image_fn = fullfile(wf_day_path,sprintf('meanImage_%s.npy',widefield_colors{curr_wf}));
    svdU_fn = fullfile(wf_day_path,sprintf('svdSpatialComponents_%s.npy',widefield_colors{curr_wf}));
    svdV_fn = fullfile(wf_rec_path,sprintf('svdTemporalComponents_%s.npy',widefield_colors{curr_wf}));

    wf_avg_all{curr_wf} = readNPY(mean_image_fn);
    wf_U_raw{curr_wf} = readNPY(svdU_fn);
    wf_V_raw{curr_wf} = readNPY(svdV_fn);

    % Timestamps: assume colors go in order (dictated by Arduino)
    wf_t_all{curr_wf} = widefield_expose_times(curr_wf:length(widefield_colors):end);
end


% %%%%%%%% TESTING: SINGLE SVD
% 
% 
% % Load widefield data for all colors
% widefield_colors = {''};
% 
% [wf_avg_all,wf_U_raw,wf_V_raw,wf_t_all] = deal(cell(length(widefield_colors),1));
% for curr_wf = 1:length(widefield_colors)
%     mean_image_fn = fullfile(wf_day_path,sprintf('meanImage%s.npy',widefield_colors{curr_wf}));
%     svdU_fn = fullfile(wf_day_path,sprintf('svdSpatialComponents%s.npy',widefield_colors{curr_wf}));
%     svdV_fn = fullfile(wf_rec_path,sprintf('svdTemporalComponents%s.npy',widefield_colors{curr_wf}));
% 
%     wf_avg_all{curr_wf} = readNPY(mean_image_fn);
%     wf_U_raw{curr_wf} = readNPY(svdU_fn);
%     wf_V_raw{curr_wf} = readNPY(svdV_fn);
% 
%     % Timestamps: assume colors go in order (dictated by Arduino)
%     wf_t_all{curr_wf} = widefield_expose_times(curr_wf:length(widefield_colors):end);
% end
% 
% % Correct hemodynamics
% [V_neuro_hemocorr,hemocorr_t] = plab.wf.hemo_correct( ...
%     wf_U_raw{1},wf_V_raw{1}(:,1:2:end),wf_t_all{1}(1:2:end), ...
%     wf_U_raw{1},wf_V_raw{1}(:,2:2:end),wf_t_all{1}(2:2:end));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get framerate from timestamps
wf_framerate = mean(1./diff(wf_t_all{1}));

%% Check for timelite/frame mismatches 

% Check for timelite missed exposures
if any(cellfun(@(v,t) size(v,2) > length(t),wf_V_raw,wf_t_all))
    warning('Widefield: timelite missed exposures, cutting from end')
    wf_V_raw = cellfun(@(v,t) v(:,1:length(t)),wf_V_raw,wf_t_all,'uni',false);
end

% Check for dropped frames
if any(cellfun(@(v,t) size(v,2) < length(t),wf_V_raw,wf_t_all))
    error('Widefield: dropped frames, currently unusable');
%     % At the moment: dropped frames unusable because light order switches
%     plab.wf.find_dropped_frames(animal,rec_day,rec_time)
end

%%%%%%%%%% TESTING HEMO

% U_neuro = wf_U_raw{1};
% V_neuro = wf_V_raw{1};
% t_neuro = wf_t_all{1};
% 
% U_hemo = wf_U_raw{2};
% V_hemo = wf_V_raw{2};
% t_hemo = wf_t_all{2};


U_neuro = wf_U_raw{1};
V_neuro = wf_V_raw{1}(:,1:2:end);
t_neuro = wf_t_all{1}(1:2:end);

U_hemo = wf_U_raw{1};
V_hemo = wf_V_raw{1}(:,2:2:end);
t_hemo = wf_t_all{1}(2:2:end);


%%%%%%%%%%%%%%



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
        wf_avg = ap.wf_align(wf_avg,animal,rec_day);
        wf_U = ap.wf_align(wf_U,animal,rec_day);
        if verbose; disp('Aligned widefield U/avg...');end
    catch  me
        warning(me.identifier,'Widefield: %s',me.message);
    end
end










