% Load widefield data

if verbose; disp('Loading Widefield...'); end

% Load widefield data for all colors
widefield_colors = {'blue','violet'};
[wf_avg_all,wf_U_raw,wf_V_raw,wf_t_all] = deal(cell(length(widefield_colors),1));
for curr_wf = 1:length(widefield_colors)
    mean_image_fn = plab.locations.filename('server',animal,rec_day,[], ...
        'widefield',sprintf('meanImage_%s.npy',widefield_colors{curr_wf}));
    svdU_fn = plab.locations.filename('server',animal,rec_day,[], ...
        'widefield',sprintf('svdSpatialComponents_%s.npy',widefield_colors{curr_wf}));
    svdV_fn = plab.locations.filename('server',animal,rec_day,rec_time, ...
        'widefield',sprintf('svdTemporalComponents_%s.npy',widefield_colors{curr_wf}));

    wf_avg_all{curr_wf} = readNPY(mean_image_fn);
    wf_U_raw{curr_wf} = readNPY(svdU_fn);
    wf_V_raw{curr_wf} = readNPY(svdV_fn);

    % Timestamps: assume colors go in order (dictated by Arduino)
    wf_t_all{curr_wf} = widefield_expose_times(curr_wf:length(widefield_colors):end);
end

% % BUG? FOR NOW: if mismatching number of frames/times, cut off ends to match
if ~all(cellfun(@(v,t) size(v,2) == length(t),wf_V_raw,wf_t_all))
%     wf_V_raw = cellfun(@(v,t) v(:,1:length(t)),wf_V_raw,wf_t_all,'uni',false);
%     disp('Widefield: timeline/frame number mismatch, cutting end');
    error('Widefield: timeline/frame number mismatch');
end

% Get framerate from timestamps
wf_framerate = mean(1./diff(wf_t_all{1}));

% Correct hemodynamics
[V_neuro_hemocorr,hemocorr_t] = plab.wf.hemo_correct( ...
    wf_U_raw{1},wf_V_raw{1},wf_t_all{1}, ...
    wf_U_raw{2},wf_V_raw{2},wf_t_all{2});

% % High-pass filter (causal filter)
% highpassCutoff = 0.01; % Hz
% [b100s, a100s] = butter(2, highpassCutoff/(wf_framerate/2), 'high');
% fV_neuro_hemocorr = filter(b100s,a100s,V_neuro_hemocorr,[],2);

% Get DF/F
wf_Vdf = plab.wf.svd_dff(wf_U_raw{1},V_neuro_hemocorr,wf_avg_all{1});

% Deconvolve
wf_Vdf_deconv = ap.wf_deconv(wf_Vdf,wf_framerate);

% Set final processed widefield variables
wf_U = wf_U_raw{1};
wf_V = wf_Vdf_deconv;
wf_t = hemocorr_t;
wf_avg = wf_avg_all{1};

% Align widefield (if alignment exists)
try
    wf_avg = ap.wf_align(wf_avg,animal,rec_day);
    wf_U = ap.wf_align(wf_U,animal,rec_day);
    if verbose; disp('Aligned widefield U/avg...');end
catch  me
    warning(me.identifier,'Widefield: %s',me.message);
end










