
animal = 'AP005';

use_workflow = {'lcr_passive'};
% use_workflow = {'lcr_passive_fullscreen'};
% use_workflow = {'lcr_passive','lcr_passive_fullscreen'};
% use_workflow = {'stim_wheel_right_stage1','stim_wheel_right_stage2'};
% use_workflow = 'sparse_noise';

recordings = ap.find_recordings(animal,use_workflow);

% use_rec = 1;
rec_day = '2023-06-22';
use_rec = strcmp(rec_day,{recordings.day});
% use_rec = length(recordings);

rec_day = recordings(use_rec).day;
rec_time = recordings(use_rec).protocol{end};

verbose = true;

%% Define what to load

if verbose; fprintf('Loading %s, %s, Protocol %s\n', animal, rec_day, rec_time); end;

% If nothing specified, load everything (but not LFP)
if ~exist('load_parts','var')
    load_parts.mousecam = true;
    load_parts.widefield = true;
    load_parts.ephys = true;
else
    % If only some things specified, don't load others
    if ~isfield(load_parts,'mousecam')
        load_parts.mousecam = false;
    end
    if ~isfield(load_parts,'widefield')
        load_parts.widefield = false;
    end
    if ~isfield(load_parts,'ephys')
        load_parts.ephys = false;
    end
end


%% Load timelite and associated inputs

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

% Wheel position and velocity
timelite_wheel_idx = strcmp({timelite.daq_info.channel_name}, 'wheel');
wheel_position = timelite.data(:,timelite_wheel_idx);
[wheel_velocity,wheel_move] = AP_parse_wheel(wheel_position,timelite.daq_info(timelite_wheel_idx).rate);

% Screen on times
screen_idx = strcmp({timelite.daq_info.channel_name}, 'stim_screen');
screen_on = timelite.data(:,screen_idx) > ttl_thresh;

% Photodiode flips (interpolate from previous across screen flicker)
photodiode_thresh_level = 1; % low: has a relatively slow rise time
photodiode_idx = strcmp({timelite.daq_info.channel_name}, 'photodiode');
photodiode_thresh_screen_on = medfilt1(timelite.data(screen_on,photodiode_idx),3);
photodiode_thresh = interp1(timelite.timestamps(screen_on),photodiode_thresh_screen_on, ...
    timelite.timestamps,'previous','extrap') > photodiode_thresh_level;
photodiode_times = timelite.timestamps(find(diff(photodiode_thresh) ~= 0) + 1);

% Reward times (if on past certain time, valve signal flips rapidly to
% avoid burnout - take the reward onset as flipping up and staying high for
% some length of samples)
reward_idx = strcmp({timelite.daq_info.channel_name}, 'reward_valve');
reward_thresh = timelite.data(:,reward_idx) >= ttl_thresh;
reward_on_pattern = [0,ones(1,3)]; % flip up, be consecutively high
reward_times = timelite.timestamps(strfind(reward_thresh',reward_on_pattern));

%% Load Bonsai

if verbose; disp('Loading Bonsai...'); end

% Get Bonsai workflow
bonsai_dir = dir(plab.locations.make_server_filename(animal,rec_day,rec_time,'bonsai'));
bonsai_workflow = bonsai_dir([bonsai_dir.isdir] & ~contains({bonsai_dir.name},'.')).name;

% Load Bonsai events
bonsai_events_fn = plab.locations.make_server_filename( ...
    animal,rec_day,rec_time,'bonsai','bonsai_events.csv');
if exist(bonsai_events_fn,'file')
    trial_events = AP_load_bonsai(bonsai_events_fn);
end

% Sparse noise: get noise locations and times
if strcmp(bonsai_workflow,'sparse_noise')
    bonsai_noise_fn = plab.locations.make_server_filename( ...
        animal,rec_day,rec_time,'bonsai','NoiseLocations.bin');
    fid = fopen(bonsai_noise_fn);

    n_x_squares = trial_events.parameters.ScreenExtentX./trial_events.parameters.StimSize;
    n_y_squares = trial_events.parameters.ScreenExtentY./trial_events.parameters.StimSize;

    noise_locations = reshape(fread(fid),n_y_squares,n_x_squares,[]);
    fclose(fid);

    % Get stim times from photodiode (extrapolate: sparse noise photodiode
    % flips every N stim to give a more robust signal)
    photodiode_stim_idx = 1:trial_events.parameters.NthPhotodiodeFlip:size(noise_locations,3);
    % (check that the number of photodiode flips is expected)
    if length(photodiode_stim_idx) ~= length(photodiode_times)
        if length(photodiode_stim_idx) < length(photodiode_times)
            % (rarely: Bonsai square to black temporarily, don't know why)
            % (fix? find when time differences on either side are less than a
            % threshold and remove those flips)
            photodiode_diff = diff(photodiode_times);
            photodiode_diff_thresh = mean(photodiode_diff)/2;
            bad_photodiode_idx = ...
                find(photodiode_diff(1:end-1) < photodiode_diff_thresh & ...
                photodiode_diff(2:end) < photodiode_diff_thresh) + 1;

            if (length(photodiode_times) - length(photodiode_stim_idx)) == ...
                    length(bad_photodiode_idx)
                % (if detected bad flips even out the numbers, apply fix)
                photodiode_times(bad_photodiode_idx) = [];
            else
                % (otherwise, error)
                error('Sparse noise: photodiode > stim, unfixable')
            end
        else
            error('Sparse noise: photodiode < stim')
        end
    end

    stim_times = interp1(photodiode_stim_idx,photodiode_times, ...
        1:size(noise_locations,3),'linear','extrap')';
end


%%%%% testing: bonsai times into timelite times
% bonsai_reward_datetime = [trial_events.timestamps([trial_events.values.Outcome] == 1).Outcome];
% bonsai_relative_t = bonsai_reward_datetime(1);
% 
% bonsai_reward_t = seconds(bonsai_reward_datetime - bonsai_relative_t);
% 
% % event_datetime = cellfun(@(x) x(1),{trial_events.timestamps.StimOn}');
% event_datetime = vertcat(trial_events.timestamps.QuiescenceStart);
% % event_datetime = vertcat(trial_events.timestamps.QuiescenceReset);
% 
% event_t_bonsai = seconds(event_datetime - bonsai_relative_t);
% 
% event_t_tl = interp1(bonsai_reward_t,reward_times,event_t_bonsai,'linear','extrap');



%% Load mousecam

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

%% Load widefield

if load_parts.widefield && ...
        exist(plab.locations.make_server_filename(animal,rec_day,[],'widefield'),'dir')

    if verbose; disp('Loading Widefield...'); end

    % Load widefield data for all colors
    widefield_colors = {'blue','violet'};
    [wf_avg_all,wf_U_raw,wf_V_raw,wf_t_all] = deal(cell(length(widefield_colors),1));
    for curr_wf = 1:length(widefield_colors)
        mean_image_fn = plab.locations.make_server_filename(animal,rec_day,[], ...
            'widefield',sprintf('meanImage_%s.npy',widefield_colors{curr_wf}));
        svdU_fn = plab.locations.make_server_filename(animal,rec_day,[], ...
            'widefield',sprintf('svdSpatialComponents_%s.npy',widefield_colors{curr_wf}));
        svdV_fn = plab.locations.make_server_filename(animal,rec_day,rec_time, ...
            'widefield',sprintf('svdTemporalComponents_%s.npy',widefield_colors{curr_wf}));

        wf_avg_all{curr_wf} = readNPY(mean_image_fn);
        wf_U_raw{curr_wf} = readNPY(svdU_fn);
        wf_V_raw{curr_wf} = readNPY(svdV_fn);

        % Timestamps: assume colors go in order (dictated by Arduino)
        wf_t_all{curr_wf} = widefield_expose_times(curr_wf:length(widefield_colors):end);

    end

    % Correct hemodynamics
    V_neuro_hemocorr = plab.wf.hemo_correct( ...
        wf_U_raw{1},wf_V_raw{1},wf_t_all{1}, ...
        wf_U_raw{2},wf_V_raw{2},wf_t_all{2});

    % Get DF/F
    wf_Vdf = plab.wf.svd_dff(wf_U_raw{1},V_neuro_hemocorr,wf_avg_all{1});

    % Deconvolve
    wf_framerate = mean(1./diff(wf_t_all{1}));
    wf_Vdf_deconv = AP_deconv_wf(wf_Vdf,[],wf_framerate);

    % Set final processed widefield variables
    wf_U = wf_U_raw{1};
    wf_V = wf_Vdf_deconv;
    wf_times = wf_t_all{1};
    wf_avg = wf_avg_all{1};

end

%% Load ephys

ephys_quality_control = false;

if load_parts.ephys

    if verbose; disp('Loading Ephys...'); end

    % Get paths for ephys, raw data, kilosort
    ephys_path = plab.locations.make_server_filename(animal,rec_day,[],'ephys');
    kilosort_path = fullfile(ephys_path,'pykilosort');
    open_ephys_path_dir = dir(fullfile(ephys_path,'experiment*','recording*'));
    open_ephys_path = fullfile(open_ephys_path_dir.folder,open_ephys_path_dir.name);
        
    % These are the digital channels going into the FPGA
    flipper_sync_idx = 1;
    
    % Load phy sorting if it exists
    % (old = cluster_groups.csv, new = cluster_group.tsv)
    cluster_filepattern = [kilosort_path filesep 'cluster_group*'];
    cluster_filedir = dir(cluster_filepattern);
    if ~isempty(cluster_filedir)
        cluster_filename = [kilosort_path filesep cluster_filedir.name];
        fid = fopen(cluster_filename);
        cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
        fclose(fid);
    end
    
    % Read header information
    header_path = [kilosort_path filesep 'dat_params.txt'];
    header_fid = fopen(header_path);
    header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
    fclose(header_fid);
    
    header = struct;
    for i = 1:length(header_info{1})
        header.(header_info{1}{i}) = header_info{2}{i};
    end
    
    % Load raw timestamps and sync
    open_ephys_timestamps = readNPY(fullfile(open_ephys_path, ...
        'continuous','Neuropix-3a-100.Neuropix-3a-AP','timestamps.npy'));
    open_ephys_ttl_states = readNPY(fullfile(open_ephys_path, ...
        'events','Neuropix-3a-100.Neuropix-3a-AP','TTL','states.npy'));
    open_ephys_ttl_timestamps = readNPY(fullfile(open_ephys_path, ...
        'events','Neuropix-3a-100.Neuropix-3a-AP','TTL','timestamps.npy'));
    open_ephys_ttl_flipper_idx = abs(open_ephys_ttl_states) == flipper_sync_idx;
    open_ephys_flipper.value = sign(open_ephys_ttl_states(open_ephys_ttl_flipper_idx));
    open_ephys_flipper.timestamps = open_ephys_ttl_timestamps(open_ephys_ttl_flipper_idx);

    % Load kilosort data
    if isfield(header,'sample_rate')
        ephys_sample_rate = str2num(header.sample_rate);
    elseif isfield(header,'ap_sample_rate')
        ephys_sample_rate = str2num(header.ap_sample_rate);
    end
    % (spike times: index Open Ephys timestamps rather than assume constant
    % sampling rate as before, this accounts for potentially dropped data)
    spike_times = open_ephys_timestamps(readNPY(fullfile(kilosort_path,'spike_times.npy')));
    spike_templates_0idx = readNPY(fullfile(kilosort_path,'spike_templates.npy'));
    templates_whitened = readNPY(fullfile(kilosort_path,'templates.npy'));
    channel_positions = readNPY(fullfile(kilosort_path,'channel_positions.npy'));
    channel_map = readNPY(fullfile(kilosort_path,'channel_map.npy'));
    winv = readNPY(fullfile(kilosort_path,'whitening_mat_inv.npy'));
    template_amplitudes = readNPY(fullfile(kilosort_path,'amplitudes.npy'));
    
    % Default channel map/positions are from end: make from surface
    % (hardcode this: kilosort2 drops channels)
    max_depth = 3840;
    channel_positions(:,2) = max_depth - channel_positions(:,2);
    
    % Unwhiten templates
    templates = zeros(size(templates_whitened));
    for t = 1:size(templates_whitened,1)
        templates(t,:,:) = squeeze(templates_whitened(t,:,:))*winv;
    end
    
    % Get the waveform of all templates (channel with largest amplitude)
    [~,max_site] = max(max(abs(templates),[],2),[],3);
    templates_max = nan(size(templates,1),size(templates,2));
    for curr_template = 1:size(templates,1)
        templates_max(curr_template,:) = ...
            templates(curr_template,:,max_site(curr_template));
    end
    waveforms = templates_max;
    
    % Get depth of each template
    % (get min-max range for each channel)
    template_chan_amp = squeeze(range(templates,2));
    % (zero-out low amplitude channels)
    template_chan_amp_thresh = max(template_chan_amp,[],2)*0.5;
    template_chan_amp_overthresh = template_chan_amp.*(template_chan_amp >= template_chan_amp_thresh);
    % (get center-of-mass on thresholded channel amplitudes)
    template_depths = sum(template_chan_amp_overthresh.*channel_positions(:,2)',2)./sum(template_chan_amp_overthresh,2);
    
    % Get the depth of each spike (templates are zero-indexed)
    spike_depths = template_depths(spike_templates_0idx+1);
    
    % Get trough-to-peak time for each template
    templates_max_signfix = bsxfun(@times,templates_max, ...
        sign(abs(min(templates_max,[],2)) - abs(max(templates_max,[],2))));
    
    [~,waveform_trough] = min(templates_max,[],2);
    [~,waveform_peak_rel] = arrayfun(@(x) ...
        max(templates_max(x,waveform_trough(x):end),[],2), ...
        transpose(1:size(templates_max,1)));
    waveform_peak = waveform_peak_rel + waveform_trough;
    
    templateDuration = waveform_peak - waveform_trough;
    templateDuration_us = (templateDuration/ephys_sample_rate)*1e6;
    
    % Get sync points for alignment
    
    % Get index of this experiment within day 
    % (to determine which flipper boundaries match this recording)
    experiment_idx = find(strcmp(ap.find_recordings(animal,rec_day).protocol,rec_time));

    if exist('flipper_times','var')
        % (if flipper, use that)
        % (at least one experiment the acqLive connection to ephys was bad
        % so it was delayed - ideally check consistency since it's
        % redundant)
        bad_flipper = false;
        
        % Get flipper experiment differences by long delays
        % (note: this is absolute difference, if recording stopped and
        % started then the clock starts over again, although I thought it
        % wasn't supposed to when I grab the concatenated sync, so
        % something might be wrong)
        flip_diff_thresh = 10; % time between flips to define experiment gap (s)
        flipper_expt_idx = [1;find(abs(diff(open_ephys_flipper.timestamps)) > ...
            flip_diff_thresh)+1;length(open_ephys_flipper.timestamps)+1];
        
        flipper_flip_times_ephys = open_ephys_flipper.timestamps( ...
            flipper_expt_idx(experiment_idx):flipper_expt_idx(experiment_idx+1)-1);
        
        % Pick flipper times to use for alignment
        if length(flipper_flip_times_ephys) == length(flipper_times)
            % If same number of flips in ephys/timeline, use all
            sync_timeline = flipper_times;
            sync_ephys = flipper_flip_times_ephys;
        elseif length(flipper_flip_times_ephys) ~= length(flipper_times)
            % If different number of flips in ephys/timeline, best
            % contiguous set via xcorr of diff
            warning([animal ' ' rec_day ':Flipper flip times different in timeline/ephys'])
            warning(['TEMPORARY: using only first and last flip'])
            sync_ephys = flipper_flip_times_ephys([1,end]);
            sync_timeline = flipper_times([1,end]);
        end
        
    else
        bad_flipper = true;
    end
        
    % Get spike times in timeline time
    spike_times_timeline = interp1(sync_ephys,sync_timeline,spike_times,'linear','extrap');
    
    % Get "good" templates from labels if quality control selected and
    % manual labels exist
    if ephys_quality_control && exist('cluster_groups','var')
        % If there's a manual classification
        if verbose; disp('Keeping manually labelled good units...'); end
        
        % Check that all used spike templates have a label
        spike_templates_0idx_unique = unique(spike_templates_0idx);
        if ~all(ismember(spike_templates_0idx_unique,uint32(cluster_groups{1}))) || ...
                ~all(ismember(cluster_groups{2},{'good','mua','noise'}))
            warning([animal ' ' day ': not all templates labeled']);
        end
        
        % Define good units from labels
        good_templates_idx = uint32(cluster_groups{1}( ...
            strcmp(cluster_groups{2},'good') | strcmp(cluster_groups{2},'mua')));
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);        
    else
        % If no cluster groups at all, keep all
        warning([animal ' ' rec_day ' - no ephys quality control']);
        if verbose; disp('No ephys quality control, keeping all and re-indexing'); end
        good_templates_idx = unique(spike_templates_0idx);
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);
    end
    
    % Throw out all non-good template data
    templates = templates(good_templates,:,:);
    template_depths = template_depths(good_templates);
    waveforms = waveforms(good_templates,:);
    templateDuration = templateDuration(good_templates);
    templateDuration_us = templateDuration_us(good_templates);
    
    % Throw out all non-good spike data
    good_spike_idx = ismember(spike_templates_0idx,good_templates_idx);
    spike_times = spike_times(good_spike_idx);
    spike_templates_0idx = spike_templates_0idx(good_spike_idx);
    template_amplitudes = template_amplitudes(good_spike_idx);
    spike_depths = spike_depths(good_spike_idx);
    spike_times_timeline = spike_times_timeline(good_spike_idx);
    
    % Rename the spike templates according to the remaining templates
    % (and make 1-indexed from 0-indexed)
    new_spike_idx = nan(max(spike_templates_0idx)+1,1);
    new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
    spike_templates = new_spike_idx(spike_templates_0idx+1);
       
end


%% Experiment scroller


% AP_expscroll(wf_U_raw{1},AP_deconv_wf(wf_V_raw{1},[],wf_framerate),wf_t_all{1},mousecam_fn,mousecam_times)

% AP_expscroll(wf_U_raw{1},wf_V_raw{1},wf_t_all{1},mousecam_fn,mousecam_times)
% AP_expscroll(wf_U_raw{2},wf_V_raw{2},wf_t_all{2},mousecam_fn,mousecam_times)

% AP_expscroll(wf_U,wf_Vdf,wf_times,mousecam_fn,mousecam_times)

% AP_expscroll(wf_U,wf_V,wf_times,mousecam_fn,mousecam_times)

% AP_expscroll(wf_U_raw{1},wf_V_raw{1},wf_t_all{1})
% AP_expscroll(wf_U_raw{2},wf_V_raw{2},wf_t_all{2})
% AP_expscroll(wf_U,wf_V,wf_times)












