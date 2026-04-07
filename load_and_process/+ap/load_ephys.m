% Load electrophysiology data

% Set parent ephys path
ephys_path = plab.locations.filename('server',animal,rec_day,[],'ephys');

% Find latest kilosort version
kilosort_versions = ["kilosort4","pykilosort"]; % in preference order
for curr_kilosort_version = kilosort_versions
    if exist(fullfile(ephys_path,curr_kilosort_version),'dir')
        kilosort_folder = curr_kilosort_version;
        break
    end
end

if ~exist('kilosort_folder','var') || isempty(kilosort_folder)
    error('%s %s: no kilosort folder found',animal,rec_day);
end

% Get path for raw data and kilosort for given recording
kilosort_top_path = fullfile(ephys_path,kilosort_folder);

if verbose; fprintf('Loading Ephys (%s)...\n',kilosort_folder); end

%% Load and prepare Kilosort data

% Set default probe to load
if ~exist('load_probe','var')
    load_probe = 1;
end

% Get OE folders
open_ephys_data_dir = dir(fullfile(ephys_path,'/**/','continuous.dat'));
open_ephys_recording_paths = string({open_ephys_data_dir(~contains({open_ephys_data_dir.folder}, ...
    {'ADC','LFP'},'IgnoreCase',true)).folder});

oe_recordings = unique(extract(open_ephys_recording_paths,'recording'+digitsPattern));

if ~isscalar(open_ephys_recording_paths) && length(open_ephys_recording_paths) == length(oe_recordings)
    % Multiple recordings: assume concatenated, use all folders
    open_ephys_path = open_ephys_recording_paths;
else
    % Single recording or muliple non-recording folders: separate probes
    open_ephys_path = open_ephys_recording_paths(load_probe);
end

% Get kilosort output folder
kilosort_dir = dir(kilosort_top_path);
if any(contains({kilosort_dir([kilosort_dir.isdir]).name},'probe'))
    % Multi-probe: nested folder kilosort/probe_n
    kilosort_path = fullfile(kilosort_top_path,sprintf('probe_%d',load_probe));
else
    % Single probe: data in kilosort top path
    kilosort_path = kilosort_top_path;
end

%%%%% OLD - for serial multi-site
% (not currently used: multi-site with same probe = /site_n)
if any(contains({kilosort_dir.name},'site'))
    % If 'site' folders (recordings in serial), choose last before recording
    ephys_site_paths = dir(fullfile(kilosort_top_path,'site_*'));
    if ~isempty(ephys_site_paths)
        ephys_site_datetime = NaT(size(ephys_site_paths));
        for curr_site = 1:length(ephys_site_paths)
            curr_ephys_settings_file = fullfile( ...
                ephys_path,ephys_site_paths(curr_site).name, ...
                'settings.xml');
            curr_ephys_settings = readstruct(curr_ephys_settings_file,'filetype','xml');
            ephys_site_datetime(curr_site) = ...
                datetime(curr_ephys_settings.INFO.DATE, ...
                'InputFormat','dd MMM yyyy HH:mm:ss');
        end

        % (add 1 minute leeway to recording time since no seconds)
        ephys_use_site = find(ephys_site_datetime - ...
            (rec_datetime + minutes(1)) < 0,1,'last');

        kilosort_path = fullfile(kilosort_top_path, ...
            ephys_site_paths(ephys_use_site).name);

        open_ephys_recording_paths = dir(fullfile(ephys_path, ...
            ephys_site_paths(ephys_use_site).name,'experiment*', ...
            'recording*','continuous','*-AP'));
        open_ephys_path = fullfile({open_ephys_recording_paths.folder},{open_ephys_recording_paths.name});
    end
end

% % Get start time of ephys recording (unused currently)
% ephys_settings_filename = fullfile(fileparts(open_ephys_dir(1).folder),'settings.xml');
% ephys_settings = readstruct(ephys_settings_filename,'filetype','xml');
% ephys_datetime = datetime(ephys_settings.INFO.DATE, ...
%     'InputFormat','dd MMM yyyy HH:mm:ss');
%%%%%%

% Load Open Ephys metadata (for sample rate)
oe_metadata_fn = fullfile(fileparts(fileparts(open_ephys_path{1})),'structure.oebin');
oe_metadata = jsondecode(fileread(oe_metadata_fn));
oe_ap_samplerate = oe_metadata(1).continuous(1).sample_rate;

% Load kilosort data
% (spike times: load open ephys times - create if not yet created)
spike_times_openephys_filename = fullfile(kilosort_path,'spike_times_openephys.npy');
if ~exist(spike_times_openephys_filename,'file')
    ks_spike_times_fn = fullfile(kilosort_path,'spike_times.npy');

    oe_samples_dir = cellfun(@(data_path) ...
        dir(fullfile(data_path,'sample_numbers.npy')), ...
        open_ephys_path,'uni',false);
    oe_samples_fn = cellfun(@(data_dir) ...
        fullfile(data_dir.folder,data_dir.name),oe_samples_dir,'uni',false);

    plab.ephys.ks2oe_timestamps(ks_spike_times_fn,oe_samples_fn, ...
        oe_metadata(1).continuous(1).sample_rate);
end
spike_times_openephys = readNPY(spike_times_openephys_filename);

% (spike times: index Open Ephys timestamps rather than assume constant
% sampling rate as before, this accounts for potentially dropped data)
spike_templates = readNPY(fullfile(kilosort_path,'spike_clusters.npy'))+1; % (convert 0-idx to 1-idx)
templates_whitened = readNPY(fullfile(kilosort_path,'templates.npy'));
channel_positions = readNPY(fullfile(kilosort_path,'channel_positions.npy'));
channel_map = readNPY(fullfile(kilosort_path,'channel_map.npy'));
winv = readNPY(fullfile(kilosort_path,'whitening_mat_inv.npy'));
template_amplitudes = double(readNPY(fullfile(kilosort_path,'amplitudes.npy')));

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
waveforms = cell2mat(arrayfun(@(x) templates(x,:,max_site(x)), ...
    (1:size(templates,1))','uni',false));

% Get depth of each template
% (get min-max range for each channel)
template_chan_amp = squeeze(range(templates,2));
% (zero-out low amplitude channels)
template_chan_amp_thresh = max(template_chan_amp,[],2)*0.5;
template_chan_amp_overthresh = template_chan_amp.*(template_chan_amp >= template_chan_amp_thresh);
% (get center-of-mass on thresholded channel amplitudes)
template_depths = sum(template_chan_amp_overthresh.*channel_positions(:,2)',2)./sum(template_chan_amp_overthresh,2);

% Get the depth of each spike
spike_depths = template_depths(spike_templates);

% Get waveform width 
ephys_sample_rate = 30000; % (just hardcoded for now, it never changes)

% (use smoothed waveform - Kilosort often has bumps)
waveforms_movmean = movmean(waveforms,3,2);

% 1) trough-to-peak
[~,waveform_trough] = min(waveforms_movmean,[],2);
[~,waveform_peak_rel] = arrayfun(@(x) ...
    max(waveforms_movmean(x,waveform_trough(x):end),[],2), ...
    transpose(1:size(waveforms_movmean,1)));
waveform_peak = waveform_peak_rel + waveform_trough - 1;
waveform_duration_peaktrough = ...
    1e6*(waveform_peak - waveform_trough)/ephys_sample_rate;

% 2) full width half max
waveform_duration_fwhm = arrayfun(@(x) ...
    sum(waveforms_movmean(x,:) <= min(waveforms_movmean(x,:))/2) * ...
    1e6/ephys_sample_rate,1:size(templates,1));

%% Check for Open Ephys dropped data

%%%%%%%%%%%%% UNDER CONSTRUCTION
% 
% At the moment - just gives a warning
% 
% Examples where this happens: 
% Many drops: DS031 2026-03-29
% One drop: PG003 2026-03-26

oe_sync_messages_fn = fullfile(fileparts(fileparts(open_ephys_path{1})),"sync_messages.txt");
oe_sync_messages = regexp(readlines(oe_sync_messages_fn), ...
    '- (?<stream>.*) @ (?<sample_rate>.*) Hz: (?<first_sample>\d*)','names');
oe_stream_info = vertcat(oe_sync_messages{:});

oe_messages_fn = fullfile(fileparts(fileparts(open_ephys_path{1})), ...
    "events","MessageCenter","text.npy");
oe_message_center_text = readlines(oe_messages_fn);
oe_bad_samples_string = regexp(oe_message_center_text(2), ...
    'NPX TIMESTAMP JUMP: (?<time_jump>\d+).*?sample number (?<sample>\d+)','names');
if ~isempty(oe_bad_samples_string)
    % (convert strings: 'timestamp jump' is 100khz clock on Npx headstage)
    oe_bad_samples = struct('time_jump', ...
        num2cell(str2double(vertcat(oe_bad_samples_string.time_jump))/1e5), ...
        'sample', ...
        num2cell(str2double(vertcat(oe_bad_samples_string.sample))));

    warning('%s %s ephys had %d drops = %.4fs',animal,rec_day, ...
        length(oe_bad_samples),sum(vertcat(oe_bad_samples.time_jump)));

    % % drop ephys flips that happen on the bad sample
    % drop_ephys_flips = ismember(vertcat(open_ephys_ttl_sample_numbers{:}), ...
    %     str2double(vertcat(oe_bad_samples.sample)));
    %
    % % drop timelite flips that are within dropped time
    % (this is going to be harder: each drop means it needs to be re-aligned)
    % % (if I'm lucky - I could just iteratively remove these? would need to be
    % lag-aligned first though)
    % % vertcat(oe_bad_samples.sample)
end


%% Convert timestamps to timelite (with flipper)
% (allow for multiple files if multiple recordings concatenated)

% Load ephys flipper
% (use sample numbers as accurate, convert into timestamps -
% 'timestamps.npy' are recomputed across recording and not always accurate)
flipper_sync_idx = 1;

open_ephys_ttl_path = fullfile(strrep(open_ephys_path, ...
    'continuous','events'),'TTL');

% Load TTL sample numbers
open_ephys_ttl_sample_numbers = cellfun(@(data_path) ...
    readNPY(fullfile(data_path,'sample_numbers.npy')), ...
    open_ephys_ttl_path,'uni',false)';

% Check for clock resets as backwards TTL timesteps across recordings
open_ephys_ttl_sample_backstep = ...
    find(cellfun(@(x) x(1),open_ephys_ttl_sample_numbers(2:end)) - ...
    cellfun(@(x) x(end),open_ephys_ttl_sample_numbers(1:end-1)) < 0) + 1;

if isempty(open_ephys_ttl_sample_backstep)
    % If no clock resets, just concatenate
    open_ephys_ttl_timestamps = ...
        double(vertcat(open_ephys_ttl_sample_numbers{:}))/oe_ap_samplerate;
else
    % If clock resets, make pseudocontinuous to match ks2oe
    % (load OE samples)
    oe_samples_dir = cellfun(@(data_path) ...
        dir(fullfile(data_path,'/**/','sample_numbers.npy')), ...
        open_ephys_path,'uni',false);
    oe_samples_fns = cellfun(@(data_dir) ...
        fullfile(data_dir.folder,data_dir.name),oe_samples_dir,'uni',false);
    oe_samples_split = cellfun(@readNPY,oe_samples_fns,'uni',false);
    oe_recordings_last_samples = cellfun(@(x) x(end),oe_samples_split);

    open_ephys_ttl_sample_numbers_pseudocontinuous = ...
        open_ephys_ttl_sample_numbers;
    open_ephys_ttl_sample_numbers_pseudocontinuous(open_ephys_ttl_sample_backstep) = ...
        cellfun(@(x,sample_add) x + sample_add, ...
        open_ephys_ttl_sample_numbers(open_ephys_ttl_sample_backstep), ...
        num2cell(oe_recordings_last_samples(open_ephys_ttl_sample_backstep-1)), ...
        'uni',false);

    open_ephys_ttl_timestamps = ...
        double(vertcat(open_ephys_ttl_sample_numbers_pseudocontinuous{:}))/oe_ap_samplerate;
end

open_ephys_ttl_states = cell2mat(cellfun(@(data_path) ...
    readNPY(fullfile(data_path,'states.npy')),open_ephys_ttl_path,'uni',false)');

open_ephys_ttl_flipper_idx = abs(open_ephys_ttl_states) == flipper_sync_idx;
open_ephys_flipper.value = sign(open_ephys_ttl_states(open_ephys_ttl_flipper_idx));
open_ephys_flipper.timestamps = open_ephys_ttl_timestamps(open_ephys_ttl_flipper_idx);

% Resample Open Ephys flipper to DAQ sample rate
open_ephys_flipper_trace_t = open_ephys_flipper.timestamps(1): ...
    1/timelite.daq_info(1).rate: open_ephys_flipper.timestamps(end);
open_ephys_flipper_trace = logical(normalize(interp1( ...
    open_ephys_flipper.timestamps, ...
    single(open_ephys_flipper.value), ...
    open_ephys_flipper_trace_t,'previous'),'range'));

% Get Open Ephys corresponding to timelite flipper
% (get lag between timelite and ephys by correlation)
ephys_timelite_flipper_lag = finddelay(+flipper_thresh, ...
    +open_ephys_flipper_trace)/timelite.daq_info(1).rate + ...
    open_ephys_flipper.timestamps(1);
% (get all ephys flips within the matching continuous window)
curr_ephys_flipper_idx =  ...
    open_ephys_flipper.timestamps >= ephys_timelite_flipper_lag & ...
    open_ephys_flipper.timestamps <= ephys_timelite_flipper_lag + ...
    length(flipper_thresh)/timelite.daq_info(1).rate;
% (set equivalent flips for ephys/timelite)
flipper_flip_times_ephys = open_ephys_flipper.timestamps(curr_ephys_flipper_idx);

% Pick flipper times to use for alignment and handle problems
if length(flipper_flip_times_ephys) == length(flipper_times)
    % If same number of flips in ephys/timelite, use all
    sync_timelite = flipper_times;
    sync_ephys = flipper_flip_times_ephys;

elseif length(flipper_flip_times_ephys) < length(flipper_times)
    % If more flips in timelite than ephys, assume timelite got all flips
    % and only subset of flips caught in ephys: find usable subset of TL
    % flips

    % Estimate nearest flip in ephys for each flip in timelite
    flip_timelite_ephys_timediff = ...
        interp1(flipper_flip_times_ephys,flipper_flip_times_ephys, ...
        flipper_times + ephys_timelite_flipper_lag,'nearest','extrap') - ...
        (flipper_times+ephys_timelite_flipper_lag);

    % Find the first flip that was out of sync
    flip_timediff_thresh = 0.01;
    last_good_flip_idx = find(abs(flip_timelite_ephys_timediff) > flip_timediff_thresh,1)-2;

    % Use flips that occured before desync
    sync_timelite = flipper_times(1:last_good_flip_idx);
    sync_ephys = flipper_flip_times_ephys(1:last_good_flip_idx);

    %%%%%%%% UNDER CONSTRUCTION
    % DS031 2026-03-29 1249: ephys missed ~7s which it doesn't know about,
    % need generalized solution to identify chunks of missed flips

    % For future: generalized solution to look for gaps could be
    % truncating from the unmatched flips, finding delay, matching,
    % repeating. This finds the delay in the first instance, but haven't
    % gotten to specifically removing the timelite-only flips yet.
    
    % m_t = finddelay(+flipper_thresh(timelite.timestamps>flipper_times(first_bad_flip_idx)), ...
    %     +open_ephys_flipper_trace(open_ephys_flipper_trace_t > ...
    %     open_ephys_flipper.timestamps(first_bad_flip_idx)))/timelite.daq_info(1).rate;
    %%%%%%%%%%%%%%

    % Warn missed flips: this solution may not be general yet
    warning('%s %s: Ephys missed flips, using first %d/%d (%.0f%%) flips',animal,rec_day, ...
        last_good_flip_idx,length(flipper_times),100*last_good_flip_idx/length(flipper_times));

elseif length(flipper_flip_times_ephys) > length(flipper_times)
    % If more flips in ephys than timelite, assume extraneous brief upward
    % flips from a bad connection or resting voltage: remove short upward
    % flips

    bad_ephys_flips_down = open_ephys_flipper.value(curr_ephys_flipper_idx) == -1 & ...
        [0;diff(flipper_flip_times_ephys)] < 0.05;

    use_ephys_flips = ~bad_ephys_flips_down & ~circshift(bad_ephys_flips_down,-1);

    if sum(use_ephys_flips) == length(flipper_times)
        % Successful if same number of matched flips:
        % use subset of ephys flips
        sync_ephys = flipper_flip_times_ephys(use_ephys_flips);
        % use all timelite flips
        sync_timelite = flipper_times;
    else
        % If not, unclear where the matched flips are, error out
        error('%s %s: flips ephys > timelite, cannot fix',animal,rec_day);
    end

end

% Get spike times in timelite time
spike_times_timelite = interp1(sync_ephys,sync_timelite,spike_times_openephys,'linear','extrap');


%% Load probe position 

% (load histology positions if available)
histology_probe_filename = dir(...
    plab.locations.filename('server',animal,[],[], ...
    'histology','*','probe_ccf.mat'));
if length(histology_probe_filename) > 1
    histology_probe_filename = histology_probe_filename(1);
    warning('Multiple ''probe_ccf.mat'' files, using from %s',histology_probe_filename(1).folder);
end
if ~isempty(histology_probe_filename)
   probe_histology = load(fullfile(histology_probe_filename.folder, ...
        histology_probe_filename.name));
end

% (load NTE positions if available)
nte_positions_filename = dir(fullfile(ephys_path,'**','*probe_positions*.mat'));
if ~isempty(nte_positions_filename)
    probe_nte = load(fullfile(nte_positions_filename.folder,nte_positions_filename.name));
end

% (use histology if available and matching day, use NTE if not)
if exist('probe_histology','var') && isfield(probe_histology.probe_ccf,'day') && ...
        any(find(strcmp(rec_day,{probe_histology.probe_ccf.day})))
    probe_histology_day_idx = find(strcmp(rec_day,{probe_histology.probe_ccf.day}));
    probe_areas = {probe_histology.probe_ccf(probe_histology_day_idx).trajectory_areas};
    if verbose; disp('Ephys: Loaded histology positions...'); end
elseif exist('probe_nte','var')
    probe_areas = probe_nte.probe_areas(load_probe);
    if verbose; disp('Ephys: Loaded NTE positions...'); end
end

%% Remove bad units from quality control

ephys_qc_type = 'bombcell';

if strcmp(ephys_qc_type,'bombcell')
    % Load Bombcell quality metrics (if exist)
    qMetrics_path = fullfile(kilosort_path,'qMetrics');
    if  exist(qMetrics_path,'dir')
        % Bombcell function to load all metrics, only used for troubleshooting
        %     [param, qMetric] = bc.load.loadSavedMetrics(qMetrics_path);
        % (native bombcell labels)
        %  unitType = bc_getQualityUnitType(param, qMetric);
        % (extra metrics & translated bombcell labels)
        load(fullfile(qMetrics_path, 'template_qc_labels.mat'))

        % Define good units from labels
        good_templates = ismember(template_qc_labels,{'singleunit','multiunit'});
        if verbose; fprintf('Ephys: Applying Bombcell quality metrics...'); end

        % Keep only labels from good units
        template_qc_labels = template_qc_labels(good_templates);

    else
        warning('Bombcell metrics not available');
        return
    end

elseif strcmp(ephys_qc_type,'phy')
    % Load manual Phy sorting (old)
    cluster_filepattern = [kilosort_path filesep 'cluster_group*'];
    cluster_filedir = dir(cluster_filepattern);
    if ~isempty(cluster_filedir)
        cluster_filename = [kilosort_path filesep cluster_filedir.name];
        fid = fopen(cluster_filename);
        cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
        fclose(fid);
    end
    % Get "good" templates from labels if quality control selected and
    % manual labels exist
    if exist('cluster_groups','var')
        % If there's a manual classification
        if verbose; disp('Ephys: Keeping manually labelled good units...'); end

        % Check that all used spike templates have a label
        spike_templates_unique_0idx = unique(spike_templates)-1;
        if ~all(ismember(spike_templates_unique_0idx,uint32(cluster_groups{1}))) || ...
                ~all(ismember(cluster_groups{2},{'good','mua','noise'}))
            warning([animal ' ' day ': not all templates labeled']);
        end

        % Define good units from labels
        good_templates_idx = uint32(cluster_groups{1}( ...
            strcmp(cluster_groups{2},'good') | strcmp(cluster_groups{2},'mua')));
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);
    else
        % If no cluster groups at all, keep all
        warning('Phy labels not available');
        return
    end
else
    return
end

% If good templates were selected above, throw out not-good data:

% Throw out all non-good template data
templates = templates(good_templates,:,:);
template_depths = template_depths(good_templates);
waveforms = waveforms(good_templates,:);
waveform_duration_peaktrough = waveform_duration_peaktrough(good_templates);
waveform_duration_fwhm = waveform_duration_fwhm(good_templates);

% Throw out all non-good spike data
good_spikes = ismember(spike_templates,find(good_templates));
spike_templates = spike_templates(good_spikes);
template_amplitudes = template_amplitudes(good_spikes);
spike_depths = spike_depths(good_spikes);
spike_times_timelite = spike_times_timelite(good_spikes);

% Rename the remaining spike templates (1:N, to match index for template)
[~,spike_templates] = ismember(spike_templates,find(good_templates));

if verbose; fprintf('kept %d/%d "good" units...\n',sum(good_templates),length(good_templates)); end
















