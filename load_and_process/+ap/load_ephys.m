% Load electrophysiology data

if verbose; disp('Loading Ephys...'); end

%% Set paths
 
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

% Get OE folders
open_ephys_data_dir = dir(fullfile(ephys_path,'/**/','continuous.dat'));
open_ephys_recording_paths = string({open_ephys_data_dir(~contains({open_ephys_data_dir.folder}, ...
    {'ADC','LFP'},'IgnoreCase',true)).folder});

% Set default probe to load
if ~exist('load_probe','var')
    load_probe = 1;
end

%% Load Kilosort data

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

% (spike times: index Open Ephys timestamps rather than assume constant
% sampling rate as before, this accounts for potentially dropped data)
spike_templates = readNPY(fullfile(kilosort_path,'spike_clusters.npy'))+1; % (convert 0-idx to 1-idx)
templates_whitened = readNPY(fullfile(kilosort_path,'templates.npy'));
channel_positions = readNPY(fullfile(kilosort_path,'channel_positions.npy'));
channel_map = readNPY(fullfile(kilosort_path,'channel_map.npy'));
winv = readNPY(fullfile(kilosort_path,'whitening_mat_inv.npy'));
template_amplitudes = double(readNPY(fullfile(kilosort_path,'amplitudes.npy')));

%% Calculate template properties

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
% (use smoothed waveform - Kilosort often has bumps)
waveforms_movmean = movmean(waveforms,3,2);

% 1) trough-to-peak
[~,waveform_trough] = min(waveforms_movmean,[],2);
[~,waveform_peak_rel] = arrayfun(@(x) ...
    max(waveforms_movmean(x,waveform_trough(x):end),[],2), ...
    transpose(1:size(waveforms_movmean,1)));
waveform_peak = waveform_peak_rel + waveform_trough - 1;
waveform_duration_peaktrough = ...
    1e6*(waveform_peak - waveform_trough)/oe_ap_samplerate;

% 2) full width half max
waveform_duration_fwhm = arrayfun(@(x) ...
    sum(waveforms_movmean(x,:) <= min(waveforms_movmean(x,:))/2) * ...
    1e6/oe_ap_samplerate,1:size(templates,1));

%% Load (or create) spike/TTL timestamps

% Set timestamp filenames (spike times and TTL)
spike_times_openephys_filename =  fullfile(kilosort_path,'spike_times_openephys.npy');
open_ephys_ttl_filename =  fullfile(kilosort_path,'open_ephys_ttl.mat');

% If timestamps don't exist, create them
if ~exist(spike_times_openephys_filename,'file') || ~exist(open_ephys_ttl_filename,'file')
    ks_spike_times_fn = fullfile(kilosort_path,'spike_times.npy');

    oe_samples_dir = cellfun(@(data_path) ...
        dir(fullfile(data_path,'sample_numbers.npy')), ...
        open_ephys_path,'uni',false);
    oe_samples_fn = cellfun(@(data_dir) ...
        fullfile(data_dir.folder,data_dir.name),oe_samples_dir,'uni',false);

    plab.ephys.ks2oe_timestamps(ks_spike_times_fn,oe_samples_fn, ...
        oe_metadata(1).continuous(1).sample_rate);
end

% Load timestamps
spike_times_openephys = readNPY(spike_times_openephys_filename);
load(open_ephys_ttl_filename);

%% Convert ephys to timelite times

% (allow for multiple files if multiple recordings concatenated)

% Set TTL index for flipper (sync)
flipper_sync_idx = 1;
open_ephys_ttl_flipper_idx = abs(open_ephys_ttl.state) == flipper_sync_idx;

% Resample Open Ephys flipper to DAQ sample rate
open_ephys_flipper_trace_t = open_ephys_ttl.timestamps(1): ...
    1/timelite.daq_info(1).rate: open_ephys_ttl.timestamps(end);
open_ephys_flipper_trace = logical(normalize(interp1( ...
    open_ephys_ttl.timestamps(open_ephys_ttl_flipper_idx), ...
    single(sign(open_ephys_ttl.state(open_ephys_ttl_flipper_idx))), ...
    open_ephys_flipper_trace_t,'previous'),'range'));

% Get Open Ephys corresponding to timelite flipper
% (get lag between timelite and ephys by correlation)
ephys_timelite_flipper_lag = finddelay(+flipper_thresh, ...
    +open_ephys_flipper_trace)/timelite.daq_info(1).rate + ...
    open_ephys_ttl.timestamps(1);
% (get all ephys flips within the matching continuous window)
curr_ephys_flipper_idx =  ...
    open_ephys_ttl_flipper_idx & ...
    open_ephys_ttl.timestamps >= ephys_timelite_flipper_lag & ...
    open_ephys_ttl.timestamps <= ephys_timelite_flipper_lag + ...
    length(flipper_thresh)/timelite.daq_info(1).rate;
% (set equivalent flips for ephys/timelite)
flipper_flip_times_ephys = ...
    open_ephys_ttl.timestamps(curr_ephys_flipper_idx);

% Pick flipper times to use for alignment
if length(flipper_flip_times_ephys) == length(flipper_times)
    % If same number of flips in ephys/timelite, use all
    sync_timelite = flipper_times;
    sync_ephys = flipper_flip_times_ephys;

else
    % If different number of timelite/ephys flips: 
    % find usable flips as "matched" with a very short estimated time delay

    % Estimate nearest flip for timelite->ephys and ephys->timelite
    flip_timelite2ephys_timediff = ...
        interp1(flipper_flip_times_ephys,flipper_flip_times_ephys, ...
        flipper_times + ephys_timelite_flipper_lag,'nearest','extrap') - ...
        (flipper_times+ephys_timelite_flipper_lag);

    flip_ephys2timelite_timediff = ...
        interp1((flipper_times+ephys_timelite_flipper_lag),(flipper_times+ephys_timelite_flipper_lag), ...
        flipper_flip_times_ephys,'nearest','extrap') - ...
        flipper_flip_times_ephys;

    % Set cutoff to find "matched" flips
    flip_timediff_thresh = 0.01;
    use_timelite_flips = abs(flip_timelite2ephys_timediff) < flip_timediff_thresh;
    use_ephys_flips = abs(flip_ephys2timelite_timediff) < flip_timediff_thresh;

    if sum(use_timelite_flips) == sum(use_ephys_flips)
        % Successful matching if same number of usable flips
        sync_ephys = flipper_flip_times_ephys(use_ephys_flips);
        sync_timelite = flipper_times(use_timelite_flips);
    else
        % If not, unclear where the matched flips are, error out
        error('%s %s: Cannot match timelite/ephys flips',animal,rec_day);
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

        % Load unit labels

        %%%%%%%%%%%%%% UNDER CONSTRUCTION: 
        % Going back to native bombcell labels to apply params on the fly,
        % makes changing params more robust without saving new labels any
        % time a lab-wide param is changed
        % 
        % - Decide whether/how to add raw/template corr param
        % - Bombcell function to load all metrics, only used for troubleshooting
        % [bombcell_param, bombcell_qMetric] = bc.load.loadSavedMetrics(qMetrics_path);
        % - Adjust these axon parameters: 
        % param.minWidthFirstPeak_nonSomatic = Inf;
        % param.maxPeak1ToPeak2Ratio_nonSomatic = 1/4;
        % - Turn of label saving (just loading for here)
        % bombcell_param.saveAsTSV = false;
        % - get bombcell unit type: 
        % unitType = bc.qm.getQualityUnitType(bombcell_param, bombcell_qMetric);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

        % (extra metrics & translated bombcell labels created in ap.run_bombcell)
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
















