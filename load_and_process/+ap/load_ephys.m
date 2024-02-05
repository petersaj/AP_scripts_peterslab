% Load electrophysiology data

ephys_path = plab.locations.filename('server',animal,rec_day,[],'ephys');
if verbose; disp('Loading Ephys...'); end

%% Load and prepare Kilosort data

% Get path for raw data and kilosort for given recording
kilosort_top_path = fullfile(ephys_path,'pykilosort');
kilosort_dir = dir(kilosort_top_path);

if ~any(contains({kilosort_dir.name},{'site','probe'}))
    % If no site/probe subfolders, use that top path
    kilosort_path = kilosort_top_path;
    open_ephys_path_dir = dir(fullfile(ephys_path,'experiment*','recording*'));
    open_ephys_path = fullfile(open_ephys_path_dir.folder,open_ephys_path_dir.name);

elseif any(contains({kilosort_dir.name},'probe'))
    % If 'probe' folders (recordings in parallel), choose last before recording
    ephys_site_paths = dir(fullfile(kilosort_top_path,'probe_*'));
    warning('Multiple probes: just loading probe 1 for now')
    kilosort_path = fullfile(ephys_site_paths(1).folder,ephys_site_paths(1).name);
    open_ephys_path_dir = dir(fullfile(ephys_path,ephys_site_paths(1).name,'experiment*','recording*'));
    open_ephys_path = fullfile(open_ephys_path_dir.folder,open_ephys_path_dir.name);

elseif any(contains({kilosort_dir.name},'site'))
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
        open_ephys_path_dir = dir(fullfile(ephys_path, ...
            ephys_site_paths(ephys_use_site).name,'experiment*','recording*'));
        open_ephys_path = fullfile(open_ephys_path_dir.folder,open_ephys_path_dir.name);
    end
end

% Get start time of ephys recording
ephys_settings_filename = fullfile(fileparts(open_ephys_path_dir.folder),'settings.xml');
ephys_settings = readstruct(ephys_settings_filename,'filetype','xml');
ephys_datetime = datetime(ephys_settings.INFO.DATE, ...
    'InputFormat','dd MMM yyyy HH:mm:ss');

% Load kilosort data
% (spike times: load open ephys times if available, create if not)
spike_times_openephys_filename = fullfile(kilosort_path,'spike_times_openephys.npy');
if exist(spike_times_openephys_filename,'file')
    spike_times_openephys = readNPY(spike_times_openephys_filename);
else
    open_ephys_timestamps_dir = dir(fullfile(open_ephys_path, ...
        'continuous','*-AP','timestamps.npy'));
    open_ephys_timestamps = readNPY(fullfile(open_ephys_timestamps_dir.folder,...
        open_ephys_timestamps_dir.name));
    spike_times_openephys = open_ephys_timestamps( ...
        readNPY(fullfile(kilosort_path,'spike_times.npy')));
    writeNPY(spike_times_openephys,spike_times_openephys_filename);
    fprintf('Converted and saved Open Ephys spike times: %s\n',spike_times_openephys_filename);
end

% (spike times: index Open Ephys timestamps rather than assume constant
% sampling rate as before, this accounts for potentially dropped data)
spike_templates = readNPY(fullfile(kilosort_path,'spike_templates.npy'))+1; % (convert 0-idx to 1-idx)
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

% Get the depth of each spike
spike_depths = template_depths(spike_templates);

% Get trough-to-peak time for each template
templates_max_signfix = bsxfun(@times,templates_max, ...
    sign(abs(min(templates_max,[],2)) - abs(max(templates_max,[],2))));

[~,waveform_trough] = min(templates_max,[],2);
[~,waveform_peak_rel] = arrayfun(@(x) ...
    max(templates_max(x,waveform_trough(x):end),[],2), ...
    transpose(1:size(templates_max,1)));
waveform_peak = waveform_peak_rel + waveform_trough;

ephys_sample_rate = 30000; % (just hardcoded for now, it never changes)
templateDuration = waveform_peak - waveform_trough;
templateDuration_us = (templateDuration/ephys_sample_rate)*1e6;

%% Convert timestamps to timelite (with flipper)

% Load ephys flipper
flipper_sync_idx = 1;

open_ephys_ttl_dir = dir(fullfile(open_ephys_path,'events', '*-AP','TTL'));
open_ephys_ttl_path = cell2mat(unique({open_ephys_ttl_dir.folder}));

open_ephys_ttl_states = readNPY(fullfile(open_ephys_ttl_path,'states.npy'));
open_ephys_ttl_timestamps = readNPY(fullfile(open_ephys_ttl_path,'timestamps.npy'));

open_ephys_ttl_flipper_idx = abs(open_ephys_ttl_states) == flipper_sync_idx;
open_ephys_flipper.value = sign(open_ephys_ttl_states(open_ephys_ttl_flipper_idx));
open_ephys_flipper.timestamps = open_ephys_ttl_timestamps(open_ephys_ttl_flipper_idx);

% Get sync points for alignment

% Resample Open Ephys flipper to DAQ sample rate
open_ephys_flipper_trace = logical(normalize(interp1(open_ephys_flipper.timestamps, ...
    single(open_ephys_flipper.value), ...
    open_ephys_flipper.timestamps(1):1/timelite.daq_info(1).rate: ...
    open_ephys_flipper.timestamps(end),'previous'),'range'));

% Get Open Ephys corresponding to timelite flipper
% (get lag between timelite and ephys by correlation)
ephys_timelite_flipper_lag = finddelay(+flipper_thresh, ...
    +open_ephys_flipper_trace)/timelite.daq_info(1).rate + ...
    open_ephys_flipper.timestamps(1);
% (get all ephys flips within that timelite window)
curr_ephys_flipper_idx =  ...
    open_ephys_flipper.timestamps >= ephys_timelite_flipper_lag & ...
    open_ephys_flipper.timestamps <= ephys_timelite_flipper_lag + ...
    length(flipper_thresh)/timelite.daq_info(1).rate;
flipper_flip_times_ephys = open_ephys_flipper.timestamps(curr_ephys_flipper_idx);

% Pick flipper times to use for alignment
if length(flipper_flip_times_ephys) == length(flipper_times)
    % If same number of flips in ephys/timeline, use all
    sync_timeline = flipper_times;
    sync_ephys = flipper_flip_times_ephys;

elseif length(flipper_flip_times_ephys) < length(flipper_times)
    warning('%s %s flips: %d ephys/%d timelite, using first/last for now', ...
        animal, rec_day, length(flipper_flip_times_ephys),length(flipper_times));
    sync_timeline = flipper_times([1,end]);
    sync_ephys = flipper_flip_times_ephys([1,end]);

    %         sync_timeline = flipper_times;
    %         sync_ephys = flipper_flip_times_ephys;
    %
    %         % Remove flips by flip interval until same number
    %         ephys_missed_flips_n = length(flipper_times) - length(flipper_flip_times_ephys);
    %         missed_flip_t_thresh = 0.1;
    %         while length(sync_ephys) ~= length(sync_timeline)
    %             curr_missed_flip = ...
    %                 find(abs(diff(sync_ephys) - ...
    %                 diff(sync_timeline(1:length(sync_ephys)))) > ...
    %                 missed_flip_t_thresh,1);
    %             sync_timeline(curr_missed_flip) = [];
    %         end
    %
    %         % Remove remaining flips with mismatching flip intervals
    %         remove_flips = abs(diff(sync_ephys) - ...
    %             diff(sync_timeline(1:length(sync_ephys)))) > missed_flip_t_thresh;
    %         sync_timeline(remove_flips) = [];
    %         sync_ephys(remove_flips) = [];

elseif length(flipper_flip_times_ephys) > length(flipper_times)
    % If more flips in ephys than timeline, screwed
    error([animal ' ' rec_day ': flips ephys > timelite'])
end

% Get spike times in timeline time
spike_times_timeline = interp1(sync_ephys,sync_timeline,spike_times_openephys,'linear','extrap');


%% Load probe position 

% (load histology positions if available)
histology_probe_filename = dir(...
    plab.locations.filename('server',animal,[],[], ...
    'histology','*','probe_ccf.mat'));
if ~isempty(histology_probe_filename)
   probe_histology = load(fullfile(histology_probe_filename.folder, ...
        histology_probe_filename.name));
end

% (load NTE positions if available)
nte_positions_filename = dir(fullfile( ...
    fileparts(open_ephys_path_dir.folder),'*probe_positions*.mat'));
if ~isempty(nte_positions_filename)
    probe_nte = load(fullfile(nte_positions_filename.folder,nte_positions_filename.name));
end

% (use histology if available and matching day, use NTE if not)
if exist('probe_histology','var') && isfield(probe_histology.probe_ccf,'day') && ...
        any(find(strcmp(rec_day,{probe_histology.probe_ccf.day})))
    probe_histology_day_idx = find(strcmp(rec_day,{probe_histology.probe_ccf.day}));
    probe_areas = {probe_histology.probe_ccf(probe_histology_day_idx).trajectory_areas};
    if verbose; disp('Ephys: Loaded histology positions...'); end
else
    probe_areas = probe_nte.probe_areas;
    if verbose; disp('Ephys: Loaded NTE positions...'); end
end

%% Remove bad units from quality control

ephys_qc_type = 'phy';

if strcmp(ephys_qc_type,'bombcell')
    % Load Bombcell quality metrics (if exist)
    qMetrics_path = fullfile(kilosort_path,'qMetrics');
    if  exist(qMetrics_path,'dir')
        % Bombcell function to load all metrics, only used for troubleshooting
        %     [param, qMetric] = bc_loadSavedMetrics(qMetrics_path);
        % (native bombcell labels)
        %  unitType = bc_getQualityUnitType(param, qMetric);
        % (extra metrics & translated bombcell labels)
        load(fullfile(qMetrics_path, 'template_qc_labels.mat'))

        % Define good units from labels
        good_templates = ismember(template_qc_labels,{'singleunit','multiunit'});
        if verbose; disp('Ephys: applying Bombcell quality metrics...'); end
    else
        warning('No ephys quality metrics available');
        good_templates = true(size(templates,1),1);
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
        if verbose; disp('Keeping manually labelled good units...'); end

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
        warning([animal ' ' rec_day ' - no ephys quality control']);
        if verbose; disp('No ephys quality control, keeping all and re-indexing'); end
        good_templates_idx = unique(spike_templates_0idx);
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);
    end
end

% Throw out all non-good template data
templates = templates(good_templates,:,:);
template_depths = template_depths(good_templates);
waveforms = waveforms(good_templates,:);
templateDuration = templateDuration(good_templates);
templateDuration_us = templateDuration_us(good_templates);

% Throw out all non-good spike data
good_spikes = ismember(spike_templates,find(good_templates));
spike_templates = spike_templates(good_spikes);
template_amplitudes = template_amplitudes(good_spikes);
spike_depths = spike_depths(good_spikes);
spike_times_timeline = spike_times_timeline(good_spikes);

% Rename the remaining spike templates (1:N, to match index for template)
[~,spike_templates] = ismember(spike_templates,find(good_templates));





