% Load electrophysiology data

ephys_quality_control = false;

ephys_path = plab.locations.make_server_filename(animal,rec_day,[],'ephys');

if verbose; disp('Loading Ephys...'); end

% Get path for raw data and kilosort
% TO DO: if multiple sites (in serial), grab closest to protocol

kilosort_top_path = fullfile(ephys_path,'pykilosort');
kilosort_dir = dir(kilosort_top_path);

if ~any([kilosort_dir.isdir] & ~contains({kilosort_dir.name},'.'))
    % If no kilosort subfolders, use that top path
    kilosort_path = kilosort_top_path;
    open_ephys_path_dir = dir(fullfile(ephys_path,'experiment*','recording*'));
    open_ephys_path = fullfile(open_ephys_path_dir.folder,open_ephys_path_dir.name);

elseif any(contains({kilosort_dir.name},'site'))
    % If multiple sites (in serial), choose last before protocol
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

        ephys_use_site = find(ephys_site_datetime - rec_datetime < 0,1,'last');

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

% Load probe position from trajectory explorer, if exists
probe_positions_filename = dir(fullfile(fileparts(open_ephys_path_dir.folder),'*probe_positions*.mat'));
if ~isempty(probe_positions_filename)
    probe_positions = load(fullfile(probe_positions_filename.folder,probe_positions_filename.name));
end

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

% Load ephys flipper
flipper_sync_idx = 1;
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

% (spike times: load open ephys times if available, create if not)
spike_times_openephys_filename = fullfile(kilosort_path,'spike_times_openephys.npy');
if exist(spike_times_openephys_filename,'file')
    spike_times_openephys = readNPY(spike_times_openephys_filename);
else
    open_ephys_timestamps = readNPY(fullfile(open_ephys_path, ...
        'continuous','Neuropix-3a-100.Neuropix-3a-AP','timestamps.npy'));
    spike_times_openephys = open_ephys_timestamps( ...
        readNPY(fullfile(kilosort_path,'spike_times.npy')));
    writeNPY(spike_times_openephys,spike_times_openephys_filename);
    fprintf('Saved spike times openephys: %s\n',spike_times_openephys_filename);
end

% (spike times: index Open Ephys timestamps rather than assume constant
% sampling rate as before, this accounts for potentially dropped data)
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

% Get flipper experiment boundaries by long delays
flip_diff_thresh = 5; % time between flips to define experiment gap (s)
flipper_protocol_idx = [1;find(abs(diff(open_ephys_flipper.timestamps)) > ...
    flip_diff_thresh)+1;length(open_ephys_flipper.timestamps)+1];

% Get protocols included in recording
% (first protocol after recording starts until N recorded protocols)
day_protocols = ap.find_recordings(animal,rec_day).protocol;
protocol_datetimes = ...
    cellfun(@(x) datetime(strjoin({rec_day,x}),'InputFormat','yyyy-MM-dd HHmm'), ...
    day_protocols);
% (add 1 minute leeway on recording time since seconds not included)
ephys_included_protocols = find(ephys_datetime -  ...
    (protocol_datetimes + minutes(1)) < 0, ...
    length(flipper_protocol_idx)-1,'first');
ephys_protocol_idx = strcmp(day_protocols(ephys_included_protocols),rec_time);

flipper_flip_times_ephys = open_ephys_flipper.timestamps( ...
    flipper_protocol_idx(ephys_protocol_idx):flipper_protocol_idx(find(ephys_protocol_idx)+1)-1);

% Pick flipper times to use for alignment
if length(flipper_flip_times_ephys) == length(flipper_times)
    % If same number of flips in ephys/timeline, use all
    sync_timeline = flipper_times;
    sync_ephys = flipper_flip_times_ephys;

elseif length(flipper_flip_times_ephys) < length(flipper_times)
    warning([animal ' ' rec_day ': ephys missed flips, using first/last - better fix IN PROGRESS'])
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
spike_templates_0idx = spike_templates_0idx(good_spike_idx);
template_amplitudes = template_amplitudes(good_spike_idx);
spike_depths = spike_depths(good_spike_idx);
spike_times_timeline = spike_times_timeline(good_spike_idx);

% Rename the spike templates according to the remaining templates
% (and make 1-indexed from 0-indexed)
new_spike_idx = nan(max(spike_templates_0idx)+1,1);
new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
spike_templates = new_spike_idx(spike_templates_0idx+1);



