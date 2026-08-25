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

% Get site if serial multi-site
% (multi-site with same probe = /site_n - not common/used anymore)
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
% (get tip distance center-of-mass on thresholded channel amplitudes)
template_tipdist = sum(template_chan_amp_overthresh.*channel_positions(:,2)',2)./sum(template_chan_amp_overthresh,2);
% (get shank for each unit)
shank_spacing = 250;
shank_borders = (0:4)*shank_spacing-shank_spacing/2;
template_xpos = sum(template_chan_amp_overthresh.*channel_positions(:,1)',2)./sum(template_chan_amp_overthresh,2);
template_shanks = discretize(template_xpos,shank_borders);

% Get the depth of each spike
spike_tipdist = template_tipdist(spike_templates);

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

% (set last update of function to check for current version)
ks2oe_timestamps_lastupdate = datetime('2026-06-25');

if ~(exist(spike_times_openephys_filename,'file') && exist(open_ephys_ttl_filename,'file')) || ...
        (datetime(dir(spike_times_openephys_filename).date) < ks2oe_timestamps_lastupdate)
    % If timestamps don't exist or are before last update: (re)create them
    ks_spike_times_fn = fullfile(kilosort_path,'spike_times.npy');

    oe_samples_dir = cellfun(@(data_path) ...
        dir(fullfile(data_path,'sample_numbers.npy')), ...
        open_ephys_path,'uni',false);
    oe_samples_fns = cellfun(@(data_dir) ...
        fullfile(data_dir.folder,data_dir.name),oe_samples_dir,'uni',false);

    plab.ephys.ks2oe_timestamps(ks_spike_times_fn,oe_samples_fns, ...
        oe_metadata(1).continuous(1).sample_rate,true);
end

% Load timestamps
spike_times_openephys = readNPY(spike_times_openephys_filename);
load(open_ephys_ttl_filename);

%% Convert ephys to timelite times

% (allow for multiple files if multiple recordings concatenated)

% Set TTL index for flipper (sync)
flipper_sync_idx = 1;
open_ephys_ttl_flipper_idx = find(abs(open_ephys_ttl.state) == flipper_sync_idx);

% Resample Open Ephys flipper to DAQ sample rate
open_ephys_flipper_trace_t = ...
    open_ephys_ttl.timestamps(open_ephys_ttl_flipper_idx(1)): ...
    1/timelite.daq_info(1).rate: ...
    open_ephys_ttl.timestamps(open_ephys_ttl_flipper_idx(end));
open_ephys_flipper_trace = logical(normalize(interp1( ...
    open_ephys_ttl.timestamps(open_ephys_ttl_flipper_idx), ...
    single(sign(open_ephys_ttl.state(open_ephys_ttl_flipper_idx))), ...
    open_ephys_flipper_trace_t,'previous'),'range'));

% Get Open Ephys corresponding to timelite flipper
% (get lag between timelite and ephys by correlation)
ephys_timelite_flipper_lag = finddelay(+flipper_thresh, ...
    +open_ephys_flipper_trace)/timelite.daq_info(1).rate + ...
    open_ephys_flipper_trace_t(1);
% (get all ephys flips within the matching continuous window)
curr_ephys_flipper_idx = intersect(...
    open_ephys_ttl_flipper_idx, ... % the flipper channel
    find(isbetween(open_ephys_ttl.timestamps,ephys_timelite_flipper_lag, ... % during timelite time frame
    ephys_timelite_flipper_lag + length(flipper_thresh)/timelite.daq_info(1).rate)));

% (flag very short ephys flips to remove -occsaionally from bad connection)
short_ephys_flips = reshape(open_ephys_ttl_flipper_idx(find(diff( ...
    open_ephys_ttl.timestamps(open_ephys_ttl_flipper_idx)) < 1e-4) + [0,1]),[],1);
% (set ephys flips that should correspond to timelite)
flipper_times_ephys = ...
    open_ephys_ttl.timestamps(setdiff(curr_ephys_flipper_idx,short_ephys_flips));

% % (for debugging: plot aligned timelite/ephys flipper)
% figure; hold on;
% plot((timelite.timestamps-timelite.timestamps(1))+ephys_timelite_flipper_lag,flipper_thresh);
% plot(open_ephys_flipper_trace_t,open_ephys_flipper_trace+1.2);
% plot(open_ephys_ttl.timestamps(short_ephys_flips),2.2,'.r','MarkerSize',20);
% ylim([-1,2.5]);
% legend({'Timelite','Ephys'})

% Pick flipper times to use for alignment
if length(flipper_times_ephys) == length(flipper_times)
    % If same number of flips in ephys/timelite, use all
    sync_timelite = flipper_times;
    sync_ephys = flipper_times_ephys;

else
    % If different number of timelite/ephys flips: 
    % find usable flips as "matched" with a very short estimated time delay   

    % Estimate nearest flip for timelite->ephys and ephys->timelite
    flip_timelite2ephys_timediff = ...
        interp1(flipper_times_ephys,flipper_times_ephys, ...
        flipper_times + ephys_timelite_flipper_lag,'nearest','extrap') - ...
        (flipper_times+ephys_timelite_flipper_lag);

    flip_ephys2timelite_timediff = ...
        interp1((flipper_times+ephys_timelite_flipper_lag),(flipper_times+ephys_timelite_flipper_lag), ...
        flipper_times_ephys,'nearest','extrap') - ...
        flipper_times_ephys;

    % Set cutoff to find "matched" flips
    % (remove rlowess-smoothed to account for clock drift: annoyingly slow,
    % but the most robust so far)
    % (for ephys: also remove very quick flips e.g. bad connection)
    flip_timediff_thresh = 1e-3; % empirical
    use_timelite_flips = abs(flip_timelite2ephys_timediff-smooth(flip_timelite2ephys_timediff,1000,'rlowess')) < flip_timediff_thresh;
    use_ephys_flips = abs(flip_ephys2timelite_timediff-smooth(flip_ephys2timelite_timediff,1000,'rlowess')) < flip_timediff_thresh & ...
            vertcat(true,diff(flipper_times_ephys) > 0.05);

    if sum(use_timelite_flips) == sum(use_ephys_flips)
        % Successful matching if same number of usable flips
        warning('%s %s: Unmatched timelite/ephys flips - found matches',animal,rec_day);
        sync_ephys = flipper_times_ephys(use_ephys_flips);
        sync_timelite = flipper_times(use_timelite_flips);
    else
        % If not, unclear where the matched flips are, error out
        error('%s %s: Unmatched timelite/ephys flips - cannot match',animal,rec_day);
    end

end

% Get spike times in timelite time
spike_times_timelite = interp1(sync_ephys,sync_timelite,spike_times_openephys,'linear','extrap');

%% Load probe position 

% Load histology positions (if available)
histology_dir = dir(plab.locations.filename('server',animal,[],[], ...
    'histology','**','AP_histology_processing.mat'));
if ~isempty(histology_dir)
    histology_filename = fullfile(histology_dir.folder,histology_dir.name);
    load(histology_filename);
    
    % (check if annotation ephys path matches loaded path)
    if isfield(AP_histology_processing,'annotation') && ...
            isfield(AP_histology_processing.annotation,'ephys_path')

        % (get annotation for loaded recording)
        histology_annotation_match_unsorted = find(strcmp(kilosort_path, ...
            {AP_histology_processing.annotation.ephys_path}));
        % (sort by shank index)
        [~,histology_annotation_shanksort_idx] = ...
            sort([AP_histology_processing.annotation(histology_annotation_match_unsorted).ephys_shank]);
        histology_annotation_match = histology_annotation_match_unsorted(histology_annotation_shanksort_idx);

        % (get probe points in CCF)
        probe_vector_histology = cat(3, ...
            ap_histology.fit_probe_line(histology_filename, ...
            histology_annotation_match).ccf);

        % Check for adjusted areas
        if isfield(AP_histology_processing.annotation,'probe_areas')
            probe_histology_adjusted = ...
                {AP_histology_processing.annotation(histology_annotation_match).probe_areas};
        end
        
        % Set histology areas
        if exist('probe_histology_adjusted','var') && ...
                all(~cellfun(@isempty,probe_histology_adjusted))
            % If all areas adjusted, use those as-is
            probe_histology = vertcat(probe_histology_adjusted{:});
        else
            % Grab probe areas from histology
            probe_histology = plab.histology.grab_probe_areas(probe_vector_histology);
            if exist('probe_histology_adjusted','var')
                % Adjust areas by shank, if present
                for adjusted_shank = find(~cellfun(@isempty,probe_histology_adjusted))
                    probe_histology(probe_histology.probe_shank == adjusted_shank,:) = ...
                        probe_histology_adjusted{adjusted_shank};
                end
            end
        end
    end
end

% Load NTE positions (if available)
nte_positions_filename = dir(fullfile(erase(kilosort_path,[filesep,'kilosort4']),'**','*probe_positions*.mat'));
if ~isempty(nte_positions_filename)
    probe_nte = load(fullfile(nte_positions_filename.folder,nte_positions_filename.name));
end

% Set probe areas (adjusted > histology > NTE)
if exist('probe_histology','var')
    probe_areas = probe_histology;
    if verbose; disp('Ephys: Loaded histology positions...'); end
elseif exist('probe_nte','var')
    probe_areas = probe_nte.probe_areas{load_probe};
    if verbose; disp('Ephys: Loaded NTE positions...'); end
end

% NTE LEGACY SUPPORT
% - If `probe_depth` instead of `tip_distance`, calculate
if ~any(strcmp(probe_areas.Properties.VariableNames,'tip_distance')) && ...
        any(strcmp(probe_areas.Properties.VariableNames,'probe_depth'))
    probe_areas.tip_distance = 3840 - probe_areas.probe_depth;
end
% - if no `probe_shank`, add 1's
if ~any(strcmp(probe_areas.Properties.VariableNames,'probe_shank'))
    probe_areas.probe_shank = ones(height(probe_areas),1);
end


%% Get CCF positions of units 

% Loop through shanks, interpolate template position relative to shank in
% CCF space: [AP,DV,ML]
% (assumes CCF = real distance, which is only approximate)

if exist('probe_areas','var') && ...
        any(ismember(probe_areas.Properties.VariableNames,'ccf'))
    % From histology, if available

    ccf2um = 10; % conversion factor: CCF is in 10um voxels (untransformed)     
    template_ccf = nan(size(templates,1),3);

    % (loop through shanks)
    for curr_shank = unique(template_shanks)'

        use_area_idx = ...
            probe_areas.probe_shank == curr_shank & ... & on current shank
            -diff(probe_areas.tip_distance,[],2) > 0; % area size is > 0

        probe_setpoints_tipdist = 1000 * ... % (convert mm to um)
            vertcat(probe_areas.tip_distance(use_area_idx,1), ...
            probe_areas.tip_distance(end,2));

        probe_setpoints_ccf = cell2mat( ...
            vertcat(probe_areas.ccf(use_area_idx,1), ...
            probe_areas.ccf(end,2)));

        template_ccf(template_shanks == curr_shank,:) = ...
            interp1(probe_setpoints_tipdist,probe_setpoints_ccf, ...
            template_tipdist(template_shanks==curr_shank));
    end
elseif exist('probe_nte','var')
    % From NTE otherwise
    % (saves positions in as probe_nte.probe_positions_ccf [top1,bottom1,top2...])

    template_ccf = nan(size(templates,1),3);
    for curr_shank = unique(template_shanks)'
        curr_shank_vector = probe_nte.probe_positions_ccf{load_probe}(:,2*curr_shank+[-1,0]);
        template_ccf(template_shanks==curr_shank,:) = ...
            interp1([10*norm(diff(curr_shank_vector,[],2)),0],curr_shank_vector', ...
            template_tipdist(template_shanks==curr_shank));
    end
end

% % % (DEBUG) plot units in CCF space
% ap.ccf_outline_3d([],{'brain','cp','aca','snr'});
% scatter3(template_ccf(:,1),template_ccf(:,3),template_ccf(:,2),10,'k','filled');


%% Quality control units (bombcell)

% Load Bombcell quality metrics
qMetrics_path = fullfile(kilosort_path,'qMetrics');
if  exist(qMetrics_path,'dir')
    % (keep axons if flagged)
    if isfield(load_parts,'ephys_axons') && ...
            load_parts.ephys_axons
        bombcell_keep = {'singleunit','multiunit','axon'};
        if verbose; fprintf('Ephys: Keeping axons...\n'); end
    else
        bombcell_keep = {'singleunit','multiunit'};
    end

    % Load unit labels (from ap.run_bombcell / ap.rerun_bombcell)
    load(fullfile(qMetrics_path, 'template_qc_labels.mat'))

    % Define good units from labels
    good_templates = ismember(template_qc_labels,bombcell_keep);

    % Keep only labels from good units
    template_qc_labels = template_qc_labels(good_templates);
else
    warning('Bombcell metrics not available');
    return
end

% If good templates were selected above, throw out not-good data:

% Throw out all non-good template data
templates = templates(good_templates,:,:);
template_tipdist = template_tipdist(good_templates);
template_xpos = template_xpos(good_templates);
template_shanks = template_shanks(good_templates);
waveforms = waveforms(good_templates,:);
waveform_duration_peaktrough = waveform_duration_peaktrough(good_templates);
waveform_duration_fwhm = waveform_duration_fwhm(good_templates);
template_ccf = template_ccf(good_templates,:);

% Throw out all non-good spike data
good_spikes = ismember(spike_templates,find(good_templates));
spike_templates = spike_templates(good_spikes);
template_amplitudes = template_amplitudes(good_spikes);
spike_tipdist = spike_tipdist(good_spikes);
spike_times_timelite = spike_times_timelite(good_spikes);

% Rename the remaining spike templates (1:N, to match index for template)
[~,spike_templates] = ismember(spike_templates,find(good_templates));

if verbose; fprintf('Ephys: Bombcell kept %d/%d units...\n',sum(good_templates),length(good_templates)); end
















