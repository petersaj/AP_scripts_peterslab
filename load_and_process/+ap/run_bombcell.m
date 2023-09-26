function run_bombcell(animal,day)
% run_bombcell(animal,day)
%
% Run Bombcell (JF quality metrics) and supplemental metrics

% Get ephys recording paths
% probe_n = multiple probes simultaneously
% site_n = multiple sites recorded in serial
ephys_path = plab.locations.filename('server',animal,day,[],'ephys');
ephys_path_dir = ...
    [dir(fullfile(ephys_path, 'probe_*')), ...
    dir(fullfile(ephys_path, 'site_*'))];
if ~isempty(ephys_path_dir)
    data_paths = cellfun(@(x) [ephys_path filesep x],{ephys_path_dir.name},'uni',false);
else
    data_paths = {ephys_path};
end

for curr_data = 1:length(data_paths)

    % Get folders/filenames
    curr_ephys_path = data_paths{curr_data};
    fprintf('Running bombcell: %s\n',curr_ephys_path);

    % (recording path as experimentX/recordingX)
    curr_ephys_rec_dir = dir(fullfile(curr_ephys_path, 'experiment*', 'recording*'));
    curr_ephys_rec_path = fullfile(curr_ephys_rec_dir.folder,curr_ephys_rec_dir.name);
    % (oebin metadata filename)
    ephys_meta_dir = dir(fullfile(curr_ephys_rec_path,'*.oebin'));
    % (AP-band data filename)
    ap_data_filename = fullfile(curr_ephys_rec_path, ...
        'continuous', 'Neuropix-3a-100.Neuropix-3a-AP', 'continuous.dat');
    % (kilosort folder)
    curr_kilosort_path = strrep(curr_ephys_path,'ephys',fullfile('ephys','pykilosort'));

    % Set save location as the kilosort path
    savePath = fullfile(curr_kilosort_path,'qMetrics');

    % Set parameters (load default, overwrite custom)
    param = bc_qualityParamValues(ephys_meta_dir, ap_data_filename);
    param.nChannels = 384;
    param.nSyncChannels = 0;
    param.extractRaw = 1;
    param.plotGlobal = false;
    param.plotDetails = false;

    % Load data
    [spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, ...
        pcFeatureIdx, channelPositions] = bc_loadEphysData(curr_kilosort_path);

    % Run quality metrics
    [qMetric, unitType] = bc_runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
        templateWaveforms, templateAmplitudes,pcFeatures,pcFeatureIdx,channelPositions, savePath);

    % Load raw waveforms to run extra template vs. raw quality metrics
    rawWaveforms.average = readNPY([fullfile(savePath, 'templates._bc_rawWaveforms.npy')]);
    rawWaveforms.peakChan = readNPY([fullfile(savePath, 'templates._bc_rawWaveformPeakChannels.npy')]);

    % Correlate template/raw waveforms concatenated across non-zero channels
    new_raw_waveforms = permute(rawWaveforms.average, [1 3 2]);
    all_unit_non_zero_corr = nan(size(new_raw_waveforms,1), 1);
    for unit_idx=1:size(new_raw_waveforms, 1)

        % get non-zero channels in template waveforms
        temp_template_waveforms = squeeze(templateWaveforms(unit_idx,:,:));

        [template_max_amplitude, template_max_channel] = max(max(abs(templateWaveforms(unit_idx, :, :)), [], 2), [], 3);

        zero_channel_idx = arrayfun(@(X) all(temp_template_waveforms(:,X) < 0.05* template_max_amplitude...
            & temp_template_waveforms(:,X) > -0.05* template_max_amplitude), [1:param.nChannels]);

        non_zero_template_waveforms = temp_template_waveforms(:, ~zero_channel_idx);

        % concatenate for comparison
        vector_non_zero_template_waveforms = non_zero_template_waveforms(:);

        % get raw waveforms on the same channels
        temp_raw_waveforms = squeeze(new_raw_waveforms(unit_idx,:,:));
        non_zero_raw_waveforms = temp_raw_waveforms(:, ~zero_channel_idx);
        vector_non_zero_raw_waveforms = non_zero_raw_waveforms(:);

        all_unit_non_zero_corr(unit_idx) = corr(vector_non_zero_raw_waveforms, vector_non_zero_template_waveforms);
    end

    cutoff_all_unit_non_zero_corr = 0.6;
    noise_non_zero = all_unit_non_zero_corr < cutoff_all_unit_non_zero_corr;

    % Correlate template/raw waveforms on respective channel with max amplitude
    all_units_corr_max_chan = nan(size(new_raw_waveforms, 1), 1);
    for unit_idx=1:size(templateWaveforms,1)
        % find channel with max amplitude
        [template_max_amplitude, template_max_channel] = max(max(abs(templateWaveforms(unit_idx, :, :)), [], 2), [], 3);
        [raw_max_amplitude, raw_max_channel] = max(max(abs(new_raw_waveforms(unit_idx, :, :)), [], 2), [], 3);

        % get corr with max raw channel
        raw_wvf = squeeze(new_raw_waveforms(unit_idx, :, raw_max_channel));
        template_wvf = squeeze(templateWaveforms(unit_idx, :, template_max_channel));
        all_units_corr_max_chan(unit_idx) = corr(raw_wvf', template_wvf');
    end

    cutoff_all_units_corr_max_chan = 0;
    noise_max_chan = all_units_corr_max_chan < cutoff_all_units_corr_max_chan;

    % Find number of channels with over-threshold amplitude
    all_units_amp_decay = nan(size(new_raw_waveforms, 1), 1);
    for unit_idx=1:size(new_raw_waveforms,1)
        % find channel with max amplitude
        [raw_max_amplitude, raw_max_channel] = max(max(abs(new_raw_waveforms(unit_idx, :, :)), [], 2), [], 3);
        raw_max_ampl_all_channels = max(abs(new_raw_waveforms(unit_idx, :, :)), [], 2);
        all_units_amp_decay(unit_idx) = sum(raw_max_ampl_all_channels>0.7*raw_max_amplitude);
    end

    cutoff_all_units_amp_decay = 10;
    noise_amp_decay = all_units_amp_decay > cutoff_all_units_amp_decay;

    % Store and save extra quality metrics / classification
    extra_qmetric = struct;
    extra_qmetric.all_unit_non_zero_corr = all_unit_non_zero_corr;
    extra_qmetric.cutoff_all_unit_non_zero_corr = cutoff_all_unit_non_zero_corr;
    extra_qmetric.all_units_corr_max_chan = all_units_corr_max_chan;
    extra_qmetric.cutoff_all_units_corr_max_chan = cutoff_all_units_corr_max_chan;
    extra_qmetric.all_units_amp_decay = all_units_amp_decay;
    extra_qmetric.cutoff_all_units_amp_decay = cutoff_all_units_amp_decay';

    % Set units that fail extra quality metrics to noise
    extra_qmetric_noise = noise_max_chan | noise_non_zero | noise_amp_decay;

    unitType_extra_qc = unitType;
    unitType_extra_qc(extra_qmetric_noise) = 0;

    % Set names for numbered unit types
    unit_type_labels = { ...
        0,'noise'; ...
        1,'singleunit';...
        2,'multiunit';...
        3,'axon'};

    [~,unitType_idx] = ismember(unitType_extra_qc,cell2mat(unit_type_labels(:,1)));
    template_qc_labels = unit_type_labels(unitType_idx,2);

    % Save extra quality metrics and final template quality control labels
    save(fullfile(savePath, 'extra_qmetric.mat'),'extra_qmetric');
    save(fullfile(savePath, 'template_qc_labels.mat'),'template_qc_labels');

end

disp('Finished.');




