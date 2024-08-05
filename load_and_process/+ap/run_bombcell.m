function run_bombcell(rawFile,ephysKilosortPath,meta_filename,kilosortVersion)
% run_bombcell(ap_band_filename,kilosort_path,meta_filename)
%
% Run Bombcell (JF quality metrics) and supplemental metrics
%
% ap_band_filename - filename of AP-band data used to run kilosort
% kilosort_path - path with kilosort output
% meta_filename - filename for recording metadata

%% Run bombcell 

% Get folders/filenames
disp('Running bombcell...');

% Set save location as the kilosort path
savePath = fullfile(ephysKilosortPath,'qMetrics');

% Get file info of metadata
ephysMetaDir = dir(meta_filename);

% Grab bit-volts from metadata
ephys_metadata = jsondecode(fileread(meta_filename));
gain_to_uV = unique([ephys_metadata.continuous(1).channels.bit_volts]);
if length(gain_to_uV) ~= 1
    error('More than 1 unique bit-volt value');
end

% Load data
[spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, ...
    pcFeatureIdx, channelPositions] = bc.load.loadEphysData(ephysKilosortPath);

% Set parameters (load default, overwrite custom)
param = bc.qm.qualityParamValues(ephysMetaDir, rawFile, ephysKilosortPath, gain_to_uV, kilosortVersion);
param.nChannels = 384;
param.nSyncChannels = 0;
param.extractRaw = 1;
param.plotGlobal = false;
param.plotDetails = false;
% param.firstPeakRatio = 2; %  possibly for future? more sensitive to potential axons

% Run quality metrics
rerun = 0;
qMetricsExist = ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(savePath, 'templates._bc_qMetrics.parquet')));

if qMetricsExist == 0 || rerun
    [qMetric, unitType] = bc.qm.runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
        templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);
else
    [param, qMetric] = bc.load.loadSavedMetrics(savePath); 
    unitType = bc.qm.getQualityUnitType(param, qMetric, savePath);
end

% % (GUI for viewing units/quality metrics)
% % load data for GUI
% loadRawTraces = 0; % default: don't load in raw data (this makes the GUI significantly faster)
% bc.load.loadMetricsForGUI;
% 
% % GUI guide: 
% % left/right arrow: toggle between units 
% % g : go to next good unit 
% % m : go to next multi-unit 
% % n : go to next noise unit
% % up/down arrow: toggle between time chunks in the raw data
% % u: brings up a input dialog to enter the unit you want to go to
% 
% % currently this GUI works best with a screen in portrait mode - we are
% % working to get it to handle screens in landscape mode better. 
% unitQualityGuiHandle = bc.viz.unitQualityGUI(memMapData, ephysData, qMetric, forGUI, rawWaveforms, ...
%     param, probeLocation, unitType, loadRawTraces);


%% Extra quality metric: threshold correlation between template v. average

% Load raw waveforms to run extra template vs. raw quality metrics
% (load channel map to use only channels used by kilosort)
channel_map = readNPY(fullfile(ephysKilosortPath,'channel_map.npy'))+1; % 0-idx > 1-idx
rawWaveforms.average_full = readNPY([fullfile(savePath, 'templates._bc_rawWaveforms.npy')]);
rawWaveforms.average = rawWaveforms.average_full(:,channel_map,:);

% Get full-probe correlation between template and raw waveform
n_templates = size(templateWaveforms,1);
template_waveform_flat = normalize(reshape(permute(templateWaveforms,[2,3,1]),[],n_templates),'center');
raw_waveform_flat = normalize(reshape(permute(rawWaveforms.average,[3,2,1]),[],n_templates),'center');
template_raw_waveform_corr = sum( ...
    (template_waveform_flat./vecnorm(template_waveform_flat,2,1)).* ...
    (raw_waveform_flat./vecnorm(raw_waveform_flat,2,1)),1);

% Set and apply threshold, set bad units in final unitType
template_raw_waveform_corr_thresh = 0.3;
bad_waveform_corr = template_raw_waveform_corr < template_raw_waveform_corr_thresh;

unitType_extra_qc = unitType;
unitType_extra_qc(bad_waveform_corr) = 0;

%% Get labels and save

% Set names for numbered unit types
unit_type_labels = { ...
    0,'noise'; ...
    1,'singleunit';...
    2,'multiunit';...
    3,'axon'};

[~,unitType_idx] = ismember(unitType_extra_qc,cell2mat(unit_type_labels(:,1)));
template_qc_labels = unit_type_labels(unitType_idx,2);

% Save extra quality metrics and final template quality control labels
save(fullfile(savePath, 'template_qc_labels.mat'),'template_qc_labels');

disp('Finished.');





