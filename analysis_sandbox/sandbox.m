%% SANDBOX
%
% dev/test code


%% MCMS API test

mcms_token = plab.mcms.login;

weights = plab.mcms.query(mcms_token,'weight','AP032');



%% Find all animals that have done specific task
% (this takes a long time because it looks through all server folders)

% workflow = '*opacity*';
% workflow = '*no_change*';
% workflow = '*opacity*';
% workflow = '*angle*';
% workflow = 'stim_wheel_right_stage1';
workflow = 'ImageDisplay';

task_dir = dir(fullfile(plab.locations.server_data_path, ...
    '**','bonsai','**',strcat(workflow,'.bonsai')));

animals = unique(extractBetween({task_dir.folder}, ...
    [plab.locations.server_data_path,filesep], ...
    filesep));




%% Move pupil SLEAP files

sleap_dir = dir(fullfile(plab.locations.server_path,'Users','Peter_Gorman','Pupils_temporary','Outputs_Sorted','lcr_passive','**','*.h5'));
sleap_fileparts = regexp({sleap_dir.name},'(?<animal>.*?)_(?<rec_day>.*?)_(?<rec_time>Recording_.*?)_mousecam','names');

for curr_file_idx = 1:length(sleap_dir)  
    curr_fn = fullfile(sleap_dir(curr_file_idx).folder,sleap_dir(curr_file_idx).name);
    target_fn = fullfile(plab.locations.server_data_path, ...
        sleap_fileparts{curr_file_idx}.animal, ...
        sleap_fileparts{curr_file_idx}.rec_day, ...
        sleap_fileparts{curr_file_idx}.rec_time, ...
        'mousecam','sleap','pupil_v1',sleap_dir(curr_file_idx).name);

    mkdir(fileparts(target_fn));
    copyfile(curr_fn,target_fn);
    fprintf('Copied: %s\n',sleap_dir(curr_file_idx).name)
end

%% bombcell axon tests

kilosortVersion = 4;

savePath = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Data\AP022\2024-05-15\ephys\kilosort4\qMetrics';
ephysKilosortPath = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Data\AP022\2024-05-15\ephys\kilosort4';
meta_filename = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Data\AP022\2024-05-15\ephys\experiment1\recording1\structure.oebin';
rawFile = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Data\AP022\2024-05-15\ephys\experiment1\recording1\continuous\Neuropix-PXI-100.ProbeA-AP\continuous.dat';

% savePath = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Data\PG001\2026-04-16\ephys\kilosort4\qMetrics';
% ephysKilosortPath = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Data\PG001\2026-04-16\ephys\kilosort4';
% meta_filename = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Data\PG001\2026-04-16\ephys\experiment1\recording2\structure.oebin';
% rawFile = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Data\PG001\2026-04-16\ephys\experiment1\recording2\continuous\Neuropix-PXI-107.ProbeA-AP\continuous.dat';

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
    pcFeatureIdx, channelPositions] = bc.load.loadEphysData(ephysKilosortPath,savePath);

% Set parameters (load default, overwrite custom)
param = bc.qm.qualityParamValues(ephysMetaDir, rawFile, ephysKilosortPath, gain_to_uV, kilosortVersion);
param.nChannels = 384;
param.nSyncChannels = 0;
param.extractRaw = 1;
param.plotGlobal = false;
param.plotDetails = false;

% Run
rerun = true;
qMetricsExist = ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(savePath, 'templates._bc_qMetrics.parquet')));

% %%%%%%%%%% NEW PARAMS TO TEST FOR AXONS
% param.minWidthFirstPeak_nonSomatic = Inf;
% param.maxPeak1ToPeak2Ratio_nonSomatic = 1/4;
% % % Julie suggestion
% param.maxPeak1ToPeak2Ratio_nonSomatic = 1/3;
% param.minTroughToPeak2Ratio_nonSomatic = 4;
% %%%%%%%%%%

if qMetricsExist == 0 || rerun
    [qMetric, unitType] = bc.qm.runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
        templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);
else
    [param, qMetric] = bc.load.loadSavedMetrics(savePath);
    unitType = bc.qm.getQualityUnitType(param, qMetric, savePath);
end


%%% Extra quality metric: threshold correlation between template v. average

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

%%% Get labels and save

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


%% AP_histology annotations: fix manual alignment 

clear all;

animal = 'DS029';

histology_dir = dir(fullfile(plab.locations.server_data_path,animal,'**','AP_histology_processing.mat'));
histology_fn = fullfile(histology_dir.folder,histology_dir.name);

load(histology_fn);

% Grab atlas images
[av,tv,st] = ap_histology.load_ccf;

n_slices = length(AP_histology_processing.histology_ccf.atlas2histology_tform);
slice_atlas = struct('tv',cell(n_slices,1), 'av',cell(n_slices,1));
slice_atlas_ccf = struct('ap',cell(n_slices,1),'ml',cell(n_slices,1),'dv',cell(n_slices,1));
for curr_slice = 1:n_slices
    [slice_atlas(curr_slice),slice_atlas_ccf(curr_slice)] = ...
        ap_histology.grab_atlas_slice(av,tv, ...
        AP_histology_processing.histology_ccf.slice_vector, ...
        AP_histology_processing.histology_ccf.slice_points(curr_slice,:), 1);
end

for annotation_idx = 1:length(AP_histology_processing.annotation)
    for curr_im_idx = 1:n_slices

        curr_vertices = AP_histology_processing.annotation(annotation_idx).vertices_histology{curr_im_idx};

        if ~any(curr_vertices)
            continue
        end

        % Get transform from slice to atlas
        % (manual if >3 paired control points, automatic otherwise)
        if isfield(AP_histology_processing.histology_ccf,'control_points') && ...
                (size(AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx},1) == ...
                size(AP_histology_processing.histology_ccf.control_points.atlas{curr_im_idx},1)) && ...
                size(AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx},1) >= 3
            % Manual alignment
            atlas_tform = fitgeotform2d( ...
                AP_histology_processing.histology_ccf.control_points.atlas{curr_im_idx}, ...
                AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx},'pwl');
        elseif isfield(AP_histology_processing.histology_ccf,'atlas2histology_tform')
            % Automatic alignment
            atlas_tform = AP_histology_processing.histology_ccf.atlas2histology_tform{curr_im_idx};
        end

        % Convert vertices to CCF coordinates
        annotation_vertices_histology_subscript = round(transformPointsInverse( ...
            atlas_tform,curr_vertices));

        annotation_vertices_histology_idx = ...
            sub2ind(size(slice_atlas(curr_im_idx).tv), ...
            annotation_vertices_histology_subscript(:,2), ...
            annotation_vertices_histology_subscript(:,1));

        [AP_histology_processing.annotation(annotation_idx).vertices_ccf(curr_im_idx).ap, ...
            AP_histology_processing.annotation(annotation_idx).vertices_ccf(curr_im_idx).ml, ...
            AP_histology_processing.annotation(annotation_idx).vertices_ccf(curr_im_idx).dv] = ...
            deal(slice_atlas_ccf(curr_im_idx).ap(annotation_vertices_histology_idx), ...
            slice_atlas_ccf(curr_im_idx).ml(annotation_vertices_histology_idx), ...
            slice_atlas_ccf(curr_im_idx).dv(annotation_vertices_histology_idx));

    end
end

save(histology_fn,'AP_histology_processing');
fprintf('Saved %s\n',histology_fn)

 
% %  (plotting tests)
% figure; hold on;
% 
% curr_probe = 5;
% 
% ap = {AP_histology_processing.annotation(curr_probe).vertices_ccf.ap};
% dv = {AP_histology_processing.annotation(curr_probe).vertices_ccf.dv};
% ml = {AP_histology_processing.annotation(curr_probe).vertices_ccf.ml};
% 
% slice = 10;
% 
% % slice = slice+1;
% plot3(ap{slice},ml{slice},dv{slice},'.k');
% 
% x = vertcat(AP_histology_processing.annotation(curr_probe).vertices_histology{:});




