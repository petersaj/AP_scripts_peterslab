%% SANDBOX
%
% Temporary code

%% Load data

animal = 'AP007';
rec_day = '2023-06-14';

workflow = 'lcr_passive';
% workflow = 'lcr_passive_fullscreen';
% workflow = 'stim_wheel_right*';
% workflow = 'sparse_noise';

rec_time = ap.find_recordings(animal,rec_day,workflow).recording{end};

verbose = true;
ap.load_recording;


%% Testing MCMS API

% MCMS API documentation is here:
% Production: https://oxford.colonymanagement.org/api/swagger-ui/index.html
% Test: https://oxford-uat.colonymanagement.org/api/swagger-ui/index.html

% Get authentication token

% (production)
basicUrl = 'https://oxford.colonymanagement.org/api';
% (test database)
% basicUrl = 'https://oxford-uat.colonymanagement.org/api';
authenticateEndpoint = [basicUrl '/authenticate'];

usr = 'ap7';
psw = 'Bluecookiejar';

headers = struct;
headers.Accept = '*/*';
headers.username = usr;
headers.password = psw;

header_cell = [fieldnames(headers),struct2cell(headers)];

options = weboptions( ...
    'MediaType','application/json', ...
    'ContentType','json', ...
    'RequestMethod','post', ...
    'HeaderFields',header_cell);
mcms_token = webread(authenticateEndpoint,options);


% Get procedure list

proceduresEndpoint = [basicUrl '/procedures'];
headers = struct;
headers.Accept = 'application/json';
headers.Authorization = ['Bearer ' mcms_token.token];

header_cell = [fieldnames(headers),struct2cell(headers)];

options = weboptions( ...
    'MediaType','application/json', ...
    'ContentType','json', ...
    'HeaderFields',header_cell);

data = webread(proceduresEndpoint,options);

% Get weights

curr_animal = '02150140';

endpoint = [basicUrl '/animalweights/animal/' curr_animal];
headers = struct;
headers.Accept = 'application/json';
headers.Authorization = ['Bearer ' mcms_token.token];

header_cell = [fieldnames(headers),struct2cell(headers)];

options = weboptions( ...
    'MediaType','application/json', ...
    'ContentType','json', ...
    'HeaderFields',header_cell);

data = webread(endpoint,options);


data_timestamps = datetime({data.sampleDate},'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSZ','TimeZone','local');

[~,sort_idx] = sort(data_timestamps);
[data(sort_idx).weightValue]


% Get animal via name
curr_animal = 'TOAA2.1d';
endpoint = [basicUrl '/animals/name/' curr_animal];
headers = struct;
headers.Accept = 'application/json';
headers.Authorization = ['Bearer ' mcms_token.token];

header_cell = [fieldnames(headers),struct2cell(headers)];

options = weboptions( ...
    'MediaType','application/json', ...
    'ContentType','json', ...
    'HeaderFields',header_cell);

data = webread(endpoint,options);


% Get project licenses

endpoint = [basicUrl '/projectlicenses'];
headers = struct;
headers.Accept = 'application/json';
headers.Authorization = ['Bearer ' mcms_token.token];

header_cell = [fieldnames(headers),struct2cell(headers)];

options = weboptions( ...
    'MediaType','application/json', ...
    'ContentType','json', ...
    'HeaderFields',header_cell);

data = webread(endpoint,options);

% Cohort history
curr_animal = '2152600';
endpoint = [basicUrl '/animalcohorthistory/animal/' curr_animal];
headers = struct;
headers.Accept = 'application/json';
headers.Authorization = ['Bearer ' mcms_token.token];

header_cell = [fieldnames(headers),struct2cell(headers)];

options = weboptions( ...
    'MediaType','application/json', ...
    'ContentType','json', ...
    'HeaderFields',header_cell);

data = webread(endpoint,options);


%% temp histology: combine jpg channels

hist_path = "P:\Data\AP004\histology";
save_path = "C:\Users\petersa\Desktop\test_histology";
for slice = 1:8
    curr_g_filename = fullfile(hist_path,sprintf('thalamus_%d_g.jpg',slice));
    curr_r_filename = fullfile(hist_path,sprintf('thalamus_%d_r.jpg',slice));

    curr_g = imread(curr_g_filename);
    curr_r = imread(curr_r_filename);

    curr_gr = curr_g + curr_r*3;

    curr_gr_filename = fullfile(save_path,sprintf('%d.tif',slice));
    imwrite(curr_gr,curr_gr_filename);
end

%% temp histology: convert single jpeg to tiff

hist_path = "P:\Data\AP009\histology";
save_path = "C:\Users\petersa\Desktop\test_histology";
for slice = 1:18
    curr_filename = fullfile(hist_path,sprintf('%d.jpg',slice));
    curr_im = imread(curr_filename);

    curr_im = curr_im.*5;

    curr_save_filename = fullfile(save_path,sprintf('%d.tif',slice));
    imwrite(curr_im,curr_save_filename);
end

%% temp: plot probe histology and intended together

nte_filename = "P:\Data\AP004\2023-06-22\ephys\AP004_2023-06-22_probe_positions.mat";
hist_filename = "P:\Data\AP004\histology\slices\probe_ccf.mat";


%% Change server filenames
% currently: move qMetrics folder into pykilosort

server_path = 'P:\Data';
animal_paths = dir(server_path);

for animal_idx = 3:length(animal_paths)

    animal_dir = dir(fullfile(server_path,animal_paths(animal_idx).name));

    date_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
    recording_day_idx = matches({animal_dir.name},date_pattern) & [animal_dir.isdir];

    for day_idx = find(recording_day_idx)

        curr_qmetrics_path = fullfile(animal_dir(day_idx).folder,animal_dir(day_idx).name,'ephys','qMetrics');
        if exist(curr_qmetrics_path,'dir')
            new_qmetrics_path = fullfile(animal_dir(day_idx).folder,animal_dir(day_idx).name,'ephys','pykilosort','qMetrics');

            movefile(curr_qmetrics_path,new_qmetrics_path);
            fprintf('%s --> %s\n',curr_qmetrics_path,new_qmetrics_path);
        end


        %         curr_protocols = dir(fullfile(curr_day_path,'Protocol_*'));
        %
        %         for protocol_idx = 1:length(curr_protocols)
        %             curr_protocol_path = fullfile(curr_protocols(protocol_idx).folder,curr_protocols(protocol_idx).name);
        %             curr_recording_path = strrep(curr_protocol_path,'Protocol','Recording');
        %
        %             movefile(curr_protocol_path,curr_recording_path);
        %             fprintf('%s --> %s\n',curr_protocol_path,curr_recording_path);
        %         end

    end
end
disp('done');

%% Bombcell and supplemental quality metrics

animal = 'AP007';
day = '2023-06-14';

% Get ephys recording paths
% probe_n = multiple probes simultaneously
% site_n = multiple sites recorded in serial
ephys_path = plab.locations.make_server_filename(animal,day,[],'ephys');
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
    fprintf('Bombcell: %s\n',curr_ephys_path);

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
    potential_noise_non_zero = all_unit_non_zero_corr < cutoff_all_unit_non_zero_corr;

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
    potential_noise_max_chan = all_units_corr_max_chan < cutoff_all_units_corr_max_chan;

    % Find number of channels with over-threshold amplitude
    all_units_amp_decay = nan(size(new_raw_waveforms, 1), 1);
    for unit_idx=1:size(new_raw_waveforms,1)
        % find channel with max amplitude
        [raw_max_amplitude, raw_max_channel] = max(max(abs(new_raw_waveforms(unit_idx, :, :)), [], 2), [], 3);
        raw_max_ampl_all_channels = max(abs(new_raw_waveforms(unit_idx, :, :)), [], 2);
        all_units_amp_decay(unit_idx) = sum(raw_max_ampl_all_channels>0.7*raw_max_amplitude);
    end

    cutoff_all_units_amp_decay = 10;
    potential_noise_amp_decay = all_units_amp_decay > cutoff_all_units_amp_decay;

    % Set units that fail extra quality metrics to noise
    potential_noise = potential_noise_max_chan | potential_noise_non_zero | potential_noise_amp_decay;

    am_bc_units = unitType;
    am_bc_units(potential_noise) = 0;

    % Store and save extra quality metrics / classification
    extra_qmetric = struct;
    extra_qmetric.all_unit_non_zero_corr = all_unit_non_zero_corr;
    extra_qmetric.cutoff_all_unit_non_zero_corr = cutoff_all_unit_non_zero_corr;
    extra_qmetric.all_units_corr_max_chan = all_units_corr_max_chan;
    extra_qmetric.cutoff_all_units_corr_max_chan = cutoff_all_units_corr_max_chan;
    extra_qmetric.all_units_amp_decay = all_units_amp_decay;
    extra_qmetric.cutoff_all_units_amp_decay = cutoff_all_units_amp_decay';

    save(fullfile(savePath, 'am_bc_unit_type.mat'), ...
        'unitType', 'am_bc_units','potential_noise', 'extra_qmetric', '-v7.3')

end

disp('Finished.');












