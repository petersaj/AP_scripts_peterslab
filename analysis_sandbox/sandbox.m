%% SANDBOX
%
% dev/test code

%% Load data (specific day)

animal = 'AM014';
rec_day = '2024-01-24';

workflow = 'lcr_passive';
% workflow = 'lcr_passive_fullscreen';
% workflow = 'stim_wheel_right*';
% workflow = 'sparse_noise';
% workflow = 'black_screen';
% workflow = 'gray_screen';

rec_time = plab.find_recordings(animal,rec_day,workflow).recording{end};

% load_parts.widefield = true;
load_parts.ephys = true;
% load_parts.mousecam = true;

verbose = true;
ap.load_recording;


%% Load data (relative day)

animal = 'AM017';

workflow = 'lcr_passive';
% workflow = 'lcr_passive_fullscreen';
% workflow = 'stim_wheel_right*';
% workflow = 'sparse_noise';

recordings = plab.find_recordings(animal,[],workflow);

% use_day = 1;
use_day = length(recordings);

rec_day = recordings(use_day).day;
rec_time = recordings(use_day).recording{end};

verbose = true;

% load_parts.behavior = true;
% load_parts.widefield = true;
load_parts.ephys = true;

ap.load_recording;


%% Experiment scroller

ap.expscroll;


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

%% temp histology: combine tiff channels

hist_path = 'P:\Data\AM010\histology\raw';
save_path = 'P:\Data\AM010\histology\raw_combined';

im_filenames = dir(fullfile(hist_path,'*.tif'));
if ~exist(save_path,'dir')
    mkdir(save_path)
end

slice_name_split = regexp({im_filenames.name},'_','split');
slice_name_split_cat = vertcat(slice_name_split{:});

if size(slice_name_split_cat,2) == 2

    slice_num = unique(slice_name_split_cat(:,1));
    im_colors = unique(slice_name_split_cat(:,2));

    for curr_slice = 1:length(slice_num)
        curr_save_filename = fullfile(save_path,sprintf('%s.tif',slice_num{curr_slice}));
        for curr_color = 1:length(im_colors)
            curr_im_filename = fullfile(hist_path,sprintf('%s_%s', ...
                slice_num{curr_slice},im_colors{curr_color}));
            curr_im = imread(curr_im_filename);
            imwrite(curr_im,curr_save_filename,'tif','WriteMode','append');
        end
    end

elseif size(slice_name_split_cat,2) == 3
    
    animal = cell2mat(unique(slice_name_split_cat(:,1)));
    slice_num = unique(slice_name_split_cat(:,2));
    im_colors = unique(slice_name_split_cat(:,3));

    for curr_slice = 1:length(slice_num)
        curr_save_filename = fullfile(save_path,sprintf('%s_%s.tif',animal,slice_num{curr_slice}));
        for curr_color = 1:length(im_colors)
            curr_im_filename = fullfile(hist_path,sprintf('%s_%s_%s', ...
                animal,slice_num{curr_slice},im_colors{curr_color}));
            curr_im = imread(curr_im_filename);
            imwrite(curr_im,curr_save_filename,'tif','WriteMode','append');
        end
    end
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


%% Transform CCF to stereotax for plotting on NTE

% (translation values from our bregma estimate: AP/ML from Paxinos, DV from
% rough MRI estimate)
bregma_ccf = [570.5,520,44]; % [ML,AP,DV]
ccf_translation_tform = eye(4)+[zeros(3,4);-bregma_ccf,0];

% (scaling "Toronto MRI transform", reflect AP/ML, convert 10um to 1mm)
scale = [0.952,-1.031,0.885]./100; % [ML,AP,DV]
ccf_scale_tform = eye(4).*[scale,1]';

% (rotation values from IBL estimate)
ap_rotation = 5; % tilt the CCF 5 degrees nose-up
ccf_rotation_tform = ...
    [1 0 0 0; ...
    0 cosd(ap_rotation) -sind(ap_rotation) 0; ...
    0 sind(ap_rotation) cosd(ap_rotation) 0; ...
    0 0 0 1];

ccf_bregma_tform_matrix = ccf_translation_tform*ccf_scale_tform*ccf_rotation_tform;
ccf_bregma_tform = affine3d(ccf_bregma_tform_matrix);

% (transform points, a = gco);
[ml,ap,dv] = transformPointsForward(ccf_bregma_tform, ...
    a.YData,a.XData,a.ZData);


%% Plot CCF deep areas over widefield image
% (fold into ap.ccf_outline eventually)

% Load atlas
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

% Get area in search, draw

% Prompt for which structures to show (only structures which are
% labelled in the slice-spacing downsampled annotated volume)
structure_search = lower(inputdlg('Search structures'));
structure_match = find(contains(lower(st.safe_name),structure_search));

selected_structure = listdlg('PromptString','Select a structure to plot:', ...
    'ListString',st.safe_name(structure_match),'ListSize',[520,500], ...
    'SelectionMode','single');

plot_structure = structure_match(selected_structure);

% Get all areas within and below the selected hierarchy level
plot_structure_id = st.structure_id_path{plot_structure};
plot_ccf_idx = find(cellfun(@(x) contains(x,plot_structure_id), ...
    st.structure_id_path));

% Get structure color and volume
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{plot_structure},2,[])')./255;
plot_ccf_volume = ismember(av,plot_ccf_idx);

% Get top-down structure boundaries
structure_outline = bwboundaries(squeeze((max(plot_ccf_volume,[],2))));

% Transform structure coordinates
um2pixel = 20.6; % (hardcoded for now - is this necessary?)

alignment_path = fullfile(plab.locations.server_path, ...
    'Users','Andy_Peters','widefield_alignment');
ccf_tform_fn = [alignment_path filesep 'ccf_tform.mat'];
load(ccf_tform_fn);

structure_outline_aligned_x = cellfun(@(coords) ...
    [fliplr(coords*(10/um2pixel)),ones(size(coords,1),1)]*ccf_tform.T, ...
    structure_outline,'uni',false);
structure_outline_aligned = cellfun(@(coords) ...
    coords(:,[2,1]),structure_outline_aligned_x,'uni',false);

% Plot 2D structure
ccf_axes = gca;
cellfun(@(x) plot(ccf_axes(curr_view),x(:,2), ...
        x(:,1),'color',plot_structure_color,'linewidth',2), ...
        structure_outline_aligned)


%% Neuropixels preprocessing

ap.preprocess_neuropixels('AM014','2024-01-25');

%% Load raw data snippet

% animal = 'AM016';
% day = '2024-02-08';
% 
% ap_dir = dir(plab.locations.filename('server',animal,day,[],'ephys','**','*-AP','continuous.dat'));
% ap_filename = fullfile(ap_dir.folder,ap_dir.name);
% fid = fopen(ap_filename, 'r');
% fseek(fid,(30000*60*20*384),'bof');
% dat = fread(fid, [384 chunkSize], '*int16');

ap_filename = 'D:\data_temp\kilosort\continuous.dat';
ap_car_filename = 'D:\data_temp\kilosort\continuous_car.dat';

% Choose data chunk size
chunkSize = 1000000;

n_chan = 384;

ap_filename_dir = dir(ap_filename);
nSampsTotal = ap_filename_dir.bytes/n_chan/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);

% Get channel index across ADCs (16 ADCs, 12 channels each, same-index
% channels sampled together e.g. ADC1/Ch1 + ADC2/Ch1...)
% (https://open-ephys.github.io/gui-docs/User-Manual/Plugins/Neuropixels-CAR.html)
ch_adc_idx = reshape(reshape(1:384,2,[])',12,[]);

% Open original file for reading, and new file for writing
fid_raw = fopen(ap_filename, 'r');
fid_car = fopen(ap_car_filename, 'w');

% Loop through data chunks and apply common average referencing
chunkInd = 1;
while 1

    % Load in current data chunk
    curr_data = fread(fid_raw, [384 chunkSize], '*int16');

    if ~isempty(curr_data)

        % Subtract median across channels
        curr_data_centered = curr_data - median(curr_data,2);

        % Loop through ADC channel indicies
        curr_data_car = zeros(size(curr_data),'int16');
        for curr_adc_idx = 1:size(ch_adc_idx,1)
            % Get same-index channels across ADCs
            curr_channels = ch_adc_idx(curr_adc_idx,:);

            % Get median signal across channels
            curr_data_adc_median = median(curr_data_centered(curr_channels,:),1);

            % Get scaling factor of median for each channel
            curr_data_adc_median_scale = ...
                sum(curr_data_centered(curr_channels,:).*curr_data_adc_median,2)./ ...
                sum(curr_data_adc_median.^2);
            curr_data_adc_median_scaled = int16(curr_data_adc_median_scale.* ...
                double(curr_data_adc_median));

            % Subtract scaled median from data
            curr_data_car(curr_channels,:) = curr_data_centered(curr_channels,:) - ...
                curr_data_adc_median_scaled;
        end

        fwrite(fid_car, curr_data_car, 'int16');

    else
        break
    end

    AP_print_progress_fraction(chunkInd,nChunksTotal);
    chunkInd = chunkInd+1;

end

% Close files
fclose(fid_raw);
fclose(fid_car);

disp('Finished CAR');


        

%% Testing Kilosort 2

ap_sample_rate = 30000;
ap_clean_filename = 'D:\data_temp\kilosort\continuous_car.dat';

ssd_kilosort_path = 'D:\data_temp\kilosort\kilosort2';
t_range = [0,Inf];

AP_run_kilosort2(ap_clean_filename,ap_sample_rate,ssd_kilosort_path,t_range);

% errored "EIG did not converge", probably needs downgrade: 
% https://github.com/cortex-lab/KiloSort/issues/260

%% Testing pykilosort

local_kilosort_path = 'D:\data_temp\kilosort\';
apband_local_filename = fullfile(local_kilosort_path,'continuous_car.dat');

% Path to pykilosort python environment
pykilosort_python_environment = pyenv('Version','C:\Users\petersa\anaconda3\envs\pyks2\pythonw.exe');

% Set up python
% (add pykilosort environment paths to windows system path)
pre_pykilosort_syspath = getenv('PATH');
py_env_paths = {
    fullfile(char(pykilosort_python_environment.Home),'Library','bin'); ...
    fullfile(char(pykilosort_python_environment.Home),'Scripts')};
run_pykilosort_syspath = strjoin( ...
    unique(cat(1,py_env_paths, ...
    split(pre_pykilosort_syspath,pathsep)),'stable'),pathsep);
setenv('PATH',run_pykilosort_syspath);

% Run pykilosort (does common average referencing by default)
pykilosort_output_path = fullfile(local_kilosort_path,'pykilosort');
pyrunfile('AP_run_pykilosort.py', ...
    data_filename = apband_local_filename, ...
    pykilosort_output_path = pykilosort_output_path);

% Revert system paths to pre-pykilosort
% (just in case alternate python environments used elsewhere)
setenv('PATH',pre_pykilosort_syspath);






