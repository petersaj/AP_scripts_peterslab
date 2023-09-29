%% SANDBOX
%
% dev/test code

%% Load data (specific day)

animal = 'AM005';
rec_day = '2023-09-05';

workflow = 'lcr_passive';
% workflow = 'lcr_passive_fullscreen';
% workflow = 'stim_wheel_right*';
% workflow = 'sparse_noise';

rec_time = plab.find_recordings(animal,rec_day,workflow).recording{end};

verbose = true;
ap.load_recording;


%% Load data (relative day)

animal = 'AP011';

workflow = 'lcr_passive';
% workflow = 'lcr_passive_fullscreen';
% workflow = 'stim_wheel_right*';
% workflow = 'sparse_noise';

recordings = plab.find_recordings(animal,[],workflow);

use_day = length(recordings);
rec_day = recordings(use_day).day;
rec_time = recordings(use_day).recording{end};

verbose = true;

% load_parts.ephys = true;

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

%% temp histology: combine tiff channels

hist_path = 'P:\Data\AM005\histology\raw';
save_path = 'P:\Data\AM005\histology';

im_colors = {'w','r'};

im_filenames = dir(fullfile(hist_path,'*.tif'));

slice_names = unique(strtok({im_filenames.name},'_'));

for curr_slice = 1:length(slice_names)
    curr_save_filename = fullfile(save_path,sprintf('%s.tif',slice_names{curr_slice}));
    for curr_color = 1:length(im_colors)
        curr_im_filename = fullfile(hist_path,sprintf('%s_%s.tif', ...
            slice_names{curr_slice},im_colors{curr_color}));
        curr_im = imread(curr_im_filename);
        imwrite(curr_im,curr_save_filename,'tif','WriteMode','append');
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













