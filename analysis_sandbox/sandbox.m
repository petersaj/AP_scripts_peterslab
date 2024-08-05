%% SANDBOX
%
% dev/test code

%% Load data (specific day)
clear all

% animal = 'AP019';
% rec_day = '2024-05-06';

% load_parts.widefield = true;
% load_parts.ephys = true;
% load_parts.mousecam = true;

verbose = true;
ap.load_recording;


%% Load data (relative day)

animal = 'AM014';

% workflow = 'lcr_passive';
% workflow = 'lcr_passive_fullscreen';
workflow = 'stim_wheel_right*';
% workflow = 'sparse_noise';
% workflow = 'visual_conditioning*';
% workflow = 'hml_passive_audio';

recordings = plab.find_recordings(animal,[],workflow);

% (include only ephys days)
% recordings = recordings([recordings.ephys]);

use_day = 2;
% use_day = length(recordings);

rec_day = recordings(use_day).day;
rec_time = recordings(use_day).recording{end};

verbose = true;

% load_parts.mousecam = true;
% load_parts.widefield = true;
load_parts.ephys = true;


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

animal = 'AM027';

hist_path = plab.locations.filename('server',animal,[],[],'histology','raw');
save_path = plab.locations.filename('server',animal,[],[],'histology','raw_combined');

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

            fprintf('Wrote %s: %s %s\n',animal,slice_num{curr_slice},im_colors{curr_color});

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

            fprintf('Wrote %s: %s %s\n',animal,slice_num{curr_slice},im_colors{curr_color});
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





%% Load/average/align to create master average blue/violet image

retinotopy_dir = dir(fullfile('\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Andy_Peters\widefield_alignment\retinotopy','*.mat'));

animals = cellfun(@(x) x(end-8:end-4),{retinotopy_dir.name},'uni',false);

avg_im_all = cell(length(animals),2);
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};

    recordings = plab.find_recordings(animal);
    wf_days_idx = cellfun(@(x) any(x),{recordings.widefield});
    wf_recordings = recordings(wf_days_idx);

    avg_im_aligned = cell(length(wf_recordings),2);
    for curr_day = 1:length(wf_recordings)
        try
        day = wf_recordings(curr_day).day;

        img_path = plab.locations.filename('server', ...
            animal,day,[],'widefield');

        avg_im_n = readNPY([img_path filesep 'meanImage_blue.npy']);
        avg_im_h = readNPY([img_path filesep 'meanImage_violet.npy']);

        avg_im_aligned{curr_day,1} = ap.wf_align(avg_im_n,animal,day);
        avg_im_aligned{curr_day,2} = ap.wf_align(avg_im_h,animal,day);
        catch me
            continue
        end
    end

    avg_im_all{curr_animal,1} = nanmean(cat(3,avg_im_aligned{:,1}),3);
    avg_im_all{curr_animal,2} = nanmean(cat(3,avg_im_aligned{:,2}),3);

    AP_print_progress_fraction(curr_animal,length(animals));

end

blue_avg = nanmean(cat(3,avg_im_all{:,1}),3);
violet_avg = nanmean(cat(3,avg_im_all{:,2}),3);



%% Troubleshooting widefield dropped frames

% Get widefield frame times
widefield_metadata_fn = ...
        plab.locations.filename('server',animal,rec_day,[], ...
        'widefield',sprintf('widefield_%s_metadata.bin',rec_time));

fid = fopen(widefield_metadata_fn);
widefield_metadata = reshape(fread(fid,'double'),9,[]);
fclose(fid);
frame_upload_time = datetime(widefield_metadata(4:9,:)');
frame_upload_time_rel = seconds(frame_upload_time - frame_upload_time(1));

% Get times when exposures are done
widefield_exposeOff_times = timelite.timestamps(find(diff(widefield_thresh) == -1) + 1);
widefield_exposeOff_times_rel = widefield_exposeOff_times-widefield_exposeOff_times(1);

% Plot timelite and frame times (relative to first, 1 of every 100);
figure; hold on;
plot(timelite.timestamps - widefield_exposeOff_times(1),timelite.data(:,widefield_idx),'k')
plot(frame_upload_time_rel(1:100),3,'.r');


% Plot frame times - exposure offs
figure;
x = frame_upload_time_rel - widefield_exposeOff_times_rel(1:size(widefield_metadata,2));
subplot(2,1,1);
plot(x,'.k');
ylabel('Relative expose offset to frame time')

subplot(2,1,2);
plot(detrend(x)-min(detrend(x)),'.k');
ylabel('Detrended and min=0')

%% Prepping retinotopy for Jennifer

animal = 'AP014';
rec_day = '2024-01-19';
ap.load_recording;


% Get retinotopy
ap.wf_retinotopy











