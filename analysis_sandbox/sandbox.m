%% SANDBOX
% 
% Temporary code

%% Align widefield to event

% align_times = photodiode_times(1:2:end);
% align_category = vertcat(trial_events.values.TrialStimX);

align_times = reward_times;
align_category = ones(size(align_times));

surround_window = [-0.2,4];

surround_samplerate = 35;
t = surround_window(1):1/surround_samplerate:surround_window(2);
peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

use_U = wf_U;
use_V = wf_V;
use_wf_t = wf_times;

aligned_v = reshape(interp1(use_wf_t,use_V',peri_event_t,'previous'), ...
    length(align_times),length(t),[]);

align_id = findgroups(align_category);
aligned_v_avg = permute(splitapply(@nanmean,aligned_v,align_id),[3,2,1]);
aligned_v_avg_baselined = aligned_v_avg - nanmean(aligned_v_avg(:,t < 0,:),2);

aligned_px_avg = plab.wf.svd2px(use_U,aligned_v_avg_baselined);

AP_imscroll(aligned_px_avg,t);
colormap(AP_colormap('PWG'));
clim(prctile(abs(aligned_px_avg(:)),99.5).*[-1,1]);
axis image;

%% Testing MCMS API

% MCMS API documentation is here: 
% https://oxford-uat.colonymanagement.org/api/swagger-ui/index.html

% Get authentication token 

basicUrl = 'https://oxford-uat.colonymanagement.org/api';
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


% Get name

curr_animal = 'TOAA1.1a';
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
















