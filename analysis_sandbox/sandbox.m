%% SANDBOX
%
% dev/test code

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

animal = 'AP026';

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

%% Find all animals that have done specific task
% (this takes a long time because it looks through all server folders)

% workflow = '*opacity*';
% workflow = '*no_change*';
% workflow = '*opacity*';
% workflow = '*angle*';
workflow = 'stim_wheel_right_stage1';

task_dir = dir(fullfile(plab.locations.server_data_path, ...
    '**','bonsai','**',strcat(workflow,'.bonsai')));

animals = unique(extractBetween({task_dir.folder}, ...
    [plab.locations.server_data_path,filesep], ...
    filesep));


%% Testing: within-day dRxn

animals = ...
    {'AM001','AM002','AM004','AM005','AM006','AM007','AM008','AM009','AM011','AM012', ...
    'AM013','AM014','AM015','AM016','AM017','AM018','AM019','AM021','AM022','AM023', ...
    'AM024','AM025','AM026','AM027','AM029','AM030','AP004','AP005','AP006','AP007', ...
    'AP008','AP009','AP010','AP011','AP012','AP013','AP014','AP015','AP016','AP017', ...
    'AP018','AP020','AP021','AP022','AP023','AP025','DS001','DS007','DS009','DS010', ...
    'DS011'};

% Set reaction statistic to use
use_stat = 'mean';

% Grab reaction stats and learning day for each mouse
bhv = struct;

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};
    fprintf('%s (%d/%d)\n',animal,curr_animal_idx,length(animals))

    task_workflow = 'stim_wheel_right_stage\d';
    task_recordings = plab.find_recordings(animal,[],task_workflow);

    % Exclude recordings after any alternative workflows
    passive_workflows = {'lcr_passive','hml_passive','sparse_noise'};
    passive_recordings = plab.find_recordings(animal,[],passive_workflows);

    all_recordings = plab.find_recordings(animal);

    alt_task_days = setdiff({all_recordings.day},unique({task_recordings.day,passive_recordings.day}));

    if ~isempty(alt_task_days)
        use_task_recordings_idx = datetime({task_recordings.day}) < min(datetime(alt_task_days));
    else
        use_task_recordings_idx = true(size(task_recordings));
    end

    recordings = task_recordings(use_task_recordings_idx);

    stim_to_move_all = cell(length(recordings),1);
    stim_to_lastmove_all = cell(length(recordings),1);
    stim_to_outcome_all = cell(length(recordings),1);
    outcome_all = cell(length(recordings),1);

    rxn_firstmove_stat_p = nan(length(recordings),1);
    rxn_lastmove_stat_p = nan(length(recordings),1);
    rxn_stat = nan(length(recordings),1);
    rxn_null_stat = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        stim_to_move_all{curr_recording} = stim_to_move(1:n_trials);
        stim_to_move_all{curr_recording} = stim_to_lastmove(1:n_trials);
        stim_to_outcome_all{curr_recording} = stimOff_times(1:n_trials)-stimOn_times(1:n_trials);
        outcome_all{curr_recording} = trial_outcome(1:n_trials);

        % Get association stat
        % (skip if only a few trials)
        if n_trials < 10
            continue
        end
       
        % (rxn = first move)
        [rxn_firstmove_stat_p(curr_recording), ...
            rxn_stat(curr_recording),rxn_null_stat(curr_recording)] = ...
            AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

        % (rxn = last move)
        [rxn_lastmove_stat_p(curr_recording), ...
            rxn_stat(curr_recording),rxn_null_stat(curr_recording)] = ...
            AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_lastmove,use_stat);

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});
        AP_print_progress_fraction(curr_recording,length(recordings));

    end

    % Store behavior across animals
    bhv(curr_animal_idx).animal = animal;
    bhv(curr_animal_idx).stim_to_move = stim_to_move_all;
    bhv(curr_animal_idx).stim_to_lastmove = stim_to_lastmove_all;
    bhv(curr_animal_idx).stim_to_outcome = stim_to_outcome_all;
    bhv(curr_animal_idx).outcome = outcome_all;
    
    bhv(curr_animal_idx).rxn_firstmove_stat_p = rxn_firstmove_stat;
    bhv(curr_animal_idx).rxn_lastmove_stat_p = rxn_lastmove_stat;
    bhv(curr_animal_idx).rxn_stat = rxn_stat;
    bhv(curr_animal_idx).rxn_null_stat = rxn_null_stat;

end

ld_idx_cat = cell2mat(cellfun(@(x) (1:length(x))'-min([find(x,1),nan]),{bhv.learned_day}','uni',false));
ld_cat = vertcat(bhv.learned_day);
rxn_cat = vertcat(bhv.stim_to_move);

use_days = ld_idx_cat >= 0;

m1 = cellfun(@(x) x(2:end-1),rxn_cat(use_days),'uni',false); m1 = vertcat(m1{:});
m2 = cellfun(@(x) x(3:end)-x(1:end-2),rxn_cat(use_days),'uni',false); m2 = vertcat(m2{:});

figure; h = tiledlayout(1,2);


m1_bins = [0:0.05:2];
m1_bins_centers = m1_bins(1:end-1)+diff(m1_bins)/2;

m1_binned = discretize(m1,m1_bins);
r = ap.groupfun(@mean,m2,m1_binned);
s = ap.groupfun(@AP_sem,m2,m1_binned);

nexttile;ap.errorfill(m1_bins_centers,r,s);
yline(0);
ylabel('\Delta reaction (trial n-1 vs n+1)')
xlabel('Reaction (trial n)')


m2_bins = [-inf,-1:0.05:1,inf];
m2_bin_centers = m2_bins(1:end-1)+diff(m2_bins)/2;

m2_binned = discretize(m2,m2_bins);
r = ap.groupfun(@mean,m1,m2_binned);
s = ap.groupfun(@AP_sem,m1,m2_binned);

nexttile;ap.errorfill(m2_bin_centers(2:end-1),r(2:end-1),s(2:end-1));
xlabel('\Delta reaction (trial n-1 vs n+1)')
ylabel('Reaction (trial n)')
xline(0);



%% Run longstriatum package

clear
Marica_2025.package.save_wf_passive

clear
Marica_2025.package.save_wf_task






