% ap.load_recording
%
% Load and prepare all/selected parts of a recording

%% Check that the necessary parameters are included

% If no animal - choose from list of animals
if ~exist('animal','var') || isempty(animal)
    if ~exist('rec_day','var') || isempty(rec_day)
        % (if no day - return all animals)
        data_dir = dir(plab.locations.server_data_path);
        animals_all = {data_dir(3:end).name};
    else
        % (if day specified - return animals with those days)
        day_dir = dir(fullfile(plab.locations.server_data_path,'*',rec_day));
        animals_all = erase(unique({day_dir.folder})',{plab.locations.server_data_path,rec_day,filesep});
    end
    animal_idx = listdlg('PromptString','Select animal:', ...
        'ListString',animals_all,'ListSize',[300,200], ...
        'SelectionMode','single');
    animal = animals_all{animal_idx};
    verbose = true;
end

% If no day - choose from list of days
if ~exist('rec_day','var') || isempty(rec_day)
    animal_recordings = plab.find_recordings(animal);
    recording_days = {animal_recordings.day};
    day_idx = listdlg('PromptString','Select day:', ...
        'ListString',recording_days,'ListSize',[300,200], ...
        'SelectionMode','single');
    rec_day = recording_days{day_idx};
    verbose = true;
end

% If no time - choose from list of workflows
if ~exist('rec_time','var') || isempty(rec_time)
    recordings = plab.find_recordings(animal,rec_day);
    rec_idx = listdlg('PromptString','Select workflow:', ...
        'ListString',recordings.workflow,'ListSize',[300,200], ...
        'SelectionMode','single');
    rec_time = recordings.recording{rec_idx};
    verbose = true;
end


%% Define what to load

if ~exist('verbose','var')
    verbose = false;
end

if verbose; fprintf('Loading %s, %s, Recording %s\n', animal, rec_day, rec_time); end;

% If nothing specified, load everything (but not LFP)
if ~exist('load_parts','var')
    load_parts.mousecam = true;
    load_parts.widefield = true;
    load_parts.ephys = true;
else
    % If only some things specified, don't load others
    if ~isfield(load_parts,'mousecam')
        load_parts.mousecam = false;
    end
    if ~isfield(load_parts,'widefield')
        load_parts.widefield = false;
    end
    if ~isfield(load_parts,'ephys')
        load_parts.ephys = false;
    end
end

% Get datetime of selected recording
rec_datetime = datetime(strjoin({rec_day,rec_time}), ...
                'InputFormat','yyyy-MM-dd HHmm');

%% Load experiment components

% Load timelite and associated inputs
ap.load_timelite

% Load Bonsai
ap.load_bonsai

% Load mousecam
if load_parts.mousecam
    ap.load_mousecam
end

% Load widefield
if load_parts.widefield && ...
        exist(plab.locations.filename('server',animal,rec_day,rec_time,'widefield'),'dir')
    ap.load_widefield
end

% Load ephys
if load_parts.ephys && ...
        exist(plab.locations.filename('server',animal,rec_day,[],'ephys'),'dir')
    ap.load_ephys
end

if verbose; disp('Finished.'); end;













