% ap.load_recording
%
% Load and prepare all/selected parts of a recording

%% Check that the necessary parameters are included

if ~exist('animal','var') || isempty(animal)
    error('Animal not defined')
end
if ~exist('rec_day','var') || isempty(animal)
    error('Recording day (rec_day) not defined')
end
if ~exist('rec_time','var') || isempty(animal)
    error('Recording time (rec_time) not defined');
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













