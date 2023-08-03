function recordings = find_recordings(animal,recording_day,workflow)
% recordings = find_recordings(animal,recording_day,workflow)
%
% Find recordings for a given animal/day/workflow, package info by day.
%
% Input: 
% recording_day - (as 'yyyy-mm-yy'), can be multiple in a cell array
% workflow - can be multiple (e.g. {'thisone','thatone'}), can include * as
% wildcard (e.g. 'this*' will return 'thisone'). 
%
%  - if day filled and workflow empty: find all recordings within one day
%  - if workflow filled and day empty: find Bonsai workflow over all days
%
% Output:
% recordings (struct: length = n days):
% .day - day of recording
% .recording - recording folder name(s) (time as HHMM)
% .index - index(ies) of recording within day (e.g. 3rd recording is [3])
% .workflow - Bonsai workflow(s)
% .mousecam - whether mousecam was recorded for recording(s)
% .widefield - whether widefield was recorded for recording(s)
%
% e.g. recording(4).workflow{3}: the 3rd matching workflow from the 4th day
% containing a matching recording.

%% Initialize parameters

% Check day pattern
if ~isempty(recording_day)
    day_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
    if ~matches(recording_day,day_pattern)
        error('Recording day (''%s'') not ''yyyy-mm-dd'' format',recording_day);
    end
end

% Set workflow empty if not defined
if ~exist('workflow','var')
    workflow = [];
elseif ~isempty(workflow)
    % If workflow defined and not a cell, make it a cell
    if ~iscell(workflow)
        workflow = {workflow};
    end
end

% Set flag for search type (day and/or workflow)
if ~isempty(recording_day)
    search_type = 'day';
elseif isempty(recording_day) && ~isempty(workflow)
    search_type = 'workflow';
end

% Set days to search (all days if empty, specified day(s) if not)
if isempty(recording_day)
    % Get contents of animal path
    animal_path = fullfile(plab.locations.server_data_path,animal);
    animal_dir = dir(animal_path);

    % Find recording paths (matches day format and is folder)
    day_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
    recording_day_idx = matches({animal_dir.name},day_pattern) & [animal_dir.isdir];
    recording_day = {animal_dir(recording_day_idx).name};
elseif ~iscell(recording_day)
    recording_day = {recording_day};
end

%% Find and package matching recordings

struct_fieldnames = {'day','index','recording','workflow','mousecam','widefield'};
recordings = cell2struct(cell(length(struct_fieldnames),0),struct_fieldnames);

for curr_day_idx = 1:length(recording_day)
    curr_day = recording_day{curr_day_idx};

    curr_path = plab.locations.make_server_filename(animal,curr_day);

    % Get recording folders within day
    curr_recording_paths = dir(fullfile(curr_path,'Recording*'));

    % Get bonsai workflows for each recording
    recording_workflows = cell(size(curr_recording_paths));
    for curr_path_idx = 1:length(curr_recording_paths)
        curr_bonsai_path = fullfile( ...
            curr_recording_paths(curr_path_idx).folder, ...
            curr_recording_paths(curr_path_idx).name,'bonsai');

        curr_bonsai_dir = dir(curr_bonsai_path);
        curr_workflow_idx = [curr_bonsai_dir.isdir] & ...
            ~contains({curr_bonsai_dir.name},'.');

        if any(curr_workflow_idx)
            curr_workflow = curr_bonsai_dir(curr_workflow_idx).name;
            recording_workflows{curr_path_idx} = curr_workflow;
        end
    end

    % Return all workflows (if unspecified), or specified workflows
    if ~isempty(workflow)
        % Regexp: workflow1$|workflow2$... 
        % ($ = string terminus, | looks for OR)
        % (converts * into .* - allows for terminus)
        workflow_regexp = strjoin(append(strrep(workflow,'*','.*'),'$'),'|');
        use_recordings = ...
            ~cellfun(@isempty,regexp(recording_workflows, ...
            workflow_regexp,'forceCellOutput'));
    else
        use_recordings = true(size(curr_recording_paths));
    end

    % If no recordings, move to next day
    if ~any(use_recordings)
        continue
    end

    % Package info
    recording_idx = length(recordings) + 1;

    recordings(recording_idx).day = curr_day;
    recordings(recording_idx).index = find(use_recordings);
    recordings(recording_idx).recording = ...
        strtok({curr_recording_paths(use_recordings).name},'Recording_');
    recordings(recording_idx).workflow = recording_workflows(use_recordings);

    recordings(recording_idx).mousecam = ...
        cellfun(@(x) any(exist(fullfile(curr_path,x,'mousecam'),'dir')), ...
        {curr_recording_paths(use_recordings).name});
    recordings(recording_idx).widefield = ...
        cellfun(@(x) any(exist(fullfile(curr_path,x,'widefield'),'dir')), ...
        {curr_recording_paths(use_recordings).name});

end







































