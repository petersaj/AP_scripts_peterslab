function recordings = find_recordings(animal,recording_date,workflow)
% recordings = find_recordings(animal,recording_date,workflow)
%
% Find recordings for a given animal
%
% Input: recording date (as 'yyyy-mm-yy'), Bonsai workflow (can be multiple {'workflow 1','workflow 2'})
%   - if date filled and workflow empty: find all recordings within one day
%   - if workflow filled and date empty: find Bonsai workflow over all days
%
% Output:
% recordings (struct):
% .day - date of recording
% .protocol - protocol folder name(s) (time as HHMM)
% .index - index of protocol within day (e.g. 3rd protocol is [3])
% .workflow - Bonsai workflow
% .mousecam - whether mousecam was recorded for each protocol
% .widefield - whether widefield was recorded for each protocol

%% Initialize parameters

% Check date pattern
if ~isempty(recording_date)
    date_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
    if ~matches(recording_date,date_pattern)
        error('Recording date (''%s'') not ''yyyy-mm-dd'' format',recording_date);
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

% Set flag for search type (date and/or workflow)
if ~isempty(recording_date)
    search_type = 'date';
elseif isempty(recording_date) && ~isempty(workflow)
    search_type = 'workflow';
end


%% Perform search depending on search type

switch search_type

    case 'date'
        %% Find workflows in specific date

        curr_path = plab.locations.make_server_filename(animal,recording_date);

        % Get protocol folders within day
        curr_protocol_paths = dir(fullfile(curr_path,'Protocol*'));

        % Get bonsai workflows for each protocol
        protocol_workflows = cell(size(curr_protocol_paths));
        for curr_path_idx = 1:length(curr_protocol_paths)
            curr_bonsai_path = fullfile( ...
                curr_protocol_paths(curr_path_idx).folder, ...
                curr_protocol_paths(curr_path_idx).name,'bonsai');

            curr_bonsai_dir = dir(curr_bonsai_path);
            curr_workflow_idx = [curr_bonsai_dir.isdir] & ...
                ~contains({curr_bonsai_dir.name},'.');

            if any(curr_workflow_idx)
                curr_workflow = curr_bonsai_dir(curr_workflow_idx).name;
                protocol_workflows{curr_path_idx} = curr_workflow;
            end
        end

        % Return all workflows (if unspecified), or specified workflows
        if ~isempty(workflow)
            use_protocols = ismember(protocol_workflows,workflow);
        else
            use_protocols = true(size(curr_protocol_paths));
        end

        % Package info
        recordings.day = recording_date;
        recordings.index = find(use_protocols);
        recordings.protocol = strtok({curr_protocol_paths(use_protocols).name},'Protocol_');
        recordings.workflow = protocol_workflows(use_protocols);

        recordings.mousecam = ...
            cellfun(@(x) any(exist(fullfile(curr_path,x,'mousecam'),'dir')), ...
            {curr_protocol_paths(use_protocols).name});
        recordings.widefield = ...
            cellfun(@(x) any(exist(fullfile(curr_path,x,'widefield'),'dir')), ...
            {curr_protocol_paths(use_protocols).name});


    case 'workflow'
        %% Find workflow across all days

        % Get contents of animal path
        animal_path = fullfile(plab.locations.server_data_path,animal);
        animal_dir = dir(animal_path);

        % Find recording paths (matches date format and is folder)
        date_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
        recording_day_idx = matches({animal_dir.name},date_pattern) & [animal_dir.isdir];
        recording_days = {animal_dir(recording_day_idx).name};

        % Loop through recordings for specified workflow and grab modalities
        struct_fieldnames = {'day','protocol','mousecam','widefield'};
        recordings = cell2struct(cell(length(struct_fieldnames),0), ...
            struct_fieldnames);

        for curr_day = 1:length(recording_days)

            curr_path = plab.locations.make_server_filename(animal,recording_days{curr_day});

            % Get protocol folders within day
            curr_protocol_paths = dir(fullfile(curr_path,'Protocol*'));

            % Look for matches with specified Bonsai workflow
            % (and store the workflow name for each protocol)
            protocol_workflows = cell(size(curr_protocol_paths));
            use_protocol_paths = false(size(curr_protocol_paths));
            for curr_path_idx = 1:length(curr_protocol_paths)
                curr_bonsai_path = fullfile( ...
                    curr_protocol_paths(curr_path_idx).folder, ...
                    curr_protocol_paths(curr_path_idx).name,'bonsai');

                curr_bonsai_dir = dir(curr_bonsai_path);
                curr_workflow_idx = [curr_bonsai_dir.isdir] & ...
                    ~contains({curr_bonsai_dir.name},'.');

                if any(curr_workflow_idx)
                    curr_workflow = curr_bonsai_dir(curr_workflow_idx).name;
                    protocol_workflows{curr_path_idx} = curr_workflow;

                    if ismember(curr_workflow,workflow)
                        use_protocol_paths(curr_path_idx) = true;
                    end
                end
            end

            % If there were Bonsai matches, package info
            if any(use_protocol_paths)

                curr_recording_idx = length(recordings) + 1;

                recordings(curr_recording_idx).day = recording_days{curr_day};
                recordings(curr_recording_idx).index = find(use_protocol_paths);
                recordings(curr_recording_idx).protocol = ...
                    strtok({curr_protocol_paths(use_protocol_paths).name},'Protocol_');
                recordings(curr_recording_idx).workflow = ...
                    protocol_workflows(use_protocol_paths);

                % Check which modalities were recorded
                recordings(curr_recording_idx).mousecam = ...
                    cellfun(@(x) any(exist(fullfile(curr_path,x,'mousecam'),'dir')), ...
                    {curr_protocol_paths(use_protocol_paths).name});

                recordings(curr_recording_idx).widefield = ...
                    cellfun(@(x) any(exist(fullfile(curr_path,x,'widefield'),'dir')), ...
                    {curr_protocol_paths(use_protocol_paths).name});
            end
        end
end

















