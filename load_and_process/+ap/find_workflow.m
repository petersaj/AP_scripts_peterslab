function recordings = find_workflow(animal,workflow_name)
% recordings = find_workflow(animal,workflow_name)
%
% Find instances of bonsai workflows for a given animal
% 
% Input: workflow_name - can be multiple {'workflow 1','workflow 2'};
%
% Output: 
% recordings (struct): 
% .day - date of recording
% .protocol - protocol folder name(s) for each workflow
% .mousecam - whether mousecam was recorded for each protocol
% .widefield - whether widefield was recorded for each protocol

%% Standardize inputs

% If workflow name isn't a cell, make it a cell
if ~iscell(workflow_name)
    workflow_name = {workflow_name};
end

%% Find recording folders

% Get contents of animal path
animal_path = fullfile(plab.locations.server_data_path,animal);
animal_dir = dir(animal_path);

% Find recording paths (matches date format and is folder)
recording_day_idx = cellfun(@(x) ...
    ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{animal_dir.name}) &...
    [animal_dir.isdir];
recording_days = {animal_dir(recording_day_idx).name};


%% Find specified workflow
% Loop through recordings for specified workflow and grab modalities

struct_fieldnames = {'day','protocol','mousecam','widefield'};
recordings = cell2struct(cell(length(struct_fieldnames),0), ...
    struct_fieldnames);

for curr_day = 1:length(recording_days)  
    
    curr_path = plab.locations.make_server_filename(animal,recording_days{curr_day});

    % Get protocol folders within day
    curr_protocol_paths = dir(fullfile(curr_path,'Protocol*'));

    % Look for matches with specified Bonsai workflow
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
            if ismember(curr_workflow,workflow_name)
                use_protocol_paths(curr_path_idx) = true;
            end
        end
    end

    % If there were Bonsai matches, store info
    if any(use_protocol_paths)

        curr_recording_idx = length(recordings) + 1;

        recordings(curr_recording_idx).day = recording_days{curr_day};
        recordings(curr_recording_idx).protocol = ...
            strtok({curr_protocol_paths(use_protocol_paths).name},'Protocol_');

        % Check which modalities were recorded
        recordings(curr_recording_idx).mousecam = ...
            cellfun(@(x) any(exist(fullfile(curr_path,x,'mousecam'),'dir')), ...
            {curr_protocol_paths(use_protocol_paths).name});

        recordings(curr_recording_idx).widefield = ...
            cellfun(@(x) any(exist(fullfile(curr_path,x,'widefield'),'dir')), ...
            {curr_protocol_paths(use_protocol_paths).name});
    end

end
    


















