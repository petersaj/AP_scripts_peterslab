% Re-run widefield preprocessing with dropped frames
%
% (correct for possible non-alternating colors, which wasn't compensated
% for in previous script, causing blue/violet to be mixed)


%% Find recordings with dropped frames
data_dir = dir(plab.locations.server_data_path);
animals_all = {data_dir(3:end).name};

reprocess_data_all = cell(0,2);
for curr_animal = 1:length(animals_all)

    animal = animals_all{curr_animal};
    animal_recordings = plab.find_recordings(animal);
    animal_days_wf = {animal_recordings(cellfun(@any,{animal_recordings.widefield})).day};

    for curr_day = 1:length(animal_days_wf)
        rec_day = animal_days_wf{curr_day};

        rec_times = plab.find_recordings(animal,rec_day).recording;

        for curr_rec = 1:length(rec_times)

            rec_time = rec_times{curr_rec};

            verbose = false;
            try
                ap.load_timelite
            catch me
                continue
            end

            widefield_metadata_fn = ...
                plab.locations.filename('server',animal,rec_day,[], ...
                'widefield',sprintf('widefield_%s_metadata.bin',rec_time));

            if ~exist(widefield_metadata_fn,'file')
                continue
            end
            fid = fopen(widefield_metadata_fn);
            widefield_metadata = reshape(fread(fid,'double'),9,[]);
            fclose(fid);

            wf_frames_rec = size(widefield_metadata,2);
            wf_frames_timelite = length(widefield_expose_times);

            if wf_frames_rec ~= wf_frames_timelite
                fprintf('Flagged %s %s (rec %s)\n',animal,rec_day,rec_time);
                reprocess_data_all(end+1,:) = {animal,rec_day};
                continue
            end

        end

        ap.print_progress_fraction(curr_day,length(animal_days_wf));

    end
    disp(['Done ' animal]);
end

[~,unique_idx] = unique(cell2mat(reprocess_data_all),'rows');
reprocess_data = reprocess_data_all(unique_idx,:);

%% Re-run preprocessing for the flagged files

clearvars -except reprocess_data

for curr_data = 1:size(reprocess_data,1)

    animal = reprocess_data{curr_data,1};
    rec_day = reprocess_data{curr_data,2};

    curr_process_path = plab.locations.filename('server',animal,rec_day,[],'widefield');

    preload_vars = who;

    fprintf('Preprocessing: %s \n',curr_process_path);

    %% SVD decomposition of widefield data

    % Currently: SVD 2 colors separately, assume alt. blue/violet
    wf_colors = {'blue','violet'};
    n_colors = length(wf_colors);

    [U,V,im_avg] = plab.wf.preprocess_widefield(curr_process_path,n_colors);

    %% Save preprocessed widefield data locally

    local_save_path = fullfile(plab.locations.local_data_path,animal,rec_day,'widefield');

    % Clear local path
    mkdir(local_save_path);

    % Set color names (empty if SVD not split by color)
    if n_colors == 1
        color_suffix = {''};
    else
        color_suffix = cellfun(@(x) sprintf('_%s',x),wf_colors,'uni',false);
    end

    % Save mean images in experiment folder by color
    for curr_color = 1:length(color_suffix)
        curr_mean_im_fn = fullfile(local_save_path, ...
            sprintf('meanImage%s.npy',color_suffix{curr_color}));
        writeNPY(im_avg{curr_color},curr_mean_im_fn);
    end

    % Save spatial components in experiment (animal/day) folder by color
    for curr_color = 1:length(color_suffix)
        curr_U_fn = fullfile(local_save_path, ...
            sprintf('svdSpatialComponents%s.npy',color_suffix{curr_color}));
        writeNPY(U{curr_color},curr_U_fn);
    end

    % Save temporal components in associated recording folders
    % (get recording times from filenames)
    data_dir = dir(fullfile(curr_process_path,'*_data.bin'));
    recording_times = extract({data_dir.name},digitsPattern);
    for curr_recording = 1:size(V,1)
        for curr_color = 1:length(color_suffix)
            % Put V's into separate recording paths
            curr_V_fn = fullfile(local_save_path,recording_times{curr_recording}, ...
                sprintf('svdTemporalComponents%s.npy',color_suffix{curr_color}));
            % Make recording path
            if ~exist(fileparts(curr_V_fn),'dir')
                mkdir(fileparts(curr_V_fn));
            end
            % Write V to file
            writeNPY(V{curr_recording,curr_color},curr_V_fn);
        end
    end

    %% Move data onto server
    % Check if the server is available
    if ~exist(plab.locations.server_data_path,'dir')
        warning('Server not accessible at %s',plab.locations.server_data_path)
        return
    end

    % Move local data to server:
    local_data_dir = dir(local_save_path);

    % Move day-level files to day widefield folder
    disp('Moving widefield files to server...')
    local_data_dayfiles_idx = ~[local_data_dir.isdir];
    for curr_file_idx = find(local_data_dayfiles_idx)

        % Local filename
        curr_local_filename = fullfile(local_data_dir(curr_file_idx).folder, ...
            local_data_dir(curr_file_idx).name);

        % Server filename: replace path from local data to server
        curr_server_filename = strrep(curr_local_filename, ...
            plab.locations.local_data_path,plab.locations.server_data_path);

        % Make server path (if it doesn't exist) and move
        if ~exist(fileparts(curr_server_filename),'dir')
            mkdir(fileparts(curr_server_filename));
        end

        [status,message] = movefile(curr_local_filename,curr_server_filename);
        if ~status
            warning('Failed moving to server: %s',message);
        else
            fprintf('Moved %s --> %s \n',curr_local_filename,curr_server_filename);
        end
    end

    % Move protocol-level files to protocol folders
    local_data_protocolfiles_idx = [local_data_dir.isdir] & ~contains({local_data_dir.name},'.');

    for curr_recording_idx = find(local_data_protocolfiles_idx)

        % Local path (name is recording time)
        curr_local_path = fullfile(local_data_dir(curr_recording_idx).folder, ...
            local_data_dir(curr_recording_idx).name);
        rec_time = local_data_dir(curr_recording_idx).name;

        % Server path
        curr_server_path = plab.locations.filename('server', ...
            animal,rec_day,rec_time,'widefield');

        % Make server path (if it doesn't exist) and move
        if ~exist(fileparts(curr_server_path),'dir')
            mkdir(fileparts(curr_server_path));
        end

        [status,message] = copyfile(curr_local_path,curr_server_path);
        if ~status
            warning('Failed moving to server: %s',message);
        else
            fprintf('Moved %s --> %s \n',curr_local_path,curr_server_path);
        end

    end

    %% Clear local data
    rmdir(fileparts(fileparts(local_save_path)),'s')

end













