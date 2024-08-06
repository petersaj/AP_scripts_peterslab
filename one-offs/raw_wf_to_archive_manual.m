% NOTE: 
% This was used when manually transferring backlog. 
% This script was converted into PowerShell and set to run every day.

%% Set destination folders
archive_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Archive\to_archive';

%% Widefield: select files

% Set widefield file formats
% ('widefield_(rec time)_data.bin is from plab.rig.widefield)
% (tif/dcimg is from HCImage)
raw_wf_format = {'*.tif','*.dcimg','widefield_*_data.bin'};

% Find all raw widefield files
raw_wf_files = cellfun(@(ext) dir(fullfile(plab.locations.server_data_path, ...
    '**','widefield',sprintf('%s',ext))),raw_wf_format,'uni',false);
raw_wf_files_all = vertcat(raw_wf_files{:});

fprintf('Total raw size: %.2f TB \n',sum([raw_wf_files_all.bytes])/1e12);

%%%%%% SELECT SUBSET OF FILES BY RELATIVE DATE (2 weeks before today)
end_date = datetime('today')-14;
move_file_idx = datetime({raw_wf_files_all.date}) < end_date;
raw_size_total = sum([raw_wf_files_all(move_file_idx).bytes])/1e12;
raw_size_cumulative = cumsum([raw_wf_files_all(move_file_idx).bytes]/1e12);

fprintf('Selected raw size: %.2f TB \n',raw_size_total);
%%%%%%%%%%%%

% %%%%%% SELECT SUBSET OF FILES BY DATE
% end_date = datetime('15-Jul-2024');
% move_file_idx = datetime({raw_wf_files_all.date}) < end_date;
% raw_size_total = sum([raw_wf_files_all(move_file_idx).bytes])/1e12;
% raw_size_cumulative = cumsum([raw_wf_files_all(move_file_idx).bytes]/1e12);
% 
% fprintf('Selected raw size: %.2f TB \n',raw_size_total);
% %%%%%%%%%%%%

% %%%%%% SELECT SUBSET OF FILES BY DATE-SORTED SIZE
% max_size = 17.6; % TB
% [~,date_sort_idx] = sort(datetime({raw_wf_files_all.date}));
% date_size_cumulative = cumsum([raw_wf_files_all(date_sort_idx).bytes]/1e12);
% last_max_file = find(date_size_cumulative < max_size,1,'last');
% 
% end_date = datetime(raw_wf_files_all(date_sort_idx(last_max_file)).date);
% move_file_idx = datetime({raw_wf_files_all.date}) < end_date;
% raw_size_total = sum([raw_wf_files_all(move_file_idx).bytes])/1e12;
% raw_size_cumulative = cumsum([raw_wf_files_all(move_file_idx).bytes]/1e12);
% 
% fprintf('Selected raw size: %.2f TB \n',raw_size_total);
% %%%%%%%%%%%%%

%% Widefield: move files to archive folder

move_files = raw_wf_files_all(move_file_idx);

for curr_move_file_idx = 1:length(move_files)

    % Get data folder
    source_folder = move_files(curr_move_file_idx).folder;

    % Pull out subfolders (e.g. animal,day,'widefield')
    source_subfolder = split(erase(source_folder, ...
        [plab.locations.server_data_path,filesep]),filesep);

    % Construct destination folder from subfolders
    destination_folder = fullfile(archive_path,strjoin(source_subfolder,'_'));
    if ~exist(destination_folder,'dir')
        mkdir(destination_folder)
    end

    curr_source = fullfile(source_folder,move_files(curr_move_file_idx).name);
    curr_destination = fullfile(destination_folder,move_files(curr_move_file_idx).name);

    % Move file (if it exists), print status
    if exist(curr_source,'file')
        fprintf('Moving: %s --> \n%s\n',curr_source,curr_destination);
        tic
        [move_status,move_message] = movefile(curr_source,curr_destination);
        move_time = toc;
        if move_status
            curr_file_size = move_files(curr_move_file_idx).bytes/1e9;
            curr_copy_rate = curr_file_size/(move_time/60);
            fprintf('Last: %.2fGB, %.2fGB/min\n',curr_file_size,curr_copy_rate);
            fprintf('Total: %.2fTB/%.2fTB, %.1f hours remaining\n',raw_size_cumulative(curr_move_file_idx),raw_size_total, ...
                (raw_size_total - raw_size_cumulative(curr_move_file_idx))*1e3/curr_copy_rate/60);
        else
            disp(move_message)
            continue
        end
    end
end

disp('Finished transfer to archive folder');


%% TEMPORARY: copy backup to external hard drive
% 
% backup_server_path = '\\qnap-ap002.dpag.ox.ac.uk\temp';
% external_hd_path = 'E:\raw_widefield_backup';
% 
% widefield_backup_folders = dir(fullfile(backup_server_path,'*widefield*'));
% 
% for curr_folder = 1:length(widefield_backup_folders)
%     curr_source = fullfile(widefield_backup_folders(curr_folder).folder, ...
%         widefield_backup_folders(curr_folder).name);
%     curr_destination = fullfile(external_hd_path, ...
%         widefield_backup_folders(curr_folder).name);
% 
%     curr_size = sum([dir(curr_source).bytes])/1e9;
% 
%     disp('Copying folder to external hd....')
%     tic
%     [status,message] = copyfile(curr_source,curr_destination);
%     copy_time = toc;
%     if status
%         curr_copy_time = curr_size/(copy_time/60);
%         fprintf('Copied (%.2fGB,%.2fGB/min): %s --> %s\n',curr_size,curr_copy_time,curr_source,curr_destination);
%     else
%         disp(message)
%         continue
%     end
% 
% end




















