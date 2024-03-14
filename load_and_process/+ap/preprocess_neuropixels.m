function preprocess_neuropixels(animal,day)
% preprocess_neuropixels(animal,day)
% Currently assumes Neuropixels 3A or 1.0 recorded with Open Ephys

%% Set user locations

% Path to copy and sort locally
local_kilosort_path = 'D:\data_temp\kilosort';

% Path to kilosort python environment
system_path = split(getenv("Path"),';');
anaconda_path = system_path{endsWith(system_path,'anaconda3')};
kilosort_environment_path = fullfile(anaconda_path,'envs','kilosort','pythonw.exe');
kilosort_python_environment = pyenv('Version',kilosort_environment_path,'ExecutionMode','OutOfProcess');

% If no day selected, default to today
if nargin < 2 || isempty(day)
    day = char(datetime('today','Format','yyyy-MM-dd'));
end

%% Get paths and filenames

ephys_path = plab.locations.filename('server',animal,day,[],'ephys');

ephys_exists = exist(ephys_path,'dir');

if ~ephys_exists
    error([animal ' ' day ': No ephys data found']);
end

save_paths = {[ephys_path filesep 'kilosort4']};
data_paths = {ephys_path};

% Get ephys recording paths
% probe_n = multiple probes simultaneously
% site_n = multiple sites recorded in serial
data_path_dir = ...
    [dir(fullfile(data_paths{1}, 'probe_*')), ...
    dir(fullfile(data_paths{1}, 'site_*'))];
if ~isempty(data_path_dir)
    data_paths = cellfun(@(x) [data_paths{1} filesep x],{data_path_dir.name},'uni',false);
    save_paths = cellfun(@(x) [save_paths{1} filesep x],{data_path_dir.name},'uni',false);
end

for curr_data = 1:length(data_paths)  
    
    % Get experiments (if turned off between)
    curr_data_path = data_paths{curr_data};   
    ephys_exp_paths = dir([curr_data_path filesep 'experiment*']);
       
    for curr_exp = 1:length(ephys_exp_paths)

        curr_exp_path = fullfile(ephys_exp_paths(curr_exp).folder, ...
            ephys_exp_paths(curr_exp).name);
        
        % Update save path with experiment (only if more than one, legacy)
        if length(ephys_exp_paths) == 1
            curr_save_path = save_paths{curr_data};
        elseif length(ephys_exp_paths) > 1
            curr_save_path = [save_paths{curr_data} filesep ephys_exp_paths(curr_exp).name];
        end
        
        % Make save path
        if ~exist(curr_save_path,'dir')
            mkdir(curr_save_path)
        end
              
        % Find Open Ephys recordings
        experiment_dir = dir(fullfile(curr_exp_path,'recording*'));

        % If more than one recording, error out at the moment
        % (multiple recordings = record stopped/started, preview continuous)
        % (in future, could manually concatenate them - kilosort allows
        % for py.lists of filenames, but it assumes chronic so performs
        % some chronic drift correction which is broken)
        if length(experiment_dir) > 1
            error('%s %s: multiple recordings, not handled in kilosort code yet',animal,day);
        end

        % Get Open Ephys filenames
        ap_data_dir = dir(fullfile(experiment_dir.folder,experiment_dir.name, ...
            'continuous', '*-AP', 'continuous.dat'));
        ap_data_filename = fullfile(ap_data_dir.folder,ap_data_dir.name);


        %% Run kilosort
        
        % Clear out local kilosort directories
        if exist(local_kilosort_path,'dir')
            rmdir(local_kilosort_path,'s')
        end
        mkdir(local_kilosort_path);
        
        % Copy AP-band data locally
        disp('Copying AP data to local drive...') 
        apband_local_filename = fullfile(local_kilosort_path, ...
            sprintf('%s_%s_apband.dat',animal,day));

        copyfile(ap_data_filename, apband_local_filename);
        disp('Done');
            
        % Run common average referencing (CAR)
        apband_car_local_filename = strrep(apband_local_filename,'.dat','_car.dat');
        ap.ephys_car(apband_local_filename,apband_car_local_filename)

        % Run Kilsort 4
        disp('Running Kilosort 4...');
        kilosort_output_path = fullfile(local_kilosort_path,'kilosort_output');
        pyrunfile('AP_run_kilosort4.py', ...
            data_filename = apband_car_local_filename, ...
            kilosort_output_path = kilosort_output_path);

        %% Convert spike times to Open Ephys timestamps
        % This has two advantages: 
        % - sync timestamps can be used as-is without applying offset
        % - compensates for dropped samples

        % Load AP-band timestamps
        openephys_ap_timestamps_filename = fullfile(fileparts(ap_data_filename),'timestamps.npy');
        openephys_ap_timestamps = readNPY(openephys_ap_timestamps_filename);

        % Convert kilosort spike time ouput (as sample index) into open
        % ephys timestamp (in seconds)
        spike_times_kilosort_filename = fullfile(kilosort_output_path,'spike_times.npy');
        spike_times_kilosort = readNPY(spike_times_kilosort_filename);

        % NOTE: sometimes kilsort outputs indicies of spike times which are
        % not in the range of recording?! Give a warning and calculate
        % those times with the sample rate
        spike_times_kilosort_validtime = ...
            spike_times_kilosort >= 1 & ...
            spike_times_kilosort <= length(openephys_ap_timestamps);

        if all(spike_times_kilosort_validtime)
            spike_times_openephys = openephys_ap_timestamps(spike_times_kilosort);
        else
            ap_sample_time = median(diff(openephys_ap_timestamps));
            warning('Kilosort %s %s: %d spike times out of time range of data, interpolating', ....
                animal,day,sum(~spike_times_kilosort_validtime));
            spike_times_openephys = nan(size(spike_times_kilosort));
            spike_times_openephys(spike_times_kilosort_validtime) = ...
                openephys_ap_timestamps(spike_times_kilosort(spike_times_kilosort_validtime));

            % (interpolate from first and last sample times given rate)
            spike_times_openephys(~spike_times_kilosort_validtime) = ...
                interp1([0,1,length(openephys_ap_timestamps),length(openephys_ap_timestamps)+1], ...
                sort(reshape([openephys_ap_timestamps([1,end]),openephys_ap_timestamps([1,end])+([-1;1].*ap_sample_time)],1,[])), ...
                double(spike_times_kilosort(~spike_times_kilosort_validtime)),'linear','extrap');
        end

        % Save open ephys spike times into kilosort output folder
        spike_times_openephys_filename = fullfile(kilosort_output_path,'spike_times_openephys.npy');
        writeNPY(spike_times_openephys,spike_times_openephys_filename);

        %% Run bombcell (using CAR data)

        % Get metadata filename
        ephys_meta_dir = dir(fullfile(experiment_dir.folder,experiment_dir.name,'**','*.oebin'));
        ephys_meta_fn = fullfile(ephys_meta_dir.folder,ephys_meta_dir.name);

        % Run bombcell
        ap.run_bombcell(apband_car_local_filename,kilosort_output_path,ephys_meta_fn);
        
        %% Copy kilosort results to server
                
        disp('Copying sorted data to server...');
        copyfile(kilosort_output_path,curr_save_path);
        
        %% Delete all temporary local data
        
        rmdir(local_kilosort_path,'s');
        mkdir(local_kilosort_path);
        
    end
    
end

%% Print end message
fprintf('\nDone preprocessing Neuropixels: %s %s\n',animal,day);


