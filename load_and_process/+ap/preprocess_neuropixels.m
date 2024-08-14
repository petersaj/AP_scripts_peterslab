function preprocess_neuropixels(animal,day)
% preprocess_neuropixels(animal,day)
% Currently assumes Neuropixels 3A or 1.0 recorded with Open Ephys

%% Set user locations

% Path to copy and sort locally
local_kilosort_path = fullfile(plab.locations.local_data_path,'kilosort');

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

%% Run kilosort on all datasets

for curr_data = 1:length(data_paths)  
    
    % Get experiments (= preview off in between)
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
        if ~exist(kilosort_output_path,'dir')
            mkdir(kilosort_output_path)
        end
        
        pyrunfile('AP_run_kilosort4.py', ...
            data_filename = apband_car_local_filename, ...
            kilosort_output_path = kilosort_output_path);

        % Terminate Python process
        terminate(kilosort_python_environment);

        %% Convert spike times to Open Ephys timestamps

        % Get metadata filename
        ephys_meta_dir = dir(fullfile(experiment_dir.folder,experiment_dir.name,'**','*.oebin'));
        ephys_meta_fn = fullfile(ephys_meta_dir.folder,ephys_meta_dir.name);
        ephys_metadata = jsondecode(fileread(ephys_meta_fn));

        % Convert kilosort "spike times" (samples) into timestamps
        ks_spike_times_fn = fullfile(kilosort_output_path,'spike_times.npy');
        oe_ap_samples_fn = fullfile(fileparts(ap_data_filename),'sample_numbers.npy');

        plab.ephys.ks2oe_timestamps(ks_spike_times_fn,oe_ap_samples_fn, ...
            ephys_metadata(1).continuous(1).sample_rate);
        
        %% Run bombcell (using CAR data)
   
        % Run bombcell
        kilosort_version = 4;
        ap.run_bombcell(apband_car_local_filename,kilosort_output_path,ephys_meta_fn,kilosort_version);
        
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


