function preprocess_neuropixels(animal,day)
% preprocess_neuropixels(animal,day)
% Currently assumes Neuropixels Phase 3A recorded with Open Ephys

%% Set user locations

% Path to copy and sort locally
local_kilosort_path = 'D:\data_temp\kilosort';

% Path to pykilosort python environment
pykilosort_python_environment = pyenv('Version','C:\Users\petersa\anaconda3\envs\pyks2\pythonw.exe');


%% Get paths and filenames

ephys_path = plab.locations.filename('server',animal,day,[],'ephys');

ephys_exists = exist(ephys_path,'dir');

if ~ephys_exists
    error([animal ' ' day ': No ephys data found']);
end

save_paths = {[ephys_path filesep 'pykilosort']};
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
        % (in future, could manually concatenate them - pykilosort allows
        % for py.lists of filenames, but it assumes chronic so performs
        % some chronic drift correction which is broken)
        if length(experiment_dir) > 1
            error('%s %s: multiple recordings, not handled in kilosort code yet',animal,day);
        end

        % Get Open Ephys filenames
        ap_data_filename = fullfile(experiment_dir.folder,experiment_dir.name, ...
            'continuous', 'Neuropix-3a-100.Neuropix-3a-AP', 'continuous.dat');

%         % (sync not used - loaded when loading experiment)
%         sync_filename = fullfile(experiment_dir.folder,experiment_dir.name, ...
%             'events', 'Neuropix-3a-100.Neuropix-3a-AP', 'TTL', 'states.npy');
%         sync_timestamps_filename = fullfile(experiment_dir.folder,experiment_dir.name, ...
%             'events', 'Neuropix-3a-100.Neuropix-3a-AP', 'TTL', 'timestamps.npy');
        

        %% Get and save recording parameters
        
        % The gains and filter cuts aren't recorded anymore?!
        ap_gain = {500};
        lfp_gain = {125};
        filter_cut = {300};
        
        % (0.195x for int16 to uV? how's this change with gain, just another x?)
        
        % Hard-coded parameters
        n_channels = 384;
        ap_sample_rate = 30000;
        lfp_sample_rate = 2500;
        
        params = {'raw_path',['''' curr_data_path '''']; ...
            'n_channels',num2str(n_channels); ...
            'ap_sample_rate',num2str(ap_sample_rate); ... % this should be 30000 AP, 2500 LFP
            'lfp_sample_rate',num2str(lfp_sample_rate);
            'ap_gain',num2str(ap_gain{1}); ...
            'lfp_gain',num2str(lfp_gain{1})
            'filter_cutoff',num2str(filter_cut{1})};
        
        param_filename = [curr_save_path filesep 'dat_params.txt'];
        
        formatSpec = '%s = %s\r\n';
        fid = fopen(param_filename,'w');
        for curr_param = 1:size(params,1)
            fprintf(fid,formatSpec,params{curr_param,:});
        end
        fclose(fid);

        %% Get and save digital input events
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%         Obsolete: original values loaded in loading script
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         % NOTE: open ephys gives timestamps relative to start of PREVIEW,
%         % kilosort only gives timestamps relative to beginning of RECORD.
%         % This means that to match up kilosort (t(1) = 0) to OE (t(1) =
%         % arbitrary), the sync timestamps are saved as sync_timestamps -
%         % ap_timestamps(1)
%               
%         % Load spike timestamps
%         ap_timestamps_filename = fullfile(fileparts(ap_data_filename),'timestamps.npy');
%         ap_timestamps = readNPY(ap_timestamps_filename);
%                     
%         % Load digital input event times
%         sync_data = readNPY(sync_filename);
%         sync_timestamps = readNPY(sync_timestamps_filename);
%         
%         sync_channels = unique(abs(sync_data));
%         open_ephys_ttl = struct('timestamps',cell(size(sync_channels)),'values',cell(size(sync_channels)));
%         for curr_sync = 1:length(sync_channels)
%             sync_events = abs(sync_data) == (sync_channels(curr_sync));
%             open_ephys_ttl(curr_sync).timestamps = sync_timestamps(sync_events);
%             open_ephys_ttl(curr_sync).values = sign(sync_data(sync_events)) == 1;
%         end
%         
%         sync_save_filename = [curr_save_path filesep 'open_ephys_ttl.mat'];
%         save(sync_save_filename,'open_ephys_ttl');
        
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
        
        % Set up python
        % (add pykilosort environment paths to windows system path)
        pre_pykilosort_syspath = getenv('PATH');
        py_env_paths = {
            fullfile(char(pykilosort_python_environment.Home),'Library','bin'); ...
            fullfile(char(pykilosort_python_environment.Home),'Scripts')};
        run_pykilosort_syspath = strjoin( ...
            unique(cat(1,py_env_paths, ...
            split(pre_pykilosort_syspath,pathsep)),'stable'),pathsep);
        setenv('PATH',run_pykilosort_syspath);

        % Run pykilosort (does common average referencing by default)
        pykilosort_output_path = fullfile(local_kilosort_path,'pykilosort');
        pyrunfile('AP_run_pykilosort.py', ...
            data_filename = apband_local_filename, ...
            pykilosort_output_path = pykilosort_output_path);

        % Revert system paths to pre-pykilosort
        % (just in case alternate python environments used elsewhere)
        setenv('PATH',pre_pykilosort_syspath);

        % Grab the path with results (pykilosort makes this)
        pykilosort_results_path = fullfile(pykilosort_output_path,'output');

        % Delete TSV files (KSLabel, group, ContamPct, Amplitude), these
        % cluster groups are not used and would otherwise be loaded by
        % default into Phy
        delete(fullfile(pykilosort_results_path,'*.tsv'))

        %% Convert spike times to Open Ephys timestamps
        % This has two advantages: 
        % - sync timestamps can be used as-is without applying offset
        % - compensates for dropped samples

        % Load AP-band timestamps
        openephys_ap_timestamps_filename = fullfile(fileparts(ap_data_filename),'timestamps.npy');
        openephys_ap_timestamps = readNPY(openephys_ap_timestamps_filename);

        % Convert kilosort spike time ouput (as sample index) into open
        % ephys timestamp (in seconds)
        spike_times_kilosort_filename = fullfile(pykilosort_results_path,'spike_times.npy');
        spike_times_kilosort = readNPY(spike_times_kilosort_filename);

        % NOTE: sometimes kilsort outputs indicies of spike times which are
        % past the length of the recording??! Give a warning and calculate
        % those times with the sample rate
        spike_times_kilosort_validtime = spike_times_kilosort <= length(openephys_ap_timestamps);
        if all(spike_times_kilosort_validtime)
            spike_times_openephys = openephys_ap_timestamps(spike_times_kilosort);
        else
            kilosort_times_overshoot = double(max(spike_times_kilosort) - ...
                length(openephys_ap_timestamps))/ap_sample_rate;
            warning('Kilosort %s %s: Kilosort spike times exceed data length by %.2fs', ....
                animal,day,kilosort_times_overshoot);
            spike_times_openephys = nan(size(spike_times_kilosort));
            spike_times_openephys(spike_times_kilosort_validtime) = ...
                openephys_ap_timestamps(spike_times_kilosort(spike_times_kilosort_validtime));
            % (get samples from end, add to last timestamp - can't use
            % offset from first timestamp if timestamps aren't consistent)
            spike_times_openephys(~spike_times_kilosort_validtime) = ...
                double(spike_times_kilosort(~spike_times_kilosort_validtime) - ...
                length(openephys_ap_timestamps))/ap_sample_rate + ...
                openephys_ap_timestamps(end);
        end

        % Save open ephys spike times into kilosort output folder
        spike_times_openephys_filename = fullfile(pykilosort_results_path,'spike_times_openephys.npy');
        writeNPY(spike_times_openephys,spike_times_openephys_filename);
        
        %% Copy kilosort results to server
                
        disp('Copying sorted data to server...');
        copyfile(pykilosort_results_path,curr_save_path);
        
        %% Delete all temporary local data
        
        rmdir(local_kilosort_path,'s');
        mkdir(local_kilosort_path);
        
    end
    
end

%% Run Bombcell
ap.run_bombcell(animal,day);

%% Print end message
fprintf('\nDone preprocessing Neuropixels: %s %s\n',animal,day);


