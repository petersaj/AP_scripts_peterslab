function preprocess_neuropixels(animal,day)
% preprocess_neuropixels(animal,day)
%
% Spike sorts with Pykilosort (into 'pykilosort' folder)
%
% Currently assumes Neuropixels Phase 3A recorded with Open Ephys

%% Get paths and filenames

ephys_path = plab.locations.make_server_filename(animal,day,[],'ephys');

ephys_exists = exist(ephys_path,'dir');

if ~ephys_exists
    error([animal ' ' day ': No ephys data found']);
end

save_paths = {[ephys_path filesep 'pykilosort']};
data_paths = {ephys_path};

% Check for multiple sites (assume sites are marked as site#)
data_path_dir = dir([data_paths{1} filesep 'probe_*']);
if ~isempty(data_path_dir)
    data_paths = cellfun(@(x) [data_paths{1} filesep x],{data_path_dir.name},'uni',false);
    save_paths = cellfun(@(x) [save_paths{1} filesep x],{data_path_dir.name},'uni',false);
end

for curr_site = 1:length(data_paths)  
    
    % Get experiments (if turned off between)
    curr_data_path = data_paths{curr_site};   
    ephys_exp_paths = dir([curr_data_path filesep 'experiment*']);
       
    for curr_exp = 1:length(ephys_exp_paths)
        
        % Update save path with experiment (only if more than one, legacy)
        if length(ephys_exp_paths) == 1
            curr_save_path = save_paths{curr_site};
        elseif length(ephys_exp_paths) > 1
            curr_save_path = [save_paths{curr_site} filesep ephys_exp_paths(curr_exp).name];
        end
        
        % Make save path
        if ~exist(curr_save_path,'dir')
            mkdir(curr_save_path)
        end
              
        % Get OE filenames (check for multiple experiments, do separately)
        exp_rec_dir = [ephys_exp_paths(curr_exp).name filesep 'recording1'];
        ap_data_filename = [curr_data_path filesep exp_rec_dir filesep 'continuous' filesep 'Neuropix-3a-100.Neuropix-3a-AP' filesep 'continuous.dat'];
        lfp_data_filename = [curr_data_path filesep exp_rec_dir filesep 'continuous' filesep 'Neuropix-3a-100.Neuropix-3a-LFP' filesep 'continuous.dat'];
        sync_filename = [curr_data_path filesep filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.Neuropix-3a-AP' filesep 'TTL' filesep 'states.npy' ];
        sync_timestamps_filename = [curr_data_path filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.Neuropix-3a-AP' filesep 'TTL' filesep 'timestamps.npy' ];
        messages_filename = [curr_data_path filesep exp_rec_dir filesep 'sync_messages.txt'];
        settings_filename = [curr_data_path filesep exp_rec_dir filesep 'structure.oebin'];
        
        
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% 

%         Currently obsolete: loaded in loading script

        %%%%%%%%%%%%%%%%%%%%%%%%%%%


        % NOTE: open ephys gives timestamps relative to start of PREVIEW,
        % kilosort only gives timestamps relative to beginning of RECORD.
        % This means that to match up kilosort (t(1) = 0) to OE (t(1) =
        % arbitrary), the sync timestamps are saved as sync_timestamps -
        % ap_timestamps(1)
              
        % Load spike timestamps
        ap_timestamps_filename = fullfile(fileparts(ap_data_filename),'timestamps.npy');
        ap_timestamps = readNPY(ap_timestamps_filename);
        sync_timestamp_offset = ap_timestamps(1);
                    
        % Load digital input event times
        sync_data = readNPY(sync_filename);
        sync_timestamps = readNPY(sync_timestamps_filename);
        
        sync_channels = unique(abs(sync_data));
        sync = struct('timestamps',cell(size(sync_channels)),'values',cell(size(sync_channels)));
        for curr_sync = 1:length(sync_channels)
            sync_events = abs(sync_data) == (sync_channels(curr_sync));
            sync(curr_sync).timestamps = sync_timestamps(sync_events) - sync_timestamp_offset;
            sync(curr_sync).values = sign(sync_data(sync_events)) == 1;
        end
        
        sync_save_filename = [curr_save_path filesep 'sync.mat'];
        save(sync_save_filename,'sync');


%         %%%%%%%% TESTING: use sample number instead of timestamp
%         ap_sample_number_filename = fullfile(fileparts(ap_data_filename),'sample_numbers.npy');
%         ap_sample_number = readNPY(ap_sample_number_filename);
%         sync_sample_number_offset = ap_sample_number(1);
%         
%         % Load digital input event times
%         sync_data = readNPY(sync_filename);
%         sync_sample_numbers_filename = [curr_data_path filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.Neuropix-3a-AP' filesep 'TTL' filesep 'sample_numbers.npy' ];
%         sync_sample_numbers = readNPY(sync_sample_numbers_filename);
% 
%         sync_channels = unique(abs(sync_data));
%         sync = struct('timestamps',cell(size(sync_channels)),'values',cell(size(sync_channels)));
%         for curr_sync = 1:length(sync_channels)
%             sync_events = abs(sync_data) == (sync_channels(curr_sync));
%             sync(curr_sync).timestamps = double((sync_sample_numbers(sync_events) - sync_sample_number_offset))/30000;
%             sync(curr_sync).values = sign(sync_data(sync_events)) == 1;
%         end
% 
%         sync_save_filename = [curr_save_path filesep 'sync.mat'];
%         save(sync_save_filename,'sync');


        %% Run kilosort
        
        % Set up local directory and clear out
        ssd_kilosort_path = 'D:\data_temp\kilosort';
        
        % Clear out local kilosort directories
        if exist(ssd_kilosort_path,'dir')
            rmdir(ssd_kilosort_path,'s')
        end
        mkdir(ssd_kilosort_path);
        
        % Copy AP-band data locally
        disp('Copying AP data to local drive...')
        apband_local_filename = [ssd_kilosort_path filesep animal '_' day  '_' 'ephys_apband.dat'];
        copyfile(ap_data_filename,apband_local_filename);
        disp('Done');
        
        % Set up python
        % (set pykilosort python environment)
        py_version = pyenv('Version','C:\Users\petersa\anaconda3\envs\pyks2\pythonw.exe');
        % (add pykilosort environment paths to windows system path)
        pre_pykilosort_syspath = getenv('PATH');
        py_env_paths = {
            fullfile(char(py_version.Home),'Library','bin'); ...
            fullfile(char(py_version.Home),'Scripts')};
        run_pykilosort_syspath = strjoin( ...
            unique(cat(1,py_env_paths, ...
            split(pre_pykilosort_syspath,pathsep)),'stable'),pathsep);
        setenv('PATH',run_pykilosort_syspath);

        % Run pykilosort
        % (directly on raw data - pykilosort does local common average
        pykilosort_output_path = fullfile(ssd_kilosort_path,'pykilosort');
        pyrunfile('AP_run_pykilosort.py', ...
            data_filename = apband_local_filename, ...
            pykilosort_output_path = pykilosort_output_path);

        % Revert system paths to pre-pykilosort
        % (just in case alternate python environments used elsewhere)
        setenv('PATH',pre_pykilosort_syspath);

        % Delete TSV files (KSLabel, group, ContamPct, Amplitude), these
        % cluster groups are not used and would otherwise be loaded by
        % default into Phy
        delete(fullfile(pykilosort_output_path,'*.tsv'))
        
        %% Copy kilosort results to server
                
        disp('Copying sorted data to server...');

        pykilosort_results_path = fullfile(pykilosort_output_path,'output');
        copyfile(pykilosort_results_path,curr_save_path);
        
        %% Delete all temporary local data
        rmdir(ssd_kilosort_path,'s');
        mkdir(ssd_kilosort_path);
        
    end
    
end

fprintf('Done preprocessing Neuropixels: %s %s',animal,day);


