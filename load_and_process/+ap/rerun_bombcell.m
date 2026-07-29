function rerun_bombcell(animal,rec_day)
% rerun_bombcell(animal,rec_day)
%
% Re-runs bombcell (in case parameters changed)
%
% INPUTS
% animal - animal name
% rec_day (optional) - day to re-bombcell (if not entered, will
% re-bombcell all recordings for that animal)

if ~exist('rec_day','var') || isempty(rec_day)
    recordings = plab.find_recordings(animal);
    rec_day = {recordings([recordings.ephys] > 0).day};
end

if ~iscell(rec_day)
    rec_day = cellstr(rec_day);
end

for curr_rec = 1:length(rec_day)

    curr_rec_day = rec_day{curr_rec};
    fprintf('Re-running Bombcell: %s %s...\n',animal,curr_rec_day)

    % Get paths and filenames
    % (ephys parent)
    ephys_path = plab.locations.filename('server',animal,curr_rec_day,[],'ephys');
    % (kilosort)
    kilosort_top_path = fullfile(ephys_path,'kilosort4');
    % (OE recordings)
    open_ephys_data_dir = dir(fullfile(ephys_path,'/**/','continuous.dat'));
    open_ephys_recording_paths = string({open_ephys_data_dir(~contains({open_ephys_data_dir.folder}, ...
        {'ADC','LFP'},'IgnoreCase',true)).folder});
    % (OE metadata)
    oe_metadata_fn = fullfile(fileparts(fileparts(open_ephys_recording_paths(1))),'structure.oebin');

    % Loop through probes, run bombcell
    for curr_probe = 1:length(open_ephys_recording_paths)

        oe_recordings = unique(extract(open_ephys_recording_paths,'recording'+digitsPattern));
        if ~isscalar(open_ephys_recording_paths) && length(open_ephys_recording_paths) == length(oe_recordings)
            % Multiple recordings: assume concatenated, use all folders
            open_ephys_path = open_ephys_recording_paths;
        else
            % Single recording or muliple non-recording folders: separate probes
            open_ephys_path = open_ephys_recording_paths(curr_probe);
        end

        open_ephys_filename = fullfile(open_ephys_path,'continuous.dat');

        % Get kilosort output folder
        kilosort_dir = dir(kilosort_top_path);
        if any(contains({kilosort_dir([kilosort_dir.isdir]).name},'probe'))
            % Multi-probe: nested folder kilosort/probe_n
            kilosort_path = fullfile(kilosort_top_path,sprintf('probe_%d',curr_probe));
        else
            % Single probe: data in kilosort top path
            kilosort_path = kilosort_top_path;
        end

        % Re-run bombcell
        kilosortVersion = 4;
        rerun = true;
        ap.run_bombcell(open_ephys_filename,kilosort_path,oe_metadata_fn,kilosortVersion,rerun);

    end
end