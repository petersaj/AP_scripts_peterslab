function plot_probe_positions(animal,plot_units,plot_nte,plot_histology)
% plot_probe_positions(animal)
%
% Plot probe positions from Trajectory Explorer and histology

arguments
    animal = [];
    plot_units logical = true;
    plot_nte logical = true;
    plot_histology logical = true;
end

%% Find probe position files (trajectory explorer and histology)

% Trajectory explorer probe position files
% (include /ephys and /ephys/** subfolders)
nte_filepattern = plab.locations.filename('server',animal,'*',[],'ephys','**','*probe_positions.mat');
nte_filenames = dir(nte_filepattern);

% Histology probe position files
% % (old version)
% histology_filepattern = plab.locations.filename('server',animal,[],[],'histology','*','probe_ccf.mat');
% histology_fn = dir(histology_filepattern);
% (new version)
histology_filepattern = plab.locations.filename('server',animal,[],[],'histology','**','AP_histology_processing.mat');
histology_dir = dir(histology_filepattern);
histology_filename = fullfile(histology_dir.folder,histology_dir.name);

% Return if no files found
if isempty(nte_filenames) && isempty(histology_filename)
    fprintf('No probe positions found: %s\n',animal);
    return;
end


%% Plot brain outline

ccf_outline = ap.ccf_outline_3d;
ccf_outline.Parent.Parent.Name = animal;
ccf_3d_axes = ccf_outline.Parent;

%% Draw probes (from Neuropixels Trajectory Explorer)

if plot_nte

    for curr_recording = 1:length(nte_filenames)

        load(fullfile(nte_filenames(curr_recording).folder, nte_filenames(curr_recording).name));

        % Loop through probes and draw
        for curr_probe = 1:length(probe_positions_ccf)
            % Draw probe in 3D view
            line(ccf_3d_axes, ...
                probe_positions_ccf{curr_probe}(1,:), ...
                probe_positions_ccf{curr_probe}(3,:), ...
                probe_positions_ccf{curr_probe}(2,:), ...
                'linewidth',2,'color','b')

            date_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
            curr_day = cell2mat(extract(nte_filenames(curr_recording).folder,date_pattern));
            text(probe_positions_ccf{curr_probe}(1,1), ...
                probe_positions_ccf{curr_probe}(3,1), ...
                probe_positions_ccf{curr_probe}(2,1), ...
                curr_day,'color','b');
            text(probe_positions_ccf{curr_probe}(1,2), ...
                probe_positions_ccf{curr_probe}(3,2), ...
                probe_positions_ccf{curr_probe}(2,2), ...
                curr_day,'color','b');

        end
    end
end

%% Draw probes (from histology)

if plot_histology
    if ~isempty(histology_filename)

        % Load histology and get probe fits
        load(histology_filename);
        probe_line_fits = ap_histology.fit_probe_line(histology_filename);

        for curr_probe = 1:length(probe_line_fits)
            % Draw probe line
            line(ccf_3d_axes, ...
                probe_line_fits(curr_probe).ccf(:,1), ...
                probe_line_fits(curr_probe).ccf(:,3), ...
                probe_line_fits(curr_probe).ccf(:,2), ...
                'linewidth',2,'color','r')

            % Label probe
            text(probe_line_fits(curr_probe).ccf(1,1), ...
                probe_line_fits(curr_probe).ccf(1,3), ...
                probe_line_fits(curr_probe).ccf(1,2), ...
                probe_line_fits(curr_probe).label,'color','r');
        end
    end
end

drawnow

%% Plot areas and unit positions/rate

if plot_units

    disp('Finding and plotting units...')

    recordings = plab.find_recordings(animal);
    ephys_recordings = recordings([recordings.ephys]>0);

    figure('color','w');
    h = tiledlayout(1,length(ephys_recordings));
    title(h,animal);
    for curr_recording = 1:length(ephys_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = ephys_recordings(curr_recording).day;
        rec_time = ephys_recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        unit_axes = nexttile(h); hold on;
        ap.plot_unit_depthrate(spike_times_timelite,spike_templates,template_depths,probe_areas,unit_axes);

        drawnow;

        ap.print_progress_fraction(curr_recording,length(ephys_recordings));

    end

end











