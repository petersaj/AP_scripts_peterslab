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
nte_fns = dir(nte_filepattern);

% Histology probe position files
histology_filepattern = plab.locations.filename('server',animal,[],[],'histology','*','probe_ccf.mat');
histology_fn = dir(histology_filepattern);

% Return if no files found
if isempty(nte_fns) && isempty(histology_fn)
    fprintf('No probe positions found: %s\n',animal);
    return;
end


%% Load atlas

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below

%% Set up axes, plot brain outline

figure('Color','w','name',animal);

% Set up 3D axes
ccf_3d_axes = axes;
set(ccf_3d_axes,'ZDir','reverse');
hold(ccf_3d_axes,'on');
axis(ccf_3d_axes,'vis3d','equal','off','manual');
view([-30,25]);
axis tight;
h = rotate3d(ccf_3d_axes);
h.Enable = 'on';

slice_spacing = 5;
brain_volume = ...
    bwmorph3(bwmorph3(av(1:slice_spacing:end, ...
    1:slice_spacing:end,1:slice_spacing:end)>1,'majority'),'majority');
brain_outline_patchdata = isosurface(permute(brain_volume,[3,1,2]),0.5);
brain_outline = patch( ...
    'Vertices',brain_outline_patchdata.vertices*slice_spacing, ...
    'Faces',brain_outline_patchdata.faces, ...
    'FaceColor',[0.7,0.7,0.7],'EdgeColor','none','FaceAlpha',0.1);


%% Draw probes (from Neuropixels Trajectory Explorer)

if plot_nte

    for curr_recording = 1:length(nte_fns)

        load(fullfile(nte_fns(curr_recording).folder, nte_fns(curr_recording).name));

        % Loop through probes and draw
        for curr_probe = 1:length(probe_positions_ccf)
            % Draw probe in 3D view
            line(ccf_3d_axes, ...
                probe_positions_ccf{curr_probe}(1,:), ...
                probe_positions_ccf{curr_probe}(3,:), ...
                probe_positions_ccf{curr_probe}(2,:), ...
                'linewidth',2,'color','b')

            date_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
            curr_day = cell2mat(extract(nte_fns(curr_recording).folder,date_pattern));
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

    if ~isempty(histology_fn)

        load(fullfile(histology_fn.folder,histology_fn.name));

        for curr_probe = 1:length(probe_ccf)

            % Get line of best fit through mean of marked points
            probe_coords_mean = mean(probe_ccf(curr_probe).points,1);
            xyz = bsxfun(@minus,probe_ccf(curr_probe).points,probe_coords_mean);
            [~,~,V] = svd(xyz,0);
            histology_probe_direction = V(:,1);

            % (make sure the direction goes down in DV - flip if it's going up)
            if histology_probe_direction(2) < 0
                histology_probe_direction = -histology_probe_direction;
            end

            % Evaluate line of best fit (length of probe to deepest point)
            [~,deepest_probe_idx] = max(probe_ccf(curr_probe).points(:,2));
            probe_deepest_point = probe_ccf(curr_probe).points(deepest_probe_idx,:);
            probe_deepest_point_com_dist = pdist2(probe_coords_mean,probe_deepest_point);
            probe_length_ccf = 3840/10; % mm / ccf voxel size

            probe_line_eval = probe_deepest_point_com_dist - [probe_length_ccf,0];
            probe_line = (probe_line_eval'.*histology_probe_direction') + probe_coords_mean;

            % Draw probe in 3D view
            line(ccf_3d_axes,probe_line(:,1),probe_line(:,3),probe_line(:,2), ...
                'linewidth',2,'color','r');

            % Write text
            if isfield(probe_ccf,'day') && ~isempty(probe_ccf(curr_probe).day)
                probe_text = sprintf('Trajectory %d [%s]',curr_probe,probe_ccf(curr_probe).day);
            else
                probe_text = sprintf('Trajectory %d',curr_probe);
            end

            text(probe_line(1,1),probe_line(1,3),probe_line(1,2), ...
                probe_text,'color','r');
            text(probe_line(2,1),probe_line(2,3),probe_line(2,2), ...
                probe_text,'color','r');

        end
    end
end

drawnow

%% Plot areas and unit positions/rate

if plot_units

    recordings = plab.find_recordings(animal);
    ephys_recordings = recordings([recordings.ephys]);

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

        unit_axes = nexttile; hold on;
        unit_axes.YDir = 'reverse';

        % Plot units (depth vs normalized rate) with areas
        if exist('probe_areas','var')
            probe_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
                probe_areas{1}.color_hex_triplet,'uni',false)),[1,3,2]);

            probe_areas_boundaries = probe_areas{1}.probe_depth;
            probe_areas_centers = mean(probe_areas_boundaries,2);

            probe_areas_image_depth = 0:1:max(probe_areas_boundaries,[],'all');
            probe_areas_image_idx = interp1(probe_areas_boundaries(:,1), ...
                1:height(probe_areas{1}),probe_areas_image_depth, ...
                'previous','extrap');
            probe_areas_image = probe_areas_rgb(probe_areas_image_idx,:,:);

            image(unit_axes,[0,1],probe_areas_image_depth,probe_areas_image);
            yline(unique(probe_areas_boundaries(:)),'color','k','linewidth',1);
            set(unit_axes,'YTick',probe_areas_centers,'YTickLabels',probe_areas{1}.acronym);
        end

        norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));

        unit_dots = scatter( ...
            norm_spike_n,template_depths(unique(spike_templates)),20,'k','filled');
        multiunit_lines = arrayfun(@(x) line(xlim,[0,0],'linewidth',2,'visible','off'),1:2);
        xlim(unit_axes,[-0.1,1]);
        ylim([-50, max(channel_positions(:,2))+50]);
        ylabel('Depth (\mum)')
        xlabel('Normalized log rate')
        title(rec_day);

        drawnow;

    end

end











