function plot_probe_positions(animal)
% plot_probe_positions(animal)
%
% Plot probe positions from Trajectory Explorer and histology

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

% Find all saved probe positions
probe_position_filepattern = plab.locations.make_server_filename(animal,'*',[],'ephys','*probe_positions.mat');
probe_position_fns = dir(probe_position_filepattern);

for curr_recording = 1:length(probe_position_fns)
    load(fullfile(probe_position_fns(curr_recording).folder, ...
        probe_position_fns(curr_recording).name));
    % Loop through probes and draw
    for curr_probe = 1:length(probe_positions_ccf)
        % Draw probe in 3D view
        line(ccf_3d_axes, ...
            probe_positions_ccf{curr_probe}(1,:), ...
            probe_positions_ccf{curr_probe}(3,:), ...
            probe_positions_ccf{curr_probe}(2,:), ...
            'linewidth',2,'color','b')

        date_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
        curr_day = cell2mat(extract(probe_position_fns(curr_recording).folder,date_pattern));
        text(probe_positions_ccf{curr_probe}(1,1), ...
            probe_positions_ccf{curr_probe}(3,1), ...
            probe_positions_ccf{curr_probe}(2,1), ...
            curr_day,'color','b');
    end
end

%% Draw probes (from histology)

% Find all saved probe positions
probe_position_filepattern = plab.locations.make_server_filename(animal,[],[],'histology','*','probe_ccf.mat');
probe_position_fn = dir(probe_position_filepattern);

load(fullfile(probe_position_fn.folder,probe_position_fn.name));

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

    text(probe_line(1,1), ...
        probe_line(1,3), ...
        probe_line(1,2), ...
        sprintf('Trajectory %d',curr_probe),'color','r');

end









