function where2slice(animal)
% where2slice(animal)
%
% Plot the anterior- and posterior-most CCF slices with probe to give a
% ballpark of where to slice for histology


%% Find probe position files (trajectory explorer and histology)

% Trajectory explorer probe position files
% (include /ephys and /ephys/** subfolders)
nte_filepattern = plab.locations.filename('server',animal,'*',[],'ephys','**','*probe_positions.mat');
nte_fns = dir(nte_filepattern);

% Return if no files found
if isempty(nte_fns) && isempty(histology_fn)
    fprintf('No probe positions found: %s\n',animal);
    return;
end

%% Load atlas

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"

%% Set up figure

figure('color','w','name',animal);colormap gray;
tiledlayout(1,3,'TileSpacing','compact');

% Set up 3D axes
ccf_3d_axes = nexttile;
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

% Draw probes and grab positions
probe_positions_ccf_all = cell(0);
for curr_recording = 1:length(nte_fns)

    load(fullfile(nte_fns(curr_recording).folder, nte_fns(curr_recording).name));
    probe_positions_ccf_all = vertcat(probe_positions_ccf_all,probe_positions_ccf);

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

%% Plot CCF slices for anterior- and posterior-most points

probe_positions_horzcat = cat(2,probe_positions_ccf_all{:});
ap_probe_extent = round(prctile(probe_positions_horzcat(1,:),[0,100]));

c = [0,300];

nexttile;
imagesc(permute(tv(ap_probe_extent(1),:,:),[2,3,1]));
axis image off; clim(c);
title('Anterior start');

nexttile;
imagesc(permute(tv(ap_probe_extent(2),:,:),[2,3,1]));
axis image off; clim(c);
title('Posterior end');




























