%% Load CCF and set paths for slide and slice images

% Load CCF atlas
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

%% Get animal slice path

% Set paths for histology images and directory to save slice/alignment
animal = 'AP008';

im_path = plab.locations.filename('server',animal,[],[],'histology');
slice_path = fullfile(im_path,'slices');


%% Preprocess slide images to produce slice images

% Set white balance and resize slide images, extract slice images
% (Note: this resizes the images purely for file size reasons - the CCF can
% be aligned to histology no matter what the scaling. If pixel size is
% available in metadata then automatically scales to CCF resolution,
% otherwise user can specify the resize factor as a second argument)

% Set resize factor
% resize_factor = []; % (slides ome.tiff: auto-resize ~CCF size 10um/px)
resize_factor = 1; % (slides tiff: resize factor)

% Set slide or slice images
% slice_images = false; % (images are slides - extract individual slices)
slice_images = true; % (images are already individual slices)

% Preprocess images
AP_process_histology(im_path,resize_factor,slice_images);

% (optional) Rotate, center, pad, flip slice images
AP_rotate_histology(slice_path);

%% Align CCF to slices

% Find CCF slices corresponding to each histology slice
AP_grab_histology_ccf(tv,av,st,slice_path);

% Align CCF slices and histology slices
% (first: automatically, by outline)
AP_auto_align_histology_ccf(slice_path);
% (second: curate manually)
AP_manual_align_histology_ccf(tv,av,st,slice_path);


%% Utilize aligned CCF

% Display aligned CCF over histology slices
AP_view_aligned_histology(st,slice_path);

% Display histology within 3D CCF
AP_view_aligned_histology_volume(tv,av,st,slice_path,1);

% Get probe trajectory from histology, convert to CCF coordinates
AP_get_probe_histology(tv,av,st,slice_path);

% Match trajectories with days
% (plot NTE/histology days, select corresponding days from list)
ap.plot_probe_positions(animal);

recordings = plab.find_recordings(animal);
ephys_days = {recordings([recordings.ephys]).day};
probe_ccf_filename = fullfile(slice_path,'probe_ccf.mat');
load(probe_ccf_filename);
for curr_probe = 1:length(probe_ccf)
    [day_idx,tf] = listdlg('PromptString',sprintf('Select day: Trajectory %d',curr_probe), ...
        'ListString',ephys_days,'SelectionMode','single');
    probe_ccf(curr_probe).day = ephys_days{day_idx};
end
save(probe_ccf_filename,'probe_ccf');

% Align histology depth to recording
probe_ccf_filename = fullfile(slice_path,'probe_ccf.mat');
load(probe_ccf_filename);

figure('Name',sprintf('%s: Trajectory areas',animal));
tiledlayout(1,length(probe_ccf));
probe_ax = gobjects(length(probe_ccf),1);

for curr_probe = 1:length(probe_ccf)

        probe_ax(curr_probe) = nexttile;    

    % Load first recording of the day
    recordings = plab.find_recordings(animal,probe_ccf(curr_probe).day);
    load_parts.ephys = true;
    rec_day = recordings.day;
    rec_time = recordings.recording{1};
    try
        ap.load_recording;
    catch me
        continue
    end

    % Plot tajectory areas
    trajectory_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
        probe_ccf(curr_probe).trajectory_areas.color_hex_triplet,'uni',false)),[1,3,2]);

    trajectory_areas_boundaries = probe_ccf(curr_probe).trajectory_areas.trajectory_depth;
    trajectory_areas_centers = mean(trajectory_areas_boundaries,2);

    trajectory_areas_image_depth = 0:1:max(trajectory_areas_boundaries,[],'all');
    trajectory_areas_image_idx = interp1(trajectory_areas_boundaries(:,1), ...
        1:height(probe_ccf(curr_probe).trajectory_areas),trajectory_areas_image_depth, ...
        'previous','extrap');
    trajectory_areas_image = trajectory_areas_rgb(trajectory_areas_image_idx,:,:);

    yyaxis left;
    image([0,1],trajectory_areas_image_depth,trajectory_areas_image);
    yline(unique(trajectory_areas_boundaries(:)),'color','k','linewidth',1);
    set(probe_ax(curr_probe),'XTick',[],'YTick',trajectory_areas_centers, ...
        'YTickLabels',probe_ccf(curr_probe).trajectory_areas.acronym);
    yline(-4000); % (draw a line to allow panning beyond data - hacky)
    ylim([0,3840]);    

    % Plot spikes normalized rate by depth
    spike_templates_unique = unique(spike_templates);
    norm_template_spike_n = mat2gray(log10(accumarray(spike_templates,1)+1));
    yyaxis right; set(gca,'YDir','reverse'); hold on;
    ylabel('Depth (\mum)'); xlabel('Normalized spike rate');
    scatter(norm_template_spike_n(spike_templates_unique), ...
        template_depths(spike_templates_unique),15,'k','filled');
    xlim([0,1]);
    ylim([0,3840]);

    title({sprintf('Trajectory %d',curr_probe),rec_day});
    pan('yon');

end

% Add probe depth information to the trajectory areas from alignment
for curr_probe = 1:length(probe_ccf)
    curr_probe_depth = round(ylim(probe_ax(curr_probe)));
    probe_ccf(curr_probe).trajectory_areas.probe_depth = ...
        probe_ccf(curr_probe).trajectory_areas.trajectory_depth - ...
        curr_probe_depth(1);
end
save(probe_ccf_filename,'probe_ccf');



%% Unused 

% Extract slices from full-resolution images
% (not worth it at the moment, each slice is 200 MB)
% AP_grab_fullsize_histology_slices(im_path)

% Convert points in histology images to CCF coordinates
ccf_points = AP_histology2ccf(histology_points,slice_path);
% Concatenate points and round to nearest integer coordinate
ccf_points_cat = round(cell2mat(ccf_points));
% Get indicies from subscripts
ccf_points_idx = sub2ind(size(av),ccf_points_cat(:,1),ccf_points_cat(:,2),ccf_points_cat(:,3));
% Find annotated volume (AV) values at points
ccf_points_av = av(ccf_points_idx);
% Get areas from the structure tree (ST) at given AV values
ccf_points_areas = st(ccf_points_areas,:).safe_name;


%% Plot probe areas

probe_ccf_filename = fullfile(slice_path,'probe_ccf.mat');
load(probe_ccf_filename);

figure('Name','Trajectory areas');
tiledlayout(1,length(probe_ccf),'TileSpacing','compact');
for curr_probe = 1:length(probe_ccf)
    
    curr_axes = nexttile;
    
    trajectory_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
        probe_ccf(curr_probe).trajectory_areas.color_hex_triplet,'uni',false)),[1,3,2]);

    trajectory_areas_boundaries = probe_ccf(curr_probe).trajectory_areas.trajectory_depth;
    trajectory_areas_centers = mean(trajectory_areas_boundaries,2);

    trajectory_areas_image_depth = 0:0.01:max(trajectory_areas_boundaries,[],'all');
    trajectory_areas_image_idx = interp1(trajectory_areas_boundaries(:,1), ...
        1:height(probe_ccf(curr_probe).trajectory_areas),trajectory_areas_image_depth, ...
        'previous','extrap');
    trajectory_areas_image = trajectory_areas_rgb(trajectory_areas_image_idx,:,:);

    image([],trajectory_areas_image_depth,trajectory_areas_image);
    yline(unique(trajectory_areas_boundaries(:)),'color','k','linewidth',1);
    set(curr_axes,'XTick',[],'YTick',trajectory_areas_centers, ...
        'YTickLabels',probe_ccf(curr_probe).trajectory_areas.acronym);
    set(curr_axes,'XTick',[]);
    if isfield(probe_ccf,'day')
        title({sprintf('Probe %d',curr_probe),probe_ccf(curr_probe).day});
    else
        title(sprintf('Probe %d',curr_probe));
    end
end












