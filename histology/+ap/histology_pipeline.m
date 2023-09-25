%% Load CCF and set paths for slide and slice images

% Load CCF atlas
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

%% Get animal slice path

% Set paths for histology images and directory to save slice/alignment
animal = 'AM005';

im_path = sprintf('P:\\Data\\%s\\histology',animal);
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

%%%%%%% TESTING: 

% Get days with ephys
recordings = plab.find_recordings(animal);
ephys_days = {recordings([recordings.ephys]).day};

% Match trajectories with days
probe_ccf_filename = fullfile(slice_path,'probe_ccf.mat');
load(probe_ccf_filename);
for curr_probe = 1:length(probe_ccf)
    [day_idx,tf] = listdlg('PromptString',sprintf('Select day: Trajectory %d',curr_probe), ...
        'ListString',ephys_days,'SelectionMode','single');
    probe_ccf(curr_probe).day = ephys_days{day_idx};
end
save(probe_ccf_filename,'probe_ccf');

% Align histology to electrophysiology
probe_ccf_filename = fullfile(slice_path,'probe_ccf.mat');
load(probe_ccf_filename);
for curr_probe = 1:length(probe_ccf)
    % Load first recording of the day
    recordings = plab.find_recordings(animal,probe_ccf(curr_probe).day);
    load_parts.ephys = true;
    rec_day = recordings.day;
    rec_time = recordings.recording{1};
    ap.load_recording;

    % Manually align depth
    gui_fig = AP_align_probe_histology(st,slice_path, ...
        spike_times_timeline,spike_templates,template_depths,curr_probe);
    waitfor(gui_fig);
end


%%%%%% TRY A NEW VERSION: just plot units/rate with areas in background,
%%%%%% slide up and down

probe_ccf_filename = fullfile(slice_path,'probe_ccf.mat');
load(probe_ccf_filename);

figure('Name','Trajectory areas');
tiledlayout(1,length(probe_ccf));
probe_ax = gobjects(length(probe_ccf),1);

for curr_probe = 1:length(probe_ccf)

    % Load first recording of the day
    recordings = plab.find_recordings(animal,probe_ccf(curr_probe).day);
    load_parts.ephys = true;
    rec_day = recordings.day;
    rec_time = recordings.recording{1};
    ap.load_recording;

    probe_ax(curr_probe) = nexttile;    

    % Plot 
    trajectory_area_boundaries = ...
        [1;find(diff(probe_ccf(curr_probe).trajectory_areas) ~= 0)+1; ...
        length(probe_ccf(curr_probe).trajectory_areas)];    
    trajectory_area_centers_um = (trajectory_area_boundaries(1:end-1) + ...
        diff(trajectory_area_boundaries)/2)*10;
    trajectory_area_labels = ...
        st.acronym(probe_ccf(curr_probe).trajectory_areas(trajectory_area_boundaries(1:end-1)));

    probe_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
        st.color_hex_triplet(probe_ccf(curr_probe).trajectory_areas),'uni',false)),[1,3,2]);

    yyaxis left;
    image([],(1:size(probe_areas_rgb,1))*10,probe_areas_rgb);
    set(gca,'YTick',trajectory_area_centers_um,'YTickLabels',trajectory_area_labels);
    set(gca,'XTick',[]);    
    yline(-4000); % (draw a line to allow panning beyond data - hacky)
    ylim([0,3840]);

    % Plot spikes normalized rate by depth
    spike_templates_unique = unique(spike_templates);
    norm_template_spike_n = mat2gray(log10(accumarray(spike_templates,1)+1));
    yyaxis right; set(gca,'YDir','reverse'); hold on;
    ylabel('Depth (/mum)'); xlabel('Normalized spike rate');
    scatter(norm_template_spike_n(spike_templates_unique), ...
        template_depths(spike_templates_unique),15,'k','filled');
    ylim([0,3840]);

    title(rec_day);
    pan('yon');
end

% Get the probe depths corresponding to the trajectory areas
trajectory_offset = cellfun(@(x) x(1),ylim(probe_ax));

probe_depths = cellfun(@(coord,offset) ...
    pdist2(coord,coord(1,:))*10 + offset, ...
    {probe_ccf.trajectory_coords},num2cell(trajectory_offset)','uni',false);

[probe_ccf.probe_depths] = probe_depths{:};

save(probe_ccf_filename,'probe_ccf');


%%%%%%%


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

% Plot probe areas
figure('Name','Trajectory areas');
tiledlayout(1,length(probe_ccf));
% (load the colormap - located in the repository, find by associated fcn)
allenCCF_path = fileparts(which('allenCCFbregma'));
cmap_filename = [allenCCF_path filesep 'allen_ccf_colormap_2017.mat'];
load(cmap_filename);
colormap(cmap);

for curr_probe = 1:length(probe_ccf)
    nexttile;    
    trajectory_area_boundaries = ...
        [1;find(diff(probe_ccf(curr_probe).trajectory_areas) ~= 0);length(probe_ccf(curr_probe).trajectory_areas)];    
    trajectory_area_centers = trajectory_area_boundaries(1:end-1) + diff(trajectory_area_boundaries)/2;
    trajectory_area_labels = st.acronym(probe_ccf(curr_probe).trajectory_areas(round(trajectory_area_centers)));
      
    image(probe_ccf(curr_probe).trajectory_areas);
    
    clim([1,size(cmap,1)])
    set(gca,'YTick',trajectory_area_centers,'YTickLabels',trajectory_area_labels);
    set(gca,'XTick',[]);
    title(['Probe ' num2str(curr_probe)]);
    
end














