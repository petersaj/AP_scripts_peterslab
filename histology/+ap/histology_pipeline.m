% (First: align/preprocess with AP_histology)
%
% This has become just a sandbox for histology

%% Match trajectories with days
% (plot NTE/histology days, select corresponding days from list)
%
% RUN THIS LINE BY LINE

% Need to plot new areas - preferably do day match and z align on same
% step?

% Set animal and plot probe positions
animal = 'AP009';
ap.plot_probe_positions(animal);

% Day selection for each histology trajectory
recordings = plab.find_recordings(animal);
ephys_days = {recordings([recordings.ephys]).day};
probe_ccf_dir = dir(plab.locations.filename('server',animal,[],[], ...
    'histology','*','probe_ccf.mat'));
probe_ccf_filename = fullfile(probe_ccf_dir.folder,probe_ccf_dir.name);
load(probe_ccf_filename);
for curr_probe = 1:length(probe_ccf)
    [day_idx,tf] = listdlg('PromptString',sprintf('Select day: Trajectory %d',curr_probe), ...
        'ListString',ephys_days,'SelectionMode','single');
    probe_ccf(curr_probe).day = ephys_days{day_idx};
end

% Save
save(probe_ccf_filename,'probe_ccf');
disp(['Saved ' probe_ccf_filename]);

%% Align histology depth to recording

animal = 'AP009';

probe_ccf_dir = dir(plab.locations.filename('server',animal,[],[], ...
    'histology','*','probe_ccf.mat'));
probe_ccf_filename = fullfile(probe_ccf_dir.folder,probe_ccf_dir.name);
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
    rec_time = recordings.recording{end};
    try
        ap.load_recording;
    catch me
        continue
    end

    % Plot trajectory areas
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

    % Interpolate depth from lowest DV point, ylim by depth estimate
    deepest_marked_depth = interp1(probe_ccf(curr_probe).trajectory_coords(:,2), ...
        probe_ccf(curr_probe).trajectory_areas.trajectory_depth([1,end]), ...
        max(probe_ccf(curr_probe).points(:,2)));
    ylim([deepest_marked_depth-3840,deepest_marked_depth]);

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

% Save probe depth
for curr_probe = 1:length(probe_ccf)
    curr_probe_depth = round(ylim(probe_ax(curr_probe)));
    probe_ccf(curr_probe).trajectory_areas.probe_depth = ...
        probe_ccf(curr_probe).trajectory_areas.trajectory_depth - ...
        curr_probe_depth(1);
end
save(probe_ccf_filename,'probe_ccf');
fprintf('Saved %s\n',probe_ccf_filename);


%% Align depth to recording for one day? 

animal = 'AP009';
rec_day = '2023-07-12';
load_parts.ephys = true;
ap.load_recording;

probe_ccf_dir = dir(plab.locations.filename('server',animal,[],[], ...
    'histology','*','probe_ccf.mat'));
probe_ccf_filename = fullfile(probe_ccf_dir.folder,probe_ccf_dir.name);
load(probe_ccf_filename);

curr_probe = 1;

% Get trajectory areas
trajectory_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
    probe_ccf(curr_probe).trajectory_areas.color_hex_triplet,'uni',false)),[1,3,2]);

trajectory_areas_boundaries = probe_ccf(curr_probe).trajectory_areas.trajectory_depth;
trajectory_areas_centers = mean(trajectory_areas_boundaries,2);

trajectory_areas_image_depth = 0:1:max(trajectory_areas_boundaries,[],'all');
trajectory_areas_image_idx = interp1(trajectory_areas_boundaries(:,1), ...
    1:height(probe_ccf(curr_probe).trajectory_areas),trajectory_areas_image_depth, ...
    'previous','extrap');
trajectory_areas_image = trajectory_areas_rgb(trajectory_areas_image_idx,:,:);

% Interpolate depth from lowest DV point, ylim by depth estimate
deepest_marked_depth = interp1(probe_ccf(curr_probe).trajectory_coords(:,2), ...
    probe_ccf(curr_probe).trajectory_areas.trajectory_depth([1,end]), ...
    max(probe_ccf(curr_probe).points(:,2)));

% Get spikes normalized rate by depth
norm_template_spike_n = mat2gray(log10(accumarray(spike_templates,1)+1));

% Get correlation of multiunit in sliding windows
depth_corr_window = 50; % MUA window in microns
depth_corr_window_spacing = 20; % MUA window spacing in microns

max_depths = 3840; % (hardcode, sometimes kilosort drops channels)

depth_corr_bins = [0:depth_corr_window_spacing:(max_depths-depth_corr_window); ...
    (0:depth_corr_window_spacing:(max_depths-depth_corr_window))+depth_corr_window];
depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

spike_binning_t = 0.05; % seconds
spike_binning_t_edges = nanmin(spike_times_timelite):spike_binning_t:nanmax(spike_times_timelite);

binned_spikes_depth = zeros(size(depth_corr_bins,2),length(spike_binning_t_edges)-1);
for curr_depth = 1:size(depth_corr_bins,2)
    curr_depth_templates_idx = ...
        find(template_depths >= depth_corr_bins(1,curr_depth) & ...
        template_depths < depth_corr_bins(2,curr_depth));
    
    binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timelite( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

mua_corr = corrcoef(binned_spikes_depth');

% Make figure
figure;
unit_ax = axes('position',[0.1,0,0.14,0.95]); axis off;
mua_corr_ax = axes('position',[0.26,0,0.74,0.95]); axis off;
area_ax = axes('position',[0.1,0,0.14,0.95],'color','none');

scatter(unit_ax,norm_template_spike_n,template_depths,15,'k','filled');
imagesc(mua_corr_ax,depth_corr_bin_centers,depth_corr_bin_centers,mua_corr)
area_image = image(area_ax,trajectory_areas_image_depth,[],trajectory_areas_image);
area_image.AlphaData = 0.5;
set(unit_ax,'YDir','reverse');
set(area_ax,'color','none');
axis([unit_ax,mua_corr_ax],'off')

linkaxes([unit_ax,mua_corr_ax],'y');
ylim([unit_ax,mua_corr_ax,area_ax],[0,3840]);
ylim(area_ax,[deepest_marked_depth-3840,deepest_marked_depth]);

colormap(mua_corr_ax,ap.colormap('BWR'))
clim(mua_corr_ax,[-0.5,0.5]);

yline(area_ax,unique(trajectory_areas_boundaries(:)),'color','k','linewidth',1);
set(area_ax(curr_probe),'XTick',[],'YTick',trajectory_areas_centers, ...
    'YTickLabels',probe_ccf(curr_probe).trajectory_areas.acronym);

title(unit_ax,'Unit depth x rate');
title(mua_corr_ax,'Multiunit correlation');

%% Miscellaneous 

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

%% Testing: depth alignment tool 

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
st = loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));

slice_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Data\AP010\histology';

rec_day = '2023-08-24';
animal = 'AP010';
load_parts.ephys = true;
ap.load_recording;

use_probe = 1;
AP_align_probe_histology(st,slice_path, ...
    spike_times_timelite,spike_templates,template_depths,1);


%% Testing: annotator saved coords

annotation_vertices_ccf = ...
    [vertcat(AP_histology_processing.annotation.vertices_ccf.ap), ...
    vertcat(AP_histology_processing.annotation.vertices_ccf.dv), ...
    vertcat(AP_histology_processing.annotation.vertices_ccf.ml)];

% Create axis
figure('color','w');
ccf_ax = axes;
set(ccf_ax,'ZDir','reverse');
axis(ccf_ax,'vis3d','equal','manual','tight');
hold(ccf_ax,'on');
view(ccf_ax,[-30,25]);
h_rot = rotate3d(ccf_ax);
h_rot.Enable = 'on';

% Plot 3D brain outlines
[av,tv,gui_data.st] = ap_histology.load_ccf;

slice_spacing = 5;
brain_volume = ...
    bwmorph3(bwmorph3(av(1:slice_spacing:end, ...
    1:slice_spacing:end,1:slice_spacing:end)>1,'majority'),'majority');
brain_outline_patchdata = isosurface(permute(brain_volume,[2,1,3]),0.5);
brain_outline = patch(ccf_ax, ...
    'Vertices',brain_outline_patchdata.vertices*slice_spacing, ...
    'Faces',brain_outline_patchdata.faces, ...
    'FaceColor',[0.7,0.7,0.7],'EdgeColor','none','FaceAlpha',0.1);

% Plot annotation points
plot3(annotation_vertices_ccf(:,1),annotation_vertices_ccf(:,2),annotation_vertices_ccf(:,3),'.r','MarkerSize',20)



%% Testing: depth alignment tool v2

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
st = loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));

animal = 'AP009';
histology_path = fullfile(plab.locations.server_path,'Data',animal,'histology');

rec_day = '2023-06-26';
rec_time = '1603';

load_parts.ephys = true;
verbose = true;
ap.load_recording;

use_probe = 1;
ap.histology_ephys_align(st,histology_path, ...
    spike_times_timelite,spike_templates,template_depths,1);


%% Grab areas along trajectory from probe annotation

slice_path = 'C:\Users\peter\OneDrive\Desktop\test_hist';

% Load probe CCF
histology_processing_fn = fullfile(slice_path,'AP_histology_processing.mat');
load(histology_processing_fn);

annotation_vertices_ccf = ...
    [vertcat(AP_histology_processing.annotation.vertices_ccf.ap), ...
    vertcat(AP_histology_processing.annotation.vertices_ccf.dv), ...
    vertcat(AP_histology_processing.annotation.vertices_ccf.ml)];

% Load CCF atlas
[av,tv,st] = ap_histology.load_ccf;

% Get line of best fit through mean of marked points
[~,~,probe_fit] = svd(annotation_vertices_ccf - mean(annotation_vertices_ccf,1),0);

probe_direction = probe_fit(:,1);
% (ensure vector goes downward in DV)
probe_direction(2) = abs(probe_direction(2));
probe_vector = mean(annotation_vertices_ccf,1) + padarray(probe_direction',[1,0],0,'pre');

% Grab AV values across probe trajectory
max_eval = round(sqrt(max(size(av)).^2*2));
eval_points_probe = (-max_eval:max_eval);
eval_points_ccf = round(interp1([0,norm(diff(probe_vector))],probe_vector, ...
    eval_points_probe,'linear','extrap'));

eval_points_ccf_valid = eval_points_ccf(all(eval_points_ccf > 0 & eval_points_ccf <= size(av),2),:);

eval_points_ccf_idx = sub2ind(size(av),eval_points_ccf_valid(:,1), ...
    eval_points_ccf_valid(:,2),eval_points_ccf_valid(:,3));

probe_trajectory_av = av(eval_points_ccf_idx);


% Get trajectory area colors and plot
trajectory_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
    st(probe_trajectory_av,:).color_hex_triplet,'uni',false)),[1,3,2]);
figure;image(trajectory_areas_rgb)


% Maybe save as probe_ccf file from before, make loadable with load_ephys?
% probe_ccf. 
% points - (not needed now? raw points)
% trajectory_coords - CCF [n x 3]
% trajectory areas - AV [n x 1]
probe_ccf = struct('trajectory_coords',eval_points_ccf, ...
    'trajectory_areas',probe_trajectory_av);

% Draw points
figure('color','w');
ccf_ax = axes;
set(ccf_ax,'ZDir','reverse');
axis(ccf_ax,'vis3d','equal','manual','tight');
hold(ccf_ax,'on');
view(ccf_ax,[-30,25]);
h_rot = rotate3d(ccf_ax);
h_rot.Enable = 'on';

[av,tv,gui_data.st] = ap_histology.load_ccf;

slice_spacing = 5;
brain_volume = ...
    bwmorph3(bwmorph3(av(1:slice_spacing:end, ...
    1:slice_spacing:end,1:slice_spacing:end)>1,'majority'),'majority');
brain_outline_patchdata = isosurface(permute(brain_volume,[2,1,3]),0.5);
brain_outline = patch(ccf_ax, ...
    'Vertices',brain_outline_patchdata.vertices*slice_spacing, ...
    'Faces',brain_outline_patchdata.faces, ...
    'FaceColor',[0.7,0.7,0.7],'EdgeColor','none','FaceAlpha',0.1);

plot3(annotation_vertices_ccf(:,1),annotation_vertices_ccf(:,2),annotation_vertices_ccf(:,3),'.r','MarkerSize',20)
line(eval_points_ccf_valid([1,end],1), ...
    eval_points_ccf_valid([1,end],2), ...
    eval_points_ccf_valid([1,end],3),'color','b','linewidth',2)




