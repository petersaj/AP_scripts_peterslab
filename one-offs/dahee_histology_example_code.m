%% Align histology to CCF atlas, plot fluorescence in coronal slices and 3D


%% Set path with images and AP histology processing file

histology_path = 'C:\Users\petersa\Desktop\Dahee_histology';

%% Get histology in CCF volume

histology_volume = ap_histology.histology_volume(histology_path);

%% Set channel to plot
% (in the example you gave me, red = channel 1)

plot_channel = 1;

%% Plot maximum-projection images within atlas bins
% Pseudo-color image, overlay atlas regions

% Get maximum projection within coronal atlas bins
atlas_bin_edges = 600:20:800; % Set bins for splitting/viewing atlas

n_channels = size(histology_volume,4);
histology_volume_binmax = nan([size(histology_volume,[2,3]),length(atlas_bin_edges)-1,n_channels]);

% (loop through bins, get max projection)
for curr_channel = 1:n_channels
    for curr_atlas_bin = 1:length(atlas_bin_edges)-1
        curr_ap = atlas_bin_edges(curr_atlas_bin):atlas_bin_edges(curr_atlas_bin+1);

        histology_volume_binmax(:,:,curr_atlas_bin,curr_channel) = ...
            squeeze(max(histology_volume(curr_ap,:,:,curr_channel), ....
            [],1,'omitnan'));
    end
end

% Plot binned max with CCF borders
% (load atlas)
av = ap_histology.load_ccf;
% (set min/max value to display)
histology_clim = [30000,max(histology_volume_binmax,[],'all')];
% (set color to pseudocolor plot)
plot_channel_color = [1,0,0];
% (set thickness of CCF border outline)
overlay_dilation = 1;
% (loop through binned images, plot with overlay)
figure; tiledlayout('TileSpacing','none');
for curr_atlas_bin = 1:size(histology_volume_binmax,3)
    % Get CCF borders from the middle of the bin
    plot_atlas_ap = round(mean(atlas_bin_edges(curr_atlas_bin+[0,1])));
    curr_ccf_borders = imdilate(boundarymask(permute(av(plot_atlas_ap,:,:),[2,3,1])),ones(overlay_dilation));

    % Create color image
    curr_im_colored = 1- ... % flip background black -> white
        (mat2gray( ... % range-normalize
        max(histology_volume_binmax(:,:,curr_atlas_bin,plot_channel),0,'omitnan'), ... % set NaNs to 0
        histology_clim) ... % clip to min/max values
        .*permute(1-plot_channel_color,[3,4,2,1])); % color (1-to flip black/white)

    % Plot max image with CCF borders overlaid
    curr_overlay = imoverlay(curr_im_colored,curr_ccf_borders,'k');
    nexttile;imagesc(curr_overlay);axis image off;
    drawnow;
end

%% Plot binarized fluorescence in 3D brain

% Set channel to plot and value to threshold 
histology_threshold = 50000;

% Moving median to fill small gaps (e.g. between slices)
histology_volume_movmed = movmedian(histology_volume(:,:,:,plot_channel),10,'omitmissing');

% Get thresholded histology 3D outline 
% (fast function extractIsosurface requires medical imaging toolbox)
[histology_faces,histology_vertices] = ...
    extractIsosurface(permute(histology_volume_movmed ...
    > histology_threshold,[3,1,2]),0.5);

% Plot histology patch data 
% (draw brain and PO outline)
ap.ccf_outline_3d([],[],["brain","PO"]);

face_alpha = 0.5;
patch('Vertices',histology_vertices, ...
    'Faces',histology_faces, ...
    'FaceColor','r','EdgeColor','none','FaceAlpha',face_alpha);