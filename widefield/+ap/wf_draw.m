function h = wf_draw(type,color,manual_bregma)
% h = wf_draw(type,color,manual_bregma)
%
% Draw CCF areas on aligned widefield image
%
% Type: 
% - scalebar
% - bregma
% - grid
% - ccf
% - point, [AP,ML] in mm
% - point_aligned, [AP,ML] in mm tranformed via CCF alignment
%
% Color: color of overlay
%
% manual_bregma: false (default) if aligned image, true to manually select

if ~exist('color','var')
    color = 'k';
end

if ~exist('manual_bregma','var')
    manual_bregma = false;
end


%% Load alignments and scale, define bregma

% Load alignment and perform CCF-widefield alignment
alignment_path = fullfile(plab.locations.server_path, ...
    'Users','Andy_Peters','widefield_alignment');
ccf_tform_fn = fullfile(alignment_path,'ccf_tform.mat');
load(ccf_tform_fn);

% Set microns per pixel scale (hard-coded)
um2pixel = 22.5; % cortexlab: 20.6

% Define CCF bregma (hard-coded)
bregma = [520,44,570];
bregma(3) = bregma(3) + 0.5;
bregma_wf = [bregma([3,1]),1]*ccf_tform.T;

% Set bregma to use (click if manual, use aligned if not)
if manual_bregma
    [bregma_offset_x,bregma_offset_y] = ginput(1);
else 
    [bregma_offset_x,bregma_offset_y] = deal(bregma_wf(1),bregma_wf(2));
end

% Load top-down cortical regions
load(fullfile(alignment_path,'dorsal_cortex_borders.mat'));

dorsal_cortex_borders_aligned_long = cellfun(@(areas) cellfun(@(coords) ...
    [fliplr(coords),ones(size(coords,1),1)]*ccf_tform.T,areas,'uni',false), ...
    dorsal_cortex_borders,'uni',false);
dorsal_cortex_borders_aligned = cellfun(@(areas) cellfun(@(coords) ...
    coords,areas,'uni',false),dorsal_cortex_borders_aligned_long,'uni',false);


%% Draw object

% Ensure image is held
hold on;

switch type   

    case 'scalebar'
        % Plot 1 mm scalebar
        y_scale_um = 1000;
        y_scale_px = y_scale_um/um2pixel;
        line(gca,repmat(min(xlim(gca)),2,1), ...
            min(ylim(gca)) + [0,y_scale_px],'color',color,'linewidth',2)

    case 'bregma'
        % Plot bregma on aligned widefield
        h.bregma = plot(bregma_offset_x,bregma_offset_y,'.','color',color,'MarkerSize',15);
    
    case 'grid'
        % Plot 500um spaced grid
                
        spacing_um = 500;
        spacing_pixels = spacing_um/um2pixel;
        
        xlines_pos = bregma_offset_y + spacing_pixels*(ceil((min(ylim)-bregma_offset_y)./spacing_pixels):floor((max(ylim)-bregma_offset_y)./spacing_pixels));
        ylines_pos = bregma_offset_x + spacing_pixels*(ceil((min(xlim)-bregma_offset_x)./spacing_pixels):floor((max(xlim)-bregma_offset_x)./spacing_pixels));
        
        h = struct;
        warning('um/px: needs more accurate measurement');

        
        for curr_xline = 1:length(xlines_pos)
            h.xlines(curr_xline) = line(xlim,repmat(xlines_pos(curr_xline),1,2),'color',color,'linestyle','-');
        end
        
        for curr_yline = 1:length(ylines_pos)
            h.ylines(curr_yline) = line(repmat(ylines_pos(curr_yline),1,2),ylim,'color',color,'linestyle','-');
        end
        
        h.bregma = plot(bregma_offset_x,bregma_offset_y,'xr','MarkerSize',10);

    case 'point'
        % Plot a point at specific coordinates (in um/px)
        [point_ap,point_ml] = deal(color(1),color(2));

        plot(bregma_offset_x+point_ml*1000/um2pixel, ...
            bregma_offset_y-point_ap*1000/um2pixel,'oy','markersize',20);       
        warning('um/px: needs more accurate measurement');

    case 'point_aligned'
        % Plot a point at specific coordinates (in CCF-aligned transform)
        point_ml_ap = fliplr(color);

        point_ccf = bregma([3,1]) + point_ml_ap*100.*[1,-1];
        point_ccf_aligned = [point_ccf,1]*ccf_tform.T;
        
        plot(point_ccf_aligned(1),point_ccf_aligned(2), ...
            '.y','markersize',20);
        
    case 'ccf'
        % Plot CCF borders aligned to master retinotopy        
        h = cellfun(@(areas) cellfun(@(outline) ...
            plot(outline(:,1),outline(:,2),'color',color),areas,'uni',false), ...
                dorsal_cortex_borders_aligned,'uni',false);

    case 'retinotopy'
        % Show the master retinotopy
        master_vfs_fn = fullfile(alignment_path,'master_vfs.mat');
        load(master_vfs_fn);
        imagesc(master_vfs);

    otherwise
        warning('Invalid reference: %s',type);
        
end

% %% Create top-down cortical boundaries from CCF (tilt CCF) (RUN ONCE)
% %
% % NOTE: CCF IN NATIVE FORMAT (NOT SCALED OR ROTATED)
% 
% % Set save path
% alignment_path = fullfile(plab.locations.server_path, ...
%     'Users','Andy_Peters','widefield_alignment');
% 
% % Load in the annotated Allen volume and names
% allen_atlas_path = fileparts(which('template_volume_10um.npy'));
% av = readNPY(fullfile(allen_atlas_path, 'annotation_volume_10um_by_index.npy'));
% st = loadStructureTree(fullfile(allen_atlas_path, 'structure_tree_safe_2017.csv'));
% 
% % %%%%%% WORKING HERE
% % % Transform CCF as the NTE/pinpoint does
% % % (translation values from our bregma estimate: AP/ML from Paxinos, DV from
% % % rough MRI estimate)
% % bregma_ccf = [570.5,520,44]; % [ML,AP,DV]
% % ccf_translation_tform = eye(4)+[zeros(3,4);-bregma_ccf,0];
% % 
% % % (scaling "Toronto MRI transform", reflect AP/ML, convert 10um to 1mm)
% % scale = [0.952,-1.031,0.885]./100; % [ML,AP,DV]
% % ccf_scale_tform = eye(4).*[scale,1]';
% % 
% % % (rotation values from IBL estimate)
% % ap_rotation = 5; % tilt the CCF 5 degrees nose-up
% % ccf_rotation_tform = ...
% %     [1 0 0 0; ...
% %     0 cosd(ap_rotation) -sind(ap_rotation) 0; ...
% %     0 sind(ap_rotation) cosd(ap_rotation) 0; ...
% %     0 0 0 1];
% % 
% % ccf_bregma_tform_matrix = ccf_translation_tform*ccf_scale_tform*ccf_rotation_tform;
% % ccf_bregma_tform = affine3d(ccf_bregma_tform_matrix);
% % %%%%%%%%%%%%%
% 
% % Get first brain pixel from top-down, get annotation at that point
% [~,top_down_depth] = max(av>1, [], 2);
% top_down_depth = squeeze(top_down_depth);
% 
% [xx,yy] = meshgrid(1:size(top_down_depth,2), 1:size(top_down_depth,1));
% dorsal_ccf_annotation = reshape(av(sub2ind(size(av),yy(:), ...
%     top_down_depth(:),xx(:))), size(av,1), size(av,3));
% 
% % Get all labelled areas
% used_areas = unique(dorsal_ccf_annotation(:));
% 
% % Restrict to only cortical areas
% structure_id_path = cellfun(@(x) textscan(x(2:end),'%d', 'delimiter',{'/'}),st.structure_id_path);
% 
% ctx_path = [997,8,567,688,695,315];
% ctx_idx = find(cellfun(@(id) length(id) > length(ctx_path) & ...
%     all(id(min(length(id),length(ctx_path))) == ctx_path(min(length(id),length(ctx_path)))),structure_id_path));
% 
% plot_areas = intersect(used_areas,ctx_idx);
% 
% % Get outlines of all areas (full/left hemi)
% dorsal_cortex_borders = cell(size(plot_areas));
% for curr_area_idx = 1:length(plot_areas)
%     dorsal_cortex_borders{curr_area_idx} = bwboundaries(dorsal_ccf_annotation == plot_areas(curr_area_idx));
% end
% 
% % Plot the borders
% ccf_color_hex = st.color_hex_triplet;
% ccf_color_hex(cellfun(@numel,ccf_color_hex)==5) = {'019399'}; % special case where leading zero was evidently dropped
% ccf_cmap_c1 = cellfun(@(x)hex2dec(x(1:2)), ccf_color_hex, 'uni', false);
% ccf_cmap_c2 = cellfun(@(x)hex2dec(x(3:4)), ccf_color_hex, 'uni', false);
% ccf_cmap_c3 = cellfun(@(x)hex2dec(x(5:6)), ccf_color_hex, 'uni', false);
% ccf_cmap = horzcat(vertcat(ccf_cmap_c1{:}),vertcat(ccf_cmap_c2{:}),vertcat(ccf_cmap_c3{:}))./255;
% 
% figure; hold on;axis image; set(gca,'ydir','reverse');
% imagesc(dorsal_ccf_annotation);
% cellfun(@(x) cellfun(@(x) plot(x(:,2),x(:,1),'k','linewidth',2),x,'uni',false), ...
%     dorsal_cortex_borders,'uni',false);
% colormap(ccf_cmap);
% title('Dorsal cortical areas');
% 
% % Save borders
% save_fn = fullfile(alignment_path,'dorsal_cortex_borders');
% save(save_fn,'dorsal_cortex_borders','dorsal_ccf_annotation');
% disp(['Saved ' save_fn]);












