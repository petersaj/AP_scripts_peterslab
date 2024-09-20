classdef ccf_draw < handle
    % ccf_draw
    %
    % Draw CCF atlas in 2D/3D views
    % 
    % Draw areas
    % Draw probes

    properties
        ccf_fig % figure
        ccf_axes % axes (1-3 = 2D, 4 = 3D)
        av % Annotated volume
        st % Structure tree
    end

    % User-called methods
    methods

        %% Constructor 
        function obj = ccf_draw

            % Load atlas
            allen_atlas_path = fileparts(which('template_volume_10um.npy'));

            obj.av = readNPY(fullfile(allen_atlas_path,'annotation_volume_10um_by_index.npy'));
            obj.st = loadStructureTree(fullfile(allen_atlas_path,'structure_tree_safe_2017.csv'));

            % Create figure and axes
            obj.ccf_fig = figure('color','w');
            h = tiledlayout('flow','Tilespacing','none');
            obj.ccf_axes = gobjects(4,1);
            for curr_ax = 1:4
                obj.ccf_axes(curr_ax) = nexttile(h);
                hold(obj.ccf_axes(curr_ax),'on');
            end

            % Properties for 2D axes
            set(obj.ccf_axes(1:3),'YDir','reverse')
            axis(obj.ccf_axes(1:3),'equal','off')

            % Properties for 3D axis
            set(obj.ccf_axes(4),'ZDir','reverse');
            axis(obj.ccf_axes(4),'vis3d','equal','off','manual','tight');
            view(obj.ccf_axes(4),[-30,25]);
            h_rot = rotate3d(obj.ccf_axes(4));
            h_rot.Enable = 'on';

            % Plot 2D brain outlines
            for curr_view = 1:3
                curr_outline = bwboundaries(squeeze((max(obj.av,[],curr_view)) > 1));
                % (only plot largest outline)
                [~,curr_outline_idx] = max(cellfun(@length,curr_outline));
                curr_outline_reduced = reducepoly(curr_outline{curr_outline_idx});
                plot(obj.ccf_axes(curr_view), ...
                    curr_outline_reduced(:,2), ...
                    curr_outline_reduced(:,1),'k','linewidth',2);
                % (draw 1mm scalebar)
                %     line(obj.ccf_axes(curr_view),[0,0],[0,100],'color','k','linewidth',2);
            end
            linkaxes(obj.ccf_axes(1:3));

            % Plot 3D brain outlines
            slice_spacing = 5;
            brain_volume = ...
                bwmorph3(bwmorph3(obj.av(1:slice_spacing:end, ...
                1:slice_spacing:end,1:slice_spacing:end)>1,'majority'),'majority');
            brain_outline_patchdata = isosurface(permute(brain_volume,[3,1,2]),0.5);
            brain_outline = patch( ...
                obj.ccf_axes(4), ...
                'Vertices',brain_outline_patchdata.vertices*slice_spacing, ...
                'Faces',brain_outline_patchdata.faces, ...
                'FaceColor',[0.7,0.7,0.7],'EdgeColor','none','FaceAlpha',0.1);
            axis tight

        end

        %% Plot areas
        % Plot area by name
        function obj = draw_name(obj,structure_name)

            % Find structure (exact match)
            plot_structure = find(strcmpi(obj.st.safe_name,structure_name));
            obj.draw_structure(plot_structure);
        end

        % Plot area by search
        function obj = draw_search(obj)
            % Prompt for which structures to show (only structures which are
            % labelled in the slice-spacing downsampled annotated volume)
            structure_search = lower(inputdlg('Search structures'));
            structure_match = find(contains(lower(obj.st.safe_name),structure_search));

            selected_structure = listdlg('PromptString','Select a structure to plot:', ...
                'ListString',obj.st.safe_name(structure_match),'ListSize',[520,500], ...
                'SelectionMode','single');

            plot_structure = structure_match(selected_structure);
            obj.draw_structure(plot_structure);
        end

        %% Plot probes

        function obj = draw_probes_nte(obj,animal,probe_color)
            % draw_probes_nte(obj,animal,probe_color)
            % Trajectory explorer probe position files
            nte_filepattern = plab.locations.filename('server',animal,'*',[],'ephys','*probe_positions.mat');
            nte_fns = dir(nte_filepattern);

            for curr_recording = 1:length(nte_fns)
                load(fullfile(nte_fns(curr_recording).folder, nte_fns(curr_recording).name));
                % Loop through probes and draw
                for curr_probe = 1:length(probe_positions_ccf)
                    % Probe line is directly stored
                    probe_line = probe_positions_ccf{curr_probe}';
                    obj.draw_probes(probe_line,probe_color);
                end
            end
        end

        function obj = draw_probes_histology(obj,animal,probe_color)
            % Histology probe position files
            histology_filepattern = plab.locations.filename('server',animal,[],[],'histology','**','probe_ccf.mat');
            histology_fn = dir(histology_filepattern);

            load(fullfile(histology_fn.folder,histology_fn.name));

            % Loop through probes and draw
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
                obj.draw_probes(probe_line,probe_color);
            end
        end

    end

    %% Internal draw methods
    methods (Access = protected)
        function draw_structure(obj,plot_structure)
            % Get all areas within and below the selected hierarchy level
            plot_structure_id = obj.st.structure_id_path{plot_structure};
            plot_ccf_idx = find(cellfun(@(x) contains(x,plot_structure_id), ...
                obj.st.structure_id_path));

            % Get structure color and volume
            slice_spacing = 5;
            plot_structure_color = hex2dec(reshape(obj.st.color_hex_triplet{plot_structure},2,[])')./255;
            plot_ccf_volume = ismember(obj.av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end),plot_ccf_idx);

            % Plot 2D structure
            for curr_view = 1:3
                curr_outline = bwboundaries(squeeze((max(plot_ccf_volume,[],curr_view))));
                cellfun(@(x) plot(obj.ccf_axes(curr_view),x(:,2)*slice_spacing, ...
                    x(:,1)*slice_spacing,'color',plot_structure_color,'linewidth',2),curr_outline)
            end

            % Plot 3D structure
            structure_3d = isosurface(permute(plot_ccf_volume,[3,1,2]),0);
            structure_alpha = 0.2;
            patch(obj.ccf_axes(4), ...
                'Vertices',structure_3d.vertices*slice_spacing, ...
                'Faces',structure_3d.faces, ...
                'FaceColor',plot_structure_color,'EdgeColor','none','FaceAlpha',structure_alpha);
        end

        function draw_probes(obj,probe_line,probe_color)
            % Draw probe in 3D view
            line(obj.ccf_axes(4),probe_line(:,1),probe_line(:,3),probe_line(:,2), ...
                'linewidth',2,'color',probe_color)

            % Draw probes on coronal + saggital
            line(obj.ccf_axes(1),probe_line(:,3),probe_line(:,2),'linewidth',2,'color',probe_color);
            line(obj.ccf_axes(2),probe_line(:,3),probe_line(:,1),'linewidth',2,'color',probe_color);
            line(obj.ccf_axes(3),probe_line(:,2),probe_line(:,1),'linewidth',2,'color',probe_color);

            % Draw probe start/end on horizontal
            plot(obj.ccf_axes(2), probe_line(1,3),probe_line(1,1), ...
                'o','MarkerSize',5,'color',probe_color);
            plot(obj.ccf_axes(2), probe_line(end,3),probe_line(end,1), ...
                '.','MarkerSize',20,'color',probe_color);
        end

    end

end



