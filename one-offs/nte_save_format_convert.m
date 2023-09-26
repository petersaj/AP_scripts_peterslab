% Convert old NTE to new NTE save format 

%% Load atlas

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

%% Find all position files

file_pattern = plab.locations.make_server_filename('*','*',[],'ephys','*probe_positions.mat');
probe_position_filenames = dir(file_pattern);

%% Convert files in loop

for curr_file = 1:length(probe_position_filenames)

    curr_filename = fullfile(probe_position_filenames(curr_file).folder, ...
        probe_position_filenames(curr_file).name);
    
    load(curr_filename);

    % Get positions along each probe
    n_probes = length(probe_positions_ccf);

    probe_areas = cell(n_probes,1);
    for curr_probe = 1:n_probes

        probe_coords_depth = [1:3840]-1;

        probe_coords_ccf = round(interp1([1,3840]-1, ...
            probe_positions_ccf{curr_probe}',probe_coords_depth));

        probe_coords_ccf_inbounds = all(probe_coords_ccf > 0 & ...
            probe_coords_ccf <= size(av),2);

        probe_location_idx = ...
            sub2ind(size(av), ...
            probe_coords_ccf(probe_coords_ccf_inbounds,1), ...
            probe_coords_ccf(probe_coords_ccf_inbounds,2), ... 
            probe_coords_ccf(probe_coords_ccf_inbounds,3));

        % Get boundaries of areas and area IDs
        probe_n_coords = length(probe_coords_depth);
        probe_area_idx_sampled = ones(probe_n_coords,1);
        probe_area_idx_sampled(probe_coords_ccf_inbounds) = av(probe_location_idx);
        probe_area_bins = [1;(find(diff(probe_area_idx_sampled)~= 0)+1);probe_n_coords];
        probe_area_boundaries = [probe_area_bins(1:end-1),probe_area_bins(2:end)];

        probe_area_idx = probe_area_idx_sampled(probe_area_boundaries(:,1));

        % Store structure tree entries for probe areas (ammend start depth for each area)
        store_areas_idx = probe_area_idx > 1; % only use areas in brain (idx > 1)
        curr_probe_areas = st(probe_area_idx(store_areas_idx),:);
        curr_probe_areas.probe_depth = probe_coords_depth(probe_area_boundaries(store_areas_idx,:));

        probe_areas{curr_probe} = curr_probe_areas;

    end

    % Save converted data
    save(curr_filename,'probe_positions_ccf','probe_areas');
    fprintf('Converted: %s\n',curr_filename)

end




