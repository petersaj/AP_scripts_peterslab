function probe_insertion_ccf = get_probe_insertion_ccf(animal)
% probe_insertion_ccf = get_probe_insertion_ccf(animal)
%
% Get probe insertion points in CCF coords from AP_histology annotations

% Load histology processing
histology_filepattern = plab.locations.filename('server',animal,[],[],'histology','**','AP_histology_processing.mat');
histology_fn = dir(histology_filepattern);
load(fullfile(histology_fn.folder,histology_fn.name));

% Get labels and vertices of all probes
probe_labels = string({AP_histology_processing.annotation.label});
probe_ccf_vertices = arrayfun(@(x) horzcat( ...
    vertcat(x.vertices_ccf.ap), ...
    vertcat(x.vertices_ccf.dv), ...
    vertcat(x.vertices_ccf.ml)), ...
    AP_histology_processing.annotation,'uni',false);

% Loop through probes, find insertion point
probe_insertion_ccf = cell(size(probe_ccf_vertices));

for curr_probe = 1:length(probe_ccf_vertices)

    % Load CCF atlas
    [av,~,~] = ap_histology.load_ccf;

    % Get line of best fit through mean of marked points
    [~,~,probe_fit] = svd(probe_ccf_vertices{curr_probe} - mean(probe_ccf_vertices{curr_probe},1),0);

    probe_direction = probe_fit(:,1);
    % (ensure vector goes downward in DV)
    probe_direction(2) = abs(probe_direction(2));
    probe_vector = mean(probe_ccf_vertices{curr_probe},1) + padarray(probe_direction',[1,0],0,'pre');

    % Grab AV values across probe trajectory
    max_eval = round(sqrt(max(size(av)).^2*2));
    eval_points_probe = (-max_eval:max_eval);
    eval_points_ccf = round(interp1([0,norm(diff(probe_vector))],probe_vector, ...
        eval_points_probe,'linear','extrap'));

    eval_points_ccf_valid = eval_points_ccf(all(eval_points_ccf > 0 & eval_points_ccf <= size(av),2),:);

    eval_points_ccf_idx = sub2ind(size(av),eval_points_ccf_valid(:,1), ...
        eval_points_ccf_valid(:,2),eval_points_ccf_valid(:,3));

    probe_trajectory_av = av(eval_points_ccf_idx);

    % Get CCF coordinate of insertion point as first point with brain label
    probe_insertion_ccf{curr_probe} = eval_points_ccf_valid(find(probe_trajectory_av>1,1),:);

end




