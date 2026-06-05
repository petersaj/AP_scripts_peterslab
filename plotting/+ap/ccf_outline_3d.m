function ccf_outline = ccf_outline_3d(ccf_3d_axes,av,plot_structures)
% Draw 3D outline of CCF
%
% ccf_3d_axes - axes to plot on (new figure if excluded) 
% av - CCF annotated volume (loaded if excluded) 
% plot_structures - string/cell array of names/acronyms, "brain" or empty
% is whole-brain outline, e.g. ["brain","caudoputamen"])

arguments 
    ccf_3d_axes = []
    av = []
    plot_structures = "brain";
end

% Ensure plot_structures is string array
plot_structures = string(plot_structures);

% Load CCF if not provided
if isempty(av)
    [av,~,st] = ap_histology.load_ccf;
end

% Create up axes if one not provided
if isempty(ccf_3d_axes)
    fig = figure('color','w');
    ccf_3d_axes = axes(fig);  
    hold(ccf_3d_axes,'on');
end

% Set up 3D axes
set(ccf_3d_axes,'ZDir','reverse');
hold(ccf_3d_axes,'on');
axis(ccf_3d_axes,'vis3d','equal','off','manual');
view([-30,25]);
axis tight;
h_rotation = rotate3d(ccf_3d_axes);
h_rotation.Enable = 'on';

% Draw CCF outline
slice_spacing = 5;
for curr_structure = reshape(plot_structures,1,[])
    if strcmpi(curr_structure,"brain")
        % Plot whole brain if `area` is empty
        plot_structure_color = [0.7,0.7,0.7];
        plot_ccf_volume = bwmorph3(bwmorph3(av(1:slice_spacing:end, ...
            1:slice_spacing:end,1:slice_spacing:end)>1,'majority'),'majority');
        face_alpha = 0.1;
    else
        % Plot `area` if specified
        plot_structure_idx = find(strcmpi(st.safe_name,curr_structure) | ...
            strcmpi(st.acronym,curr_structure));
        plot_structure_id = st.structure_id_path{plot_structure_idx};
        plot_ccf_idx = find(cellfun(@(x) contains(x,plot_structure_id), ...
            st.structure_id_path));
        plot_structure_color = hex2dec(reshape(st.color_hex_triplet{plot_structure_idx},2,[])')./255;
        plot_ccf_volume = ismember(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end),plot_ccf_idx);
        face_alpha = 0.2;
    end

    plot_ccf_patchdata = isosurface(permute(plot_ccf_volume,[3,1,2]),0.5);
    ccf_outline = patch(ccf_3d_axes, ...
        'Vertices',plot_ccf_patchdata.vertices*slice_spacing, ...
        'Faces',plot_ccf_patchdata.faces, ...
        'FaceColor',plot_structure_color,'EdgeColor','none','FaceAlpha',face_alpha);
end


