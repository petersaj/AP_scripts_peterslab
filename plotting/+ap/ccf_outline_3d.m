function ccf_outline = ccf_outline_3d(ccf_3d_axes,av)
% Draw 3D outline of CCF
%
% ccf_3d_axes - axes to plot on (new figure if excluded)
% av - CCF annotated volume (loaded if excluded)

arguments 
    ccf_3d_axes = []
    av = []
end

% Load CCF if not provided
if isempty(av)
    av = ap_histology.load_ccf;
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
ccf_volume = ...
    bwmorph3(bwmorph3(av(1:slice_spacing:end, ...
    1:slice_spacing:end,1:slice_spacing:end)>1,'majority'),'majority');
ccf_outline_patchdata = isosurface(permute(ccf_volume,[3,1,2]),0.5);
ccf_outline = patch(ccf_3d_axes, ...
    'Vertices',ccf_outline_patchdata.vertices*slice_spacing, ...
    'Faces',ccf_outline_patchdata.faces, ...
    'FaceColor',[0.7,0.7,0.7],'EdgeColor','none','FaceAlpha',0.1);
