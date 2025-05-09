function cmap = colormap(cmap_type,n_colors,cmap_gamma)
%
% cmap = ap.colormap(cmap_type,n_colors,cmap_gamma)
%
% cmap_type: color extremes as letters, or 'tube' for categorical (max n=15)
% (currently defined: W,K,G,B,P,R, base must be W or K)
% e.g. 'WR' = white->red, 'BKR' = blue->black->red
%
% n_colors: number of colors along the map (default = 2^8 = 256)
%
% cmap_gamma: value to gamma correct (resamples colormap values according to
% x^gamma, which emphasizes higher/lower values to give better contrast in
% different places) (default = 1: no gamma change)
%
% Interpolate in CIELAB colorspace to ensure linear luminance change, see:
% https://blogs.mathworks.com/steve/2006/05/09/a-lab-based-uniform-color-scale/

% Force cmap_type upper-case
cmap_type = upper(cmap_type);

%% If cmap_type is 'tube', use categorical colors
% (line colors from the London tube)

if strcmpi(cmap_type,'tube')

    tube_colors = [ ...
        0.7020    0.3882    0.0196
        0         0         0
        0.8902    0.1255    0.0902
        0    0.2118    0.5333
        1.0000    0.8275         0
        0    0.5961    0.8314
        0    0.4706    0.1647
        0.5843    0.8039    0.7294
        0.4118    0.3137    0.6314
        0    0.6431    0.6549
        0.9529    0.6627    0.7333
        0.9333    0.4863    0.0549
        0.6275    0.6471    0.6627
        0.5176    0.7216    0.0902
        0.6078         0    0.3373
        ];

    if exist('n_colors','var') && n_colors > size(tube_colors,1)
        error('%d colors requested, only %d colors in ''tube'' cmap',n_colors,size(tube_colors,1));
    elseif ~exist('n_colors','var') || isempty(n_colors)
        n_colors = size(tube_colors,1);
    end

    if exist('cmap_gamma','var') && ~isempty(cmap_gamma)
        error('Gamma specified, but not applicable for ''tube'' cmap');
    end

    cmap = tube_colors(1:n_colors,:);

    return
end

%% If color gradient chosen, create gradient

% Set default number of colors
if ~exist('n_colors','var') || isempty(n_colors)
    n_colors = 2^8; % default 8-bit color
end

% If diverging colormap, force number of colors to be odd
if length(cmap_type) == 3
    n_colors = n_colors - (1-mod(n_colors,2));
end

% Set default gamma
if ~exist('cmap_gamma','var') || isempty(cmap_gamma)
    cmap_gamma = 1; % default: no gamma change
end

% If gamma, supersample colors to resample finer gradients
% (note: randomly 100-fold now, should base on gamma...)
% (force odd number of colors if diverging colormap)
if cmap_gamma == 1
    n_colors_generate = n_colors;
elseif    cmap_gamma ~= 1
    n_colors_generate = n_colors*100 + (length(cmap_type) == 3);
end

% Set maximum color chroma
max_chroma = 80;

% Set black/white luminance
K_lum = 1;
W_lum = 99;

% Set colors in CIELAB space (angle, min chroma, max chroma)
col = struct;
col.W = 0;
col.K = 0;
col.G = 139;
col.P = 315;
col.B = 250;
col.R = 30;

% (to plot colors while writing/troubleshooting)
plot_col = false;
if plot_col
    all_col = cell2mat(cellfun(@(x) permute(lab2rgb(interp1([K_lum,W_lum], ...
        [W_lum,max_chroma*cosd(x),max_chroma*sind(x); ...
        W_lum,max_chroma*cosd(x),max_chroma*sind(x)], ...
        linspace(K_lum,W_lum,256))),[1,3,2]),struct2cell(col),'uni',false)');
    figure;imagesc(all_col);
    set(gca,'XTickLabels',fieldnames(col),'YTick','');
end

% Check base color is black or white (others not supported now)
cmap_type_iswk = ismember(cmap_type',{'W','K'})';
if ~cmap_type_iswk(ceil(length(cmap_type)/2))
    error('Base color must be white or black')
end

% Set luminance/chroma interpolation direction from base color
lum_leeway = 20; % keep extremes of colormap away from the max
switch cmap_type(ceil(length(cmap_type)/2))
    case 'W'
        lum_interp = [W_lum,K_lum+lum_leeway];
    case 'K'
        lum_interp = [K_lum,W_lum-lum_leeway];
end

% Get colormap X values to interpolate across (X^gamma)
cmap_x_defined = [0,1];
cmap_x_interp = linspace(cmap_x_defined(1),cmap_x_defined(end), ...
    ceil(n_colors_generate/(length(cmap_type)-1))).^cmap_gamma;

% Define gradients by interpolating in CIELAB space and converting to RGB
if length(cmap_type) == 2
    % 2 colors = gradient

    cmap = lab2rgb(interp1(cmap_x_defined, ...
        [lum_interp(1),0,0;
        lum_interp(2), ...
        ~cmap_type_iswk(2)*max_chroma*cosd(col.(cmap_type(2))(1)), ...
        ~cmap_type_iswk(2)*max_chroma*sind(col.(cmap_type(2))(1))], ...
        cmap_x_interp,'spline'));

elseif length(cmap_type) == 3
    % 3 colors = diverging 

    cmap_bot = lab2rgb(interp1(cmap_x_defined, ...
        [lum_interp(1),0,0;
        lum_interp(2), ...
        ~cmap_type_iswk(1)*max_chroma*cosd(col.(cmap_type(1))(1)), ...
        ~cmap_type_iswk(1)*max_chroma*sind(col.(cmap_type(1))(1))], ...
        cmap_x_interp,'spline'));

    cmap_top = lab2rgb(interp1(cmap_x_defined, ...
        [lum_interp(1),0,0;
        lum_interp(2), ...
        ~cmap_type_iswk(3)*max_chroma*cosd(col.(cmap_type(3))(1)), ...
        ~cmap_type_iswk(3)*max_chroma*sind(col.(cmap_type(3))(1))], ...
        cmap_x_interp,'spline'));
    
    % Combine bottom/top (center color is replicated - remove one)
    cmap = vertcat(flipud(cmap_bot),cmap_top(2:end,:));

end

% Force colormap into RGB gamut (compromises linearity?)
cmap(cmap > 1) = 1;
cmap(cmap < 0) = 0;






