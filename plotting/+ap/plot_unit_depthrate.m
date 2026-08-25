function plot_handles = plot_unit_depthrate(plot_axes,split_shanks,plot_shank)
% plot_handles = plot_unit_depthrate(plot_axes,split_shanks,plot_shank)
%
% % Plot unit depth vs spike rate, with areas (if available)
% (grabs variables from base workspace)
% 
% INPUTS:
% plot_axes - axis handle to plot shanks
% split_shanks - true/false split shanks on separate axes (default true without area_axes, forced false with area_axes)
% plot_shank - shank to plot (all by default, single if specified)
%
% OUTPUTS
% plot_handles - [.unit_dots, .area_rectangles] graphics handles (splits if
% shanks across axes, combines if single axis)


arguments
    plot_axes = []
    split_shanks = true
    plot_shank = [];
end

% Pull variables from base/caller workspace
try
    use_workspace = 'base';
    spike_times_openephys = evalin(use_workspace,'spike_times_openephys');
catch me
    try
        use_workspace = 'caller';
        spike_times_openephys = evalin(use_workspace,'spike_times_openephys');
    catch me
        error('Spike variables not found in workspace')
    end
end

spike_templates = evalin(use_workspace,'spike_templates');
template_tipdist = evalin(use_workspace,'template_tipdist');
template_shanks = evalin(use_workspace,'template_shanks');
load_probe = evalin(use_workspace,'load_probe');
channel_positions = evalin(use_workspace,'channel_positions');
try probe_areas = evalin(use_workspace,'probe_areas'); end
try template_qc_labels = evalin(use_workspace,'template_qc_labels'); end

% Flag somata vs axons (if included)
if exist('template_qc_labels','var')
    unit_soma = ~strcmp(template_qc_labels,'axon');
else
    unit_soma = true(size(template_tipdist));
end

% Set up axes
plot_handles = struct;

n_shanks = max(template_shanks);
shank_spacing = 1.5;
if isempty(plot_axes)
    % No axes specified: make figure and axes
    figure('Units','normalized','Position',[0.02,0.2,0.1,0.6]);
    if split_shanks
        h = tiledlayout(1,n_shanks,'TileSpacing','none');
        shank_axes = arrayfun(@(x) nexttile(h),1:n_shanks);
        shank_xoffset = zeros(1,n_shanks);
    else
        shank_axes(1:n_shanks) = deal(axes);
        shank_xoffset = (1:n_shanks)*shank_spacing;
    end
else
    % Plot axes specified: use axes, don't split shanks
    split_shanks = false; % override split shanks if axes specified
    shank_axes(1:n_shanks) = deal(plot_axes);
    shank_xoffset = (1:n_shanks)*shank_spacing;
end
linkaxes(shank_axes,'y')

% Draw areas by shank
if isempty(plot_shank)
    plot_shank = 1:n_shanks;
end
area_rectangles = cell(n_shanks,1);
if exist('probe_areas','var')
    for curr_shank = reshape(plot_shank,1,[])

        hold(shank_axes(curr_shank),'on');
        
        % Get areas on current shank
        curr_shank_areas = find(probe_areas.probe_shank==curr_shank);

        % Plot areas as rectangles
        for curr_area_idx = 1:length(curr_shank_areas)
            curr_area = curr_shank_areas(curr_area_idx);
            curr_color = ['#',probe_areas.color_hex_triplet{curr_area}];
            curr_y = probe_areas.tip_distance(curr_area,:);
            area_rectangles{curr_shank}(curr_area_idx,1) = ...
                rectangle(shank_axes(curr_shank),'Position',[shank_xoffset(curr_shank),min(curr_y), ...
                1,abs(diff(curr_y))], ...
                'FaceColor',curr_color,'EdgeColor','none');
        end

        % Label area centers
        text(shank_axes(curr_shank), ...
            repelem(shank_xoffset(curr_shank),length(curr_shank_areas),1), ...
            probe_areas.tip_distance(curr_shank_areas,1), ...
            probe_areas.acronym(curr_shank_areas));

        set(shank_axes(curr_shank),'YTick',0:0.5:max(probe_areas.tip_distance,[],'all'));
    end
end

if split_shanks
    [plot_handles(1:n_shanks).area_rectangles] = deal(area_rectangles{:});
else
    plot_handles.area_rectangles = vertcat(area_rectangles{:});
end

% Plot recorded sites
% (change x-values to be normalized by shank)
shank_spacing = 250;
shank_borders = (0:4)*shank_spacing-shank_spacing/2;
channel_shanks = discretize(channel_positions(:,1),shank_borders);

for curr_shank = reshape(plot_shank,1,[])
    plot(shank_axes(curr_shank),shank_xoffset(curr_shank)+0.1, ...
        channel_positions(channel_shanks == curr_shank,2)/1000, ...
        'squarek','MarkerSize',5);
end

% Plot units (soma)
norm_spike_count = normalize(log10(accumarray(findgroups(spike_templates),1)),'range');

unit_xplot = norm_spike_count + reshape(shank_xoffset(template_shanks),[],1);
unit_yplot = template_tipdist/1000;

soma_axon_color = [0,0,0;0.5,0.5,0.5];

if split_shanks
    unit_dots = arrayfun(@(shank) scatter(shank_axes(shank), ...
        unit_xplot(template_shanks==shank), ...
        unit_yplot(template_shanks==shank),20, ... 
        soma_axon_color(2-unit_soma(template_shanks==shank,:),:),'filled', ...
        'UserData',struct('shank',template_shanks(template_shanks==shank))), ...
        plot_shank,'uni',false);
    [plot_handles.unit_dots] = deal(unit_dots{:});
else
    plot_units = ismember(template_shanks,plot_shank);
    plot_handles.unit_dots = scatter(shank_axes(1), ...
        unit_xplot(plot_units), ...
        unit_yplot(plot_units),20, ...
        soma_axon_color(2-unit_soma(plot_units,:),:),'filled', ...
        'UserData',struct('shank',template_shanks(plot_units)));    
end

xlabel(shank_axes,'Rate')

% Set y-limit around units
ylim(shank_axes,prctile(unit_yplot,[0,100])+[-0.1,0.1]);
