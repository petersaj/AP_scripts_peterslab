function unit_dots = plot_unit_depthrate(plot_axes,split_shanks)
% plot_unit_depthrate(split_shanks)
%
% % Plot unit depth vs spike rate, with areas (if available)
% (grabs variables from base workspace)
% 
% INPUTS:
% plot_axes - axis handle to plot shanks
% split_shanks - shanks on separate axes (true) or one axis (false). True
% by default without area_axes, forced as false with area_axes.
% plot_probe - which probe to plot (if multiple)

arguments
    plot_axes = []
    split_shanks = true
end

% Pull variables from base workspace
spike_times_openephys = evalin('base','spike_times_openephys');
spike_templates = evalin('base','spike_templates');
template_tipdist = evalin('base','template_tipdist');
template_shanks = evalin('base','template_shanks');
load_probe = evalin('base','load_probe');
try
    probe_areas = evalin('base','probe_areas');
catch me
end

% Set up axes
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

% Draw areas for all shanks
if exist('probe_areas','var')
    for curr_shank = 1:n_shanks

        hold(shank_axes(curr_shank),'on');

        % LEGACY SUPPORT
        % - If `probe_depth` instead of `tip_distance`, calculate
        if ~any(strcmp(probe_areas{1}.Properties.VariableNames,'tip_distance')) && ...
                any(strcmp(probe_areas{1}.Properties.VariableNames,'probe_depth'))
            probe_areas{1}.tip_distance = 3840 - probe_areas{1}.probe_depth;
        end
        % - if no `probe_shank`, add 1's
        if ~any(strcmp(probe_areas{1}.Properties.VariableNames,'probe_shank'))
            probe_areas{1}.probe_shank = ones(height(probe_areas{1}),1);
        end
        
        % Get areas on current shank
        curr_shank_areas = find(probe_areas{1}.probe_shank==curr_shank);

        % Plot areas as rectangles
        for curr_area = reshape(curr_shank_areas,1,[])
            curr_rgb = hex2dec(mat2cell(probe_areas{1}.color_hex_triplet{curr_area},1,repmat(2,1,3)))./255;
            curr_y = probe_areas{1}.tip_distance(curr_area,:);
            rectangle(shank_axes(curr_shank),'Position',[shank_xoffset(curr_shank),min(curr_y), ...
                1,abs(diff(curr_y))], ...
                'FaceColor',curr_rgb,'EdgeColor','none');
        end

        % Label area centers
        text(shank_axes(curr_shank), ...
            repelem(shank_xoffset(curr_shank),length(curr_shank_areas),1), ...
            probe_areas{1}.tip_distance(curr_shank_areas,1), ...
            probe_areas{1}.acronym(curr_shank_areas));

        set(shank_axes(curr_shank),'YTick',0:0.5:max(probe_areas{1}.tip_distance,[],'all'));
    end
end


% Plot units
norm_spike_count = normalize(log10(accumarray(findgroups(spike_templates),1)),'range');

unit_xplot = norm_spike_count + reshape(shank_xoffset(template_shanks),[],1);
unit_yplot = template_tipdist/1000;

if split_shanks
    unit_dots = arrayfun(@(shank) scatter(shank_axes(shank), ...
        unit_xplot(template_shanks==shank), ...
        unit_yplot(template_shanks==shank),20,'k','filled', ...
        'UserData',struct('shank',template_shanks(template_shanks==shank))), ...
        1:n_shanks);
else
    unit_dots = scatter(unique(shank_axes), ...
        unit_xplot,unit_yplot,20, ...
        'k','filled', ...
        'UserData',struct('shank',template_shanks));
end
xlabel(shank_axes,'Rate')

% Set y-limit around units
ylim(prctile(unit_yplot,[0,100])+[-0.1,0.1]);
