function unit_dots = plot_unit_depthrate(spike_times,spike_templates,template_depths,probe_areas,plot_axes)
% plot_unit_depthrate(spike_templates,template_depths,probe_areas,plot_axes)
%
% Plot unit depth vs spike rate, with areas (if available)

if ~exist('plot_axes','var') || isempty(plot_axes)
    figure('Units','normalized','Position',[0.02,0.2,0.1,0.6])
    plot_axes = axes;
end

hold(plot_axes,'on');

% Plot units (depth vs normalized rate) with areas
if exist('probe_areas','var') && ~isempty(probe_areas)
    probe_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
        probe_areas{1}.color_hex_triplet,'uni',false)),[1,3,2]);

    probe_areas_boundaries = probe_areas{1}.probe_depth;
    probe_areas_centers = mean(probe_areas_boundaries,2);

    probe_areas_image_depth = 0:1:max(probe_areas_boundaries,[],'all');
    probe_areas_image_idx = interp1(probe_areas_boundaries(:,1), ...
        1:height(probe_areas{1}),probe_areas_image_depth, ...
        'previous','extrap');
    probe_areas_image = ones(length(probe_areas_image_idx),1,3);
    probe_areas_image(~isnan(probe_areas_image_idx),:,:) = ...
        probe_areas_rgb(probe_areas_image_idx(~isnan(probe_areas_image_idx)),:,:);
    
    yyaxis left; set(gca,'YDir','reverse');
    image(plot_axes,[0.5],probe_areas_image_depth,probe_areas_image);
    yline(unique(probe_areas_boundaries(:)),'color','k','linewidth',1);
    set(plot_axes,'YTick',probe_areas_centers,'YTickLabels',probe_areas{1}.acronym);
end

yyaxis right;
set(gca,'YDir','reverse');
spike_rate = (accumarray(findgroups(spike_templates),1)+1)/ ...
    diff(prctile(spike_times,[0,100]));

unit_dots = scatter( ...
    spike_rate,template_depths(unique(spike_templates)),20,'k','filled');
ylabel('Depth (\mum)')
xlabel('Spike rate')
set(gca,'XScale','log');

[plot_axes.YAxis.Color] = deal('k');

yyaxis left; 
ylim([0, 3840]);

yyaxis right;
ylim([0, 3840]);


