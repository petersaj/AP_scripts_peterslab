function unit_dots = plot_unit_depthrate(spike_templates,template_depths,probe_areas)
% plot_unit_depthrate(spike_templates,template_depths,probe_areas)
%
% Plot unit depth vs normalized log spike rate, with areas (if available)

figure('Units','normalized','Position',[0.02,0.2,0.1,0.6])
unit_axes = axes; hold on;

% Plot units (depth vs normalized rate) with areas
if exist('probe_areas','var') && ~isempty('probe_areas')
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
    image(unit_axes,[0.5],probe_areas_image_depth,probe_areas_image);
    yline(unique(probe_areas_boundaries(:)),'color','k','linewidth',1);
    set(unit_axes,'YTick',probe_areas_centers,'YTickLabels',probe_areas{1}.acronym);
end

yyaxis right; 
set(gca,'YDir','reverse');
norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));

unit_dots = scatter( ...
    norm_spike_n,template_depths(unique(spike_templates)),20,'k','filled');
multiunit_lines = arrayfun(@(x) line(xlim,[0,0],'linewidth',2,'visible','off'),1:2);
xlim(unit_axes,[-0.1,1]);
ylim([-50, 3870]);
ylabel('Depth (\mum)')
xlabel('Normalized log rate')

linkprop([unit_axes.YAxis],'limits');
[unit_axes.YAxis.Color] = deal('k');
axis tight;




