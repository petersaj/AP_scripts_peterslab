
% Set animal/day/units
animal = 'AP024';
rec_day = '2024-09-20';
plot_units = [316];

% Set raster time bins
raster_window = [-0.2,0.7];
psth_bin_size = 0.001;
raster_t_bins = raster_window(1):psth_bin_size:raster_window(2);
raster_t = raster_t_bins(1:end-1) + diff(raster_t_bins)./2;

h = gobjects(length(plot_units),1);
for curr_unit = 1:length(plot_units)
    figure('color','w','name',sprintf('Unit %d',plot_units(curr_unit)));
    h(curr_unit) = tiledlayout(3,1);
end

for workflow_idx = 1:2

    switch workflow_idx
        case 1
            rec = plab.find_recordings(animal,rec_day,'*wheel*');
            rec_time = rec.recording{end};
        case 2
            rec = plab.find_recordings(animal,rec_day,'*lcr*');
            rec_time = rec.recording{end};
    end
    verbose = true;
    load_parts.ephys = true;
    ap.load_recording;

    % Get alignment times to use
    use_align = [];
    if contains(bonsai_workflow,'lcr')
        use_align = stimOn_times([trial_events.values(1:n_trials).TrialStimX] == 90);
        trial_sort = 1:length(use_align);
    elseif contains(bonsai_workflow,'wheel')
        use_align = stimOn_times(1:n_trials);
        [~,trial_sort] = sort(stim_to_move(1:n_trials));
    end

    % If alignment times empty (e.g. not this modality), skip
    if isempty(use_align)
        continue
    end

    % Get psth/raster of spikes to plot
    [use_spikes,spike_groups] = ismember(spike_templates,plot_units);
    [psth,raster,raster_t] = ap.ephys_psth(spike_times_timelite(use_spikes), ...
        use_align,spike_groups(use_spikes));

    % Plot PSTHs and rasters
    for curr_unit_idx = 1:length(plot_units)

        % Plot PSTH
        nexttile(h(curr_unit_idx),tilenum(h(curr_unit_idx),1,1)); hold on;
        plot(raster_t,smoothdata(psth(curr_unit_idx,:),2,'gaussian',50)','linewidth',2);
        xlim(raster_window)
        axis off;
        xline(0,'r');

        % Plot raster
        nexttile(h(curr_unit_idx),tilenum(h(curr_unit_idx),min(workflow_idx+1,3),1)); hold on;
        [raster_y,raster_x] = find(raster(trial_sort,:,curr_unit_idx));
        plot(raster_t(raster_x),raster_y,'.k'); % (dots)
%         plot(raster_t(raster_x),raster_y,'|k'); % (lines)
        xlim(raster_window)
        axis off
        set(gca,'YDir','reverse');
        xline(0,'r');

        drawnow;
    end
end


% Label panels and link axes
for curr_unit = 1:length(plot_units)
    curr_h = h(curr_unit);
    title(curr_h,sprintf('%s %s %d',animal,rec_day,plot_units(curr_unit)));
    title(nexttile(curr_h,tilenum(curr_h,1,1)),'Visual');
    title(nexttile(curr_h,tilenum(curr_h,2,1)),'Task');
    title(nexttile(curr_h,tilenum(curr_h,3,1)),'Passive');
    legend(nexttile(curr_h,tilenum(curr_h,1,1)),{'Task','','Passive'});

    linkaxes(curr_h.Children([2,3]));
end









