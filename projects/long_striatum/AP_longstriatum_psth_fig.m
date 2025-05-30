
% Set animal/day/units

% % (TAN)
% animal = 'AP025';
% rec_day = '2024-09-11';
% plot_units = [270];

% % (TAN)
% animal = 'AP023';
% rec_day = '2024-09-09';
% plot_units = [289,268];

% (ad hoc)
animal = 'AM022';
rec_day = '2024-04-05';
plot_units = [187];

% Make figures
h = gobjects(length(plot_units),1);
for curr_unit = 1:length(plot_units)
    figure('color','w','name',sprintf('Unit %d',plot_units(curr_unit)));
    h(curr_unit) = tiledlayout(5,1);
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

        % (get only quiescent trials)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');

        stim_x = vertcat(trial_events.values.TrialStimX);

        use_align = cellfun(@(x) ...
            stimOn_times(stim_x(1:length(stimOn_times)) == x & quiescent_trials), ...
            num2cell(unique(stim_x)),'uni',false);
        trial_sort = cellfun(@(x) 1:length(x),use_align,'uni',false);

    elseif contains(bonsai_workflow,'wheel')
        use_align = stimOn_times(1:n_trials);
        [~,trial_sort] = sort(stim_to_move(1:n_trials));
    end

    % If alignment times empty (e.g. not this modality), skip
    if isempty(use_align)
        continue
    end

    % Get psth/raster of spikes to plot
    raster_window = [-0.5,1.5];
    psth_smoothing = 100;
    [use_spikes,spike_groups] = ismember(spike_templates,plot_units);
    [psth,raster,raster_t] = ap.psth(spike_times_timelite(use_spikes), ...
        use_align,spike_groups(use_spikes), ...
        'window',raster_window,'smoothing',psth_smoothing);

    % Plot PSTHs and rasters
    for curr_unit_idx = 1:length(plot_units)

        % Plot PSTH
        nexttile(h(curr_unit_idx),tilenum(h(curr_unit_idx),1,1)); hold on;
        plot(raster_t,permute(psth,[2,3,1]),'linewidth',2);
        xlim(raster_window)
        axis off;
        xline(0,'r');

        % Plot raster
        if ~iscell(raster)
            raster = {raster};
            trial_sort = {trial_sort};
        end
        for curr_raster = 1:length(raster)
            nexttile(h(curr_unit_idx),tilenum(h(curr_unit_idx),min(workflow_idx+curr_raster),1)); hold on;
            [raster_y,raster_x] = find(raster{curr_raster}(trial_sort{curr_raster},:,curr_unit_idx));
            plot(raster_t(raster_x),raster_y,'.k'); 
            xlim(raster_window)
            axis off
            set(gca,'YDir','reverse');
            xline(0,'r');
        end

        drawnow;
    end
end


% Label panels and link axes
for curr_unit = 1:length(plot_units)
    curr_h = h(curr_unit);
    title(curr_h,sprintf('%s %s %d',animal,rec_day,plot_units(curr_unit)));
    title(nexttile(curr_h,tilenum(curr_h,2,1)),'Task');
    title(nexttile(curr_h,tilenum(curr_h,3,1)),'Passive');
    legend(nexttile(curr_h,tilenum(curr_h,1,1)),{'Task','','Passive'});

    linkaxes(curr_h.Children(2:end-1));
end









