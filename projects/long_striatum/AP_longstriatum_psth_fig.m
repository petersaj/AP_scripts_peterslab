function AP_longstriatum_psth_fig(animal,rec_day,plot_unit,task_flag,parent_h)
% AP_longstriatum_psth_fig(animal,rec_day,plot_units,task_flag,tile_h)
%
% animal/rec_day/plot_units: recording coordinate to plot
% task_flag: plot task (true/false)
% parent_h: handle of parent container (figure or tiledlayout)
%
% Plot task and passive PSTH/rasters for longstriatum data

% Plot task by default
if ~exist('task_flag','var') || isempty(task_flag)
    task_flag = true;
end

% If no parent given, make figure
if ~exist('parent_h','var') || isempty(parent_h)
    parent_h = figure('color','w');
end

% Make tiledlayout
h = tiledlayout(parent_h,4+task_flag,1);
% (if parent is tiled layout, set sub-layout as next available)
if isa(parent_h,'matlab.graphics.layout.TiledChartLayout')
    h.Layout.Tile = length(parent_h.Children);
end

task_stim_colors = [0.7,0.7,0;ap.colormap('BKR',3)];
plot_colororder = task_stim_colors([task_flag;true(3,1)],:);

all_workflows = {'*wheel*','*lcr*'};
plot_workflows = all_workflows([task_flag,true]);

for curr_workflow = plot_workflows

    rec = plab.find_recordings(animal,rec_day,curr_workflow);
    rec_time = rec.recording{end};

    verbose = false;
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

        % % (plot quiescent trials only)
        % use_align = cellfun(@(x) ...
        %     stimOn_times(stim_x(1:length(stimOn_times)) == x & quiescent_trials), ...
        %     num2cell(unique(stim_x)),'uni',false);

        % (plot all trials)
        use_align = cellfun(@(x) ...
            stimOn_times(stim_x(1:length(stimOn_times)) == x), ...
            num2cell(unique(stim_x)),'uni',false);

        trial_sort = cellfun(@(x) 1:length(x),use_align,'uni',false);

    elseif contains(bonsai_workflow,'wheel')
        use_trials = trial_outcome == 1;

        use_align = stimOn_times(use_trials);
        [~,trial_sort] = sort(stim_to_move(use_trials));
    end

    % If alignment times empty (e.g. not this modality), skip
    if isempty(use_align)
        continue
    end

    % Get psth/raster of spikes to plot
    raster_window = [-0.5,1.5];
    psth_smoothing = 100;
    [use_spikes,spike_groups] = ismember(spike_templates,plot_unit);
    [psth,raster,raster_t] = ap.psth(spike_times_timelite(use_spikes), ...
        use_align,spike_groups(use_spikes), ...
        'window',raster_window,'smoothing',psth_smoothing);

    % Plot PSTH
    nexttile(h,tilenum(h,1,1)); hold on;
    set(gca,'ColorOrder',plot_colororder);
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
        nexttile(h); hold on;
        [raster_y,raster_x] = find(raster{curr_raster}(trial_sort{curr_raster},:));
        plot(raster_t(raster_x),raster_y,'.k');
        xlim(raster_window)
        axis off
        set(gca,'YDir','reverse');
        xline(0,'r');

        if contains(bonsai_workflow,'lcr')
            title('Passive');
        elseif contains(bonsai_workflow,'wheel')
            title('Task');
        end
    end

    drawnow;
end

% Label panels and link axes
title(h,sprintf('%s %s %d',animal,rec_day,plot_unit));
linkaxes(h.Children(1:end-1));










