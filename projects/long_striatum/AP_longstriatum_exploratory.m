%% Batch behavior

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022', ...
    'AM026','AM027'};

% Create master tiled layout
figure;
t = tiledlayout(1,length(animals),'TileSpacing','tight');

% Grab learning day for each mouse
bhv = struct;

for curr_animal_idx = 1:length(animals)

    animal = animals{curr_animal_idx};

    use_workflow = {'stim_wheel*'};
%     use_workflow = {'*audio_volume*'};
%     use_workflow = {'*no_change*'};
    recordings = plab.find_recordings(animal,[],use_workflow);

    surround_time = [-5,5];
    surround_sample_rate = 100;
    surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

    n_trials_water = nan(length(recordings),2);
    frac_move_day = nan(length(recordings),1);
    rxn_med = nan(length(recordings),1);
    frac_move_stimalign = nan(length(recordings),length(surround_time_points));
    rxn_stat_p = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get total trials/water
        n_trials_water(curr_recording,:) = [length(trial_events.timestamps), ...
            sum(([trial_events.values.Outcome] == 1)*6)];

        % Get median stim-outcome time
        n_trials = length([trial_events.timestamps.Outcome]);
        rxn_med(curr_recording) = median(seconds([trial_events.timestamps(1:n_trials).Outcome] - ...
            cellfun(@(x) x(1),{trial_events.timestamps(1:n_trials).StimOn})));

        % Align wheel movement to stim onset
        align_times = stimOn_times;
        pull_times = align_times + surround_time_points;

        frac_move_day(curr_recording) = nanmean(wheel_move);

        event_aligned_wheel_vel = interp1(timelite.timestamps, ...
            wheel_velocity,pull_times);
        event_aligned_wheel_move = interp1(timelite.timestamps, ...
            +wheel_move,pull_times,'previous');

        frac_move_stimalign(curr_recording,:) = nanmean(event_aligned_wheel_move,1);

        % Get association stat
        % (skip if only a few trials)
        if n_trials < 10
            continue
        end

        rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move);

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});
        AP_print_progress_fraction(curr_recording,length(recordings));

    end

    % Define learned day from reaction stat p-value and reaction time
    learned_day = rxn_stat_p < 0.05 & rxn_med < 5;

    relative_day = days(datetime({recordings.day}) - datetime({recordings(1).day}))+1;
    nonrecorded_day = setdiff(1:length(recordings),relative_day);

    % Draw in tiled layout nested in master
    t_animal = tiledlayout(t,4,1);
    t_animal.Layout.Tile = curr_animal_idx;
    title(t_animal,animal);

    nexttile(t_animal);
    yyaxis left; plot(relative_day,n_trials_water(:,1));
    ylabel('# trials');
    yyaxis right; plot(relative_day,frac_move_day);
    ylabel('Fraction time moving');
    xlabel('Day');
    if any(nonrecorded_day)
        xline(nonrecorded_day,'--k');
    end
    if any(learned_day)
        xline(relative_day(learned_day),'g');
    end

    nexttile(t_animal);
    yyaxis left
    plot(relative_day,rxn_med)
    set(gca,'YScale','log');
    ylabel('Med. rxn');
    xlabel('Day');
    if any(nonrecorded_day)
        xline(nonrecorded_day,'--k');
    end
    if any(learned_day)
        xline(relative_day(learned_day),'g');
    end

    yyaxis right
    prestim_max = max(frac_move_stimalign(:,surround_time_points < 0),[],2);
    poststim_max = max(frac_move_stimalign(:,surround_time_points > 0),[],2);
    stim_move_frac_ratio = (poststim_max-prestim_max)./(poststim_max+prestim_max);
    plot(relative_day,stim_move_frac_ratio);
    yline(0);
    ylabel('pre/post move idx');
    xlabel('Day');

    nexttile(t_animal);
    imagesc(surround_time_points,[],frac_move_stimalign); hold on;
    clim([0,1]);
    colormap(gca,AP_colormap('WK'));
    set(gca,'YTick',1:length(recordings),'YTickLabel', ...
        cellfun(@(day,num) sprintf('%d (%s)',num,day(6:end)), ...
        {recordings.day},num2cell(1:length(recordings)),'uni',false));
    xlabel('Time from stim');
    if any(learned_day)
        plot(0,find(learned_day),'.g')
    end

    nexttile(t_animal); hold on
    set(gca,'ColorOrder',copper(length(recordings)));
    plot(surround_time_points,frac_move_stimalign','linewidth',2);
    xline(0,'color','k');
    ylabel('Fraction moving');
    xlabel('Time from stim');
    if any(learned_day)
        AP_errorfill(surround_time_points,frac_move_stimalign(learned_day,:)', ...
            0.02,[0,1,0],0.1,false);
        
        % Store behavior across animals
        bhv(curr_animal_idx).rxn_med = rxn_med;
        bhv(curr_animal_idx).stim_move_frac_ratio = stim_move_frac_ratio;
        bhv(curr_animal_idx).learned_day = learned_day;

    end

    drawnow;

end

% Save
save_dir = 'C:\Users\petersa\Desktop';
save_fn = fullfile(save_dir,'bhv');
save(save_fn,'bhv')


%% Batch widefield (master U)

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022', ...
    'AM026','AM027'};

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.03;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

day_V_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
   
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
%     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    n_components = 2000;
    n_align = 3; % (just hardcoded for now)
    day_V = nan(n_components,length(t_centers),n_align,length(recordings));

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % If widefield isn't aligned, skip
        if ~load_parts.widefield_align
            continue
        end

        % Stim times (quiescent trials only)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        if isfield(trial_events.values,'TrialStimX')
            stim_x = vertcat(trial_events.values.TrialStimX);
            % sometimes not all trials shown?
            stim_x = stim_x(1:length(stimOn_times));
            align_times = cellfun(@(x) stimOn_times(stim_x == x & quiescent_trials), ...
                num2cell(unique(stim_x)),'uni',false);
        else
            % Stim times (task)
            align_times = {stimOn_times,stim_move_time,reward_times};
        end

        % Align ROI trace to align times
        aligned_V_mean = nan(n_components,length(t_centers),length(align_times));
        for curr_align = 1:length(align_times)
            peri_event_t = align_times{curr_align} + t_centers;

            aligned_V = interp1(wf_t,wf_V',peri_event_t);

            aligned_V_baselinesub = aligned_V - ...
                mean(aligned_V(:,t_centers < 0,:),2);

            aligned_V_mean(:,:,curr_align) = ...
                permute(nanmean(aligned_V_baselinesub,1),[3,2,1]);
        end

        % Store mean aligned ROI trace
        day_V(:,:,:,curr_recording) = aligned_V_mean;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_V_all{curr_animal} = day_V;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Save
save_dir = 'C:\Users\petersa\Desktop';
save_fn = fullfile(save_dir,'wf_passive');
save(save_fn,'day_V_all')



%% Batch MUA by depth (size depths)
% (same division as cortical maps below)

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022', ...
    'AM026','AM027'};

plot_depths = false;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set MUA length (microns)
mua_length = 200;

day_mua_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
%     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_mua = cell(length(recordings),1);

    if plot_depths
        figure;
        h = tiledlayout(1,length(recordings));
        title(h,animal);
    end

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Get striatum start = lowest unit density, end = end of probe
        unit_density_bins = 0:100:3840;
        unit_density = histcounts(template_depths,unit_density_bins);
        [~,unit_density_min_bottom_idx] = min(fliplr(unit_density));
        unit_density_min_idx = length(unit_density_bins) - unit_density_min_bottom_idx;
        template_depths_sorted = sort(template_depths);
        str_start =  template_depths_sorted(find(template_depths_sorted >= ...
            unit_density_bins(unit_density_min_idx+1),1));
        str_end = max(channel_positions(:,2));

        % Discretize spikes by depth
        depth_group_edges = str_start:mua_length:str_end;
        if length(depth_group_edges) < 2
            continue    
        end
        spike_depth_group = discretize(spike_depths,depth_group_edges);
        n_depths = max(spike_depth_group);

        % Plot units and depth groupings
        if plot_depths
            nexttile(h);
            set(gca,'YDir','reverse'); hold on
            norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
            unit_dots = scatter(norm_spike_n,template_depths( ...
                unique(spike_templates)),20,'k','filled');
            yline(depth_group_edges,'color','r');
            drawnow;
        end

        % Stim times (quiescent trials only)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        if isfield(trial_events.values,'TrialStimX')
            stim_x = vertcat(trial_events.values.TrialStimX);
            % sometimes not all trials shown?
            stim_x = stim_x(1:length(stimOn_times));
            align_times = cellfun(@(x) stimOn_times(stim_x == x & quiescent_trials), ...
                num2cell(unique(stim_x)),'uni',false);
        else
            % Stim times (task)
            % (skip first trial - somes very short iti?)
            % (also skip stim_to_move negative - bad quiescence))
            use_trials = true(n_trials,1);
            use_trials(1) = false;
            use_trials(stim_to_move < 0) = false;

            align_times = {stimOn_times(use_trials),stim_move_time(use_trials),reward_times_task};
        end

        % Get PSTH depth x t x align
        depth_psth = nan(n_depths,length(t_bins)-1,length(align_times));
        for curr_align = 1:length(align_times)

            t_peri_event = align_times{curr_align} + t_bins;

            use_spikes = spike_times_timelite >= min(t_peri_event,[],'all') & ...
                spike_times_timelite <= max(t_peri_event,[],'all') & ...
                ~isnan(spike_depth_group);

            spikes_binned_continuous = histcounts2(spike_times_timelite(use_spikes),spike_depth_group(use_spikes), ...
                reshape(t_peri_event',[],1),1:n_depths+1)./psth_bin_size;

            use_continuous_bins = reshape(padarray(true(size(t_peri_event(:,1:end-1)')),[1,0],false,'post'),[],1);
            spikes_binned = spikes_binned_continuous(use_continuous_bins,:);

            depth_psth(:,:,curr_align) = ...
                nanmean(permute(reshape(spikes_binned, ...
                size(t_peri_event,2)-1,[],n_depths),[3,1,2]),3);
        end

        smooth_size = 100;
        depth_psth_smooth = smoothdata(depth_psth,2,'gaussian',smooth_size);

        % Store MUA (raw - to combine and normalize later)
        day_mua{curr_recording} = depth_psth_smooth;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_mua_all{curr_animal} = day_mua;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Save
save_dir = 'C:\Users\petersa\Desktop';
save_fn = fullfile(save_dir,'mua_passive');
save(save_fn,'day_mua_all')


%% Get visually responsive units by size depth 

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022', ...
    'AM026','AM027'};

plot_depths = false;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set MUA length (microns)
mua_length = 200;

stim_responsive_cells = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
%     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings([recordings.ephys]);

    recording_idx = 1:length(recordings);

    stim_responsive_cells_day = cell(length(recordings),1);

    if plot_depths
        figure;
        h = tiledlayout(1,length(recordings));
        title(h,animal);
    end

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts = struct;
        load_parts.ephys = true;
        ap.load_recording;

        % Get striatum start = lowest unit density, end = end of probe
        unit_density_bins = 0:100:3840;
        unit_density = histcounts(template_depths,unit_density_bins);
        [~,unit_density_min_bottom_idx] = min(fliplr(unit_density));
        unit_density_min_idx = length(unit_density_bins) - unit_density_min_bottom_idx;
        template_depths_sorted = sort(template_depths);
        str_start =  template_depths_sorted(find(template_depths_sorted >= ...
            unit_density_bins(unit_density_min_idx+1),1));
        str_end = max(channel_positions(:,2));

        % Discretize templates by depth
        depth_group_edges = str_start:mua_length:str_end;
        if length(depth_group_edges) < 2
            continue
        end
        template_depth_group = discretize(template_depths,depth_group_edges);
        n_depths = max(template_depth_group);

        % Get responsive units
        % (get quiescent trials)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');

        % (for vis passive: right-side stim)
        stim_type = vertcat(trial_events.values.TrialStimX);
        use_align = stimOn_times(stim_type(1:length(stimOn_times)) == 90 & quiescent_trials);

        baseline_t = [-0.2,0];
        response_t = [0,0.2];

        baseline_bins = use_align + baseline_t;
        response_bins = use_align + response_t;

        event_bins = [baseline_bins,response_bins];
        spikes_binned_continuous = histcounts2(spike_times_timelite,spike_templates, ...
            reshape([baseline_bins,response_bins]',[],1),1:size(templates,1)+1);

        event_spikes = permute(reshape(spikes_binned_continuous(1:2:end,:),2, ...
            size(event_bins,1),[]),[2,1,3]);

        event_response = squeeze(mean(diff(event_spikes,[],2),1));

        n_shuff = 1000;
        event_response_shuff = cell2mat(arrayfun(@(shuff) ...
            squeeze(mean(diff(ap.shake(event_spikes,2),[],2),1)), ...
            1:n_shuff,'uni',false));

        event_response_rank = tiedrank(horzcat(event_response,event_response_shuff)')';
        event_response_p = event_response_rank(:,1)./(n_shuff+1);

        responsive_cells = event_response_p > 0.95;

        % % Plot responsive units by depth
        % unit_dots = ap.plot_unit_depthrate(spike_templates,template_depths,probe_areas);
        % unit_dots.CData = +([1,0,0].*(event_response_p > 0.95)) + ([0,0,1].*(event_response_p < 0.05));

        depth_responsive_n = accumarray(template_depth_group(~isnan(template_depth_group)), ...
            responsive_cells(~isnan(template_depth_group)),[n_depths,1],@sum);

        depth_total_n = accumarray(template_depth_group(~isnan(template_depth_group)), ...
            responsive_cells(~isnan(template_depth_group)),[n_depths,1],@length);

        % Store responsive cells
        stim_responsive_cells_day{curr_recording} = [depth_responsive_n,depth_total_n];

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    stim_responsive_cells{curr_animal} = stim_responsive_cells_day;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Save
save_dir = 'C:\Users\petersa\Desktop';
save_fn = fullfile(save_dir,'stim_cells');
save(save_fn,'stim_responsive_cells')


%% Batch cortical maps by striatal depth

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022', ...
    'AM026','AM027'};

% Set MUA length (microns)
mua_length = 200;

ctx_map_all = cell(length(animals),1);
ctx_expl_var_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
%     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    ctx_map = cell(length(recordings),1);
    ctx_expl_var = cell(length(recordings),1);
    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts.ephys = true;
        load_parts.widefield = true;
        load_parts.widefield_align = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % If not aligned, skip
        if ~load_parts.widefield_align 
            continue
        end

        % Whole thing in a try/catch
        try

            % Get striatum start = lowest unit density, end = end of probe
            unit_density_bins = 0:100:3840;
            unit_density = histcounts(template_depths,unit_density_bins);
            [~,unit_density_min_bottom_idx] = min(fliplr(unit_density));
            unit_density_min_idx = length(unit_density_bins) - unit_density_min_bottom_idx;
            template_depths_sorted = sort(template_depths);
            str_start =  template_depths_sorted(find(template_depths_sorted >= ...
                unit_density_bins(unit_density_min_idx+1),1));
            str_end = max(channel_positions(:,2));

            % Discretize spikes by depth
            depth_group_edges = str_start:mua_length:str_end;
            depth_group = discretize(spike_depths,depth_group_edges);

            % Get MUA binned by widefield
            sample_rate = (1/mean(diff(wf_t)));
            skip_seconds = 60;
            time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
            time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

            binned_spikes = zeros(max(depth_group),length(time_bins)-1);
            for curr_depth = 1:max(depth_group)
                curr_spike_times = spike_times_timelite(depth_group == curr_depth);
                binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
            end

            binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
            binned_spikes_std(isnan(binned_spikes_std)) = 0;

            use_svs = 1:500;
            kernel_t = [0,0];
            kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
            lambda = 10;
            zs = [false,false];
            cvfold = 5;
            return_constant = false;
            use_constant = true;

            fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

            % Regress cortex to spikes
            [k,predicted_spikes,explained_var] = ...
                AP_regresskernel(fVdf_deconv_resample, ...
                binned_spikes_std,kernel_frames, ...
                lambda,zs,cvfold,return_constant,use_constant);

            ctx_map{curr_recording} = permute(k,[1,3,2]);
            ctx_expl_var{curr_recording} = explained_var.total;

        catch me
            % Any errors: do nothing
        end

        % Prep for next loop
        clearvars('-except',preload_vars{:});

        ap.print_progress_fraction(curr_recording,length(recordings));

    end

    ctx_map_all{curr_animal} = ctx_map;
    ctx_expl_var_all{curr_animal} = ctx_expl_var;

    fprintf('Done %s\n',animal);

end

% Save
save_dir = 'C:\Users\petersa\Desktop';
save_fn = fullfile(save_dir,'ctx_maps_passive');
save(save_fn,'ctx_map_all','ctx_expl_var_all');


%% Batch cortical regression by striatal depth

animals = {'AM008','AP008','AP009', ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022', ...
    'AM026','AM027'};

plot_depths = false;

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set depth groups
n_depths = 3;

day_spikes_all = cell(length(animals),1);
day_predicted_spikes_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    use_workflow = 'lcr_passive';
%     use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);
    recordings = recordings(cellfun(@any,{recordings.widefield}) & [recordings.ephys]);

    recording_idx = 1:length(recordings);

    day_spikes = nan(n_depths,length(t_centers),length(recordings));
    day_predicted_spikes = nan(n_depths,length(t_centers),length(recordings));
  
    if plot_depths
        figure;
        h = tiledlayout(1,length(recordings));
        title(h,animal);
    end

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};

        load_parts.ephys = true;
        load_parts.widefield = true;
        load_parts.widefield_align = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % If not aligned, skip
        if ~load_parts.widefield_align 
            continue
        end

        % Get striatum start = lowest unit density, end = end of probe
        unit_density_bins = 0:100:3840;
        unit_density = histcounts(template_depths,unit_density_bins);
        [~,unit_density_min_bottom_idx] = min(fliplr(unit_density));
        unit_density_min_idx = length(unit_density_bins) - unit_density_min_bottom_idx;
        template_depths_sorted = sort(template_depths);
        str_start =  template_depths_sorted(find(template_depths_sorted >= ...
            unit_density_bins(unit_density_min_idx+1),1));
        str_end = max(channel_positions(:,2));

        % Discretize spikes by depth
        depth_group_edges = round(linspace(str_start,str_end,n_depths+1));
        depth_group = discretize(spike_depths,depth_group_edges);

        % Get MUA binned by widefield
        sample_rate = (1/mean(diff(wf_t)));
        skip_seconds = 60;
        time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

        binned_spikes = zeros(n_depths,length(time_bins)-1);
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timelite(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end

        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;

        use_svs = 1:100;
        kernel_t = [-0.1,0.1];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        lambda = 5;
        zs = [false,false];
        cvfold = 5;
        return_constant = true;
        use_constant = true;

        fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

        % Regress cortex to spikes
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample, ...
            binned_spikes_std,kernel_frames, ...
            lambda,zs,cvfold,return_constant,use_constant);

        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});

        % Stim times (quiescent trials only)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        if isfield(trial_events.values,'TrialStimX')
            stim_x = vertcat(trial_events.values.TrialStimX);
            % temp - what happened here, not all trials shown?
            stim_x = stim_x(1:length(stimOn_times));
            align_times = stimOn_times(stim_x == 90 & quiescent_trials);
        else
            % Stim times (task)
            align_times = stimOn_times;
        end

        % Align measured and predicted activity to align times    
        peri_event_t = align_times + t_centers;
        
        aligned_spikes = reshape(interp1(time_bin_centers,binned_spikes', ...
            peri_event_t), ...
            length(align_times),length(t_centers),[]);

        aligned_predicted_spikes = reshape(interp1(time_bin_centers, ...
            predicted_spikes',peri_event_t), ...
            length(align_times),length(t_centers),[]);


        % Store mean aligned ROI trace
        day_spikes(:,:,curr_recording) = permute(nanmean(aligned_spikes,1),[3,2,1]);
        day_predicted_spikes(:,:,curr_recording) = permute(nanmean(aligned_predicted_spikes,1),[3,2,1]);

        % Prep for next loop
        clearvars('-except',preload_vars{:});

    end

    day_spikes_all{curr_animal} = day_spikes;
    day_predicted_spikes_all{curr_animal} = day_predicted_spikes;

    ap.print_progress_fraction(curr_animal,length(animals));

end

% Plot average measured and predicted responses
softnorm = 0;
day_spikes_all_norm = cellfun(@(x) (x-nanmean(x(:,t_centers<0,:),2))./(nanmean(x(:,t_centers<0,:),2)+softnorm),day_spikes_all,'uni',false);
day_predicted_spikes_all_norm = cellfun(@(x) (x-nanmean(x(:,t_centers<0,:),2))./(nanmean(x(:,t_centers<0,:),2)+softnorm),day_predicted_spikes_all,'uni',false);

load("C:\Users\petersa\Desktop\am_bhv.mat");

use_depth = 1; 
use_t = t_centers > 0.05 & t_centers < 0.2;
a = cellfun(@(x) squeeze(mean(x(use_depth,use_t,:),2)),day_predicted_spikes_all_norm,'uni',false);
b = cell2mat(cellfun(@(x) permute(x(1,:,:),[3,2,1]),day_predicted_spikes_all_norm,'uni',false));

ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day},'uni',false);

[m,s,n] = grpstats(cell2mat(a),cell2mat(ld),{'mean','sem','numel'});
r = grpstats(b,cell2mat(ld),'mean');

ld_x = unique(cell2mat(ld));

figure;
subplot(1,2,1);
imagesc(t_centers,unique(vertcat(ld{:})),r);
subplot(1,2,2);errorbar(ld_x,m,s,'k','linewidth',2);
xline(0);
xlabel('Learned day');
ylabel('\DeltaFR/FR_0')


use_t = t_centers > 0.05 & t_centers < 0.2;
a = cellfun(@(x) squeeze(mean(x(use_depth,use_t,:),2)),day_spikes_all_norm,'uni',false);
b = cellfun(@(x) squeeze(mean(x(use_depth,use_t,:),2)),day_predicted_spikes_all_norm,'uni',false);

figure; tiledlayout(1,3);
nexttile; hold on;
cellfun(@(x,ld) plot(ld,x),a,ld); title('Measured');
nexttile; hold on;
cellfun(@(x,ld) plot(ld,x),b,ld); title('Predicted');

nexttile
hold on; axis equal
cellfun(@(x,y) plot(x,y),a,b);
cellfun(@(x,y,ld) plot(x(ld==0),y(ld==0),'.k','MarkerSize',20),a,b,ld);
axis tight;
line(xlim,xlim,'color','k');
xlabel('Measured');
ylabel('Predicted');


%% Load MUA / cortex maps and process

am_data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\longitudinal_striatum\data';

load(fullfile(am_data_path,'bhv.mat'));
load(fullfile(am_data_path,'mua_passive.mat'));
% load(fullfile(am_data_path,'mua_task.mat'));
load(fullfile(am_data_path,'ctx_maps_passive.mat'));
% load(fullfile(am_data_path,'ctx_maps_task.mat'));
load(fullfile(am_data_path,'wf_passive.mat'));
load(fullfile(am_data_path,'stim_cells.mat'));

% Load master U, convert V to px
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);
ctx_map_px = cellfun(@(x) cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:size(x,1)),x),x,'uni',false),ctx_map_all,'uni',false);

use_animals = ~cellfun(@isempty,ctx_map_all);

ld = cellfun(@(x) (1:length(x))'-find(x,1),{bhv.learned_day}','uni',false);

% (learned day)
animal_ld_idx_cell = cellfun(@(animal,ld,data) cellfun(@(ld,data) ...
    [repmat(animal,size(data,2),1),repmat(ld,size(data,2),1),(1:size(data,2))'], ...
    num2cell(ld),data,'uni',false), ...
    num2cell(find(use_animals)),ld(use_animals),ctx_map_all(use_animals),'uni',false);
animal_ld_idx = cell2mat(vertcat(animal_ld_idx_cell{:}));

% (cardinal day)
animal_day_idx_cell = cellfun(@(animal,ld,data) cellfun(@(ld,data) ...
    [repmat(animal,size(data,2),1),repmat(ld,size(data,2),1),(1:size(data,2))'], ...
    num2cell(1:length(ld))',data,'uni',false), ...
    num2cell(find(use_animals)),ld(use_animals),ctx_map_all(use_animals),'uni',false);
animal_day_idx = cell2mat(vertcat(animal_day_idx_cell{:}));


r = cell2mat(vertcat(day_mua_all{use_animals}));
softnorm = 1;
r_norm = (r-nanmean(r(:,200:400),[2,3]))./(nanmean(r(:,200:400),[2,3])+softnorm);
use_t = [500,700];
r2 = nanmean(r_norm(:,use_t(1):use_t(2)),2);


% Concatenate V's
V_cat = cat(4,day_V_all{:});

V_animal_day_idx = cell2mat(cellfun(@(x,animal,ld) ...
    [repmat(animal,size(x,4),1),ld], ...
    day_V_all(use_animals),num2cell(find(use_animals)),ld(use_animals),'uni',false));

% Concatenate maps
c2 = cell2mat(permute(vertcat(ctx_map_px{:}),[2,3,1]));

ap.imscroll(c2); 
axis image off;
clim(max(abs(clim)).*[-1,1].*0.5);
colormap(ap.colormap('BWR',[],1.5));
ap.wf_draw('ccf','k');



% (draw ROI over mVis)
vis_map_weights = roi.trace';

% (draw ROI over mpfc)
mpfc_map_weights = roi.trace';




% Get vis cortex maps, group striatal responses, plot by LD
ctx_thresh = 0.005; %0.005
curr_use_ctx = vis_map_weights >= ctx_thresh & mpfc_map_weights < ctx_thresh; % vis-only
% curr_use_ctx = vis_map_weights < ctx_thresh & mpfc_map_weights >= ctx_thresh; % mpfc-only
% curr_use_ctx = vis_map_weights >= ctx_thresh & mpfc_map_weights >= ctx_thresh; % vis+mpfc
% curr_use_ctx = vis_map_weights >= ctx_thresh; % vis (all)
% curr_use_ctx = mpfc_map_weights >= ctx_thresh; % mpfc (all)

figure;plot(vis_map_weights,mpfc_map_weights,'.k')
xline(ctx_thresh,'r');yline(ctx_thresh,'r');
xlabel('vis weights');
ylabel('mpfc weights');

[grp,~,grp_idx] = unique(animal_ld_idx(curr_use_ctx,1:2),'rows');

use_align = 3;
x1 = ap.groupfun(@sum,r(curr_use_ctx,:,use_align),grp_idx,[]);
x1_norm = (x1-nanmean(x1(:,1:500),2))./nanmean(x1(:,1:500),2);

x11t = ap.groupfun(@nanmean,x1_norm,grp(:,2),[]);
x11 = max(x1_norm(:,use_t(1):use_t(2)),[],2);

figure;imagesc([],unique(grp(:,2)),x11t)
ylabel('LD');
title('MUA');

x2 = ap.groupfun(@nanmean,c2(:,:,curr_use_ctx),[],[],grp_idx);
x22 = ap.groupfun(@nanmean,x2,[],[],grp(:,2));
figure; tiledlayout('flow','TileSpacing','none');
c = max(abs(x22),[],'all').*[-1,1].*0.8;
ld_x = unique(grp(:,2));
for i = -3:3
    nexttile; 
    imagesc(x22(:,:,ld_x==i));
    axis image off;
    clim(c);
    ap.wf_draw('ccf','k');
    title(i);
end
colormap(ap.colormap('BWR',[],1.5));

ap.imscroll(x22,ld_x);
axis image off;
clim(c);
ap.wf_draw('ccf','k');
colormap(ap.colormap('BWR',[],1.5));


figure; hold on;
set(gca,'colororder',jet(length(unique(grp(:,1)))));
arrayfun(@(x) plot(grp(grp(:,1)==x,2), x11(grp(:,1)==x)),unique(grp(:,1)));
xlabel('Learned day');
ylabel('MUA')

[m,e] = grpstats(x11,grp(:,2),{'mean','sem'});
errorbar(unique(grp(:,2)),m,e,'k','linewidth',2);


% Average widefield from animals/days used above
use_align = 3;
use_v = ismember(V_animal_day_idx,grp,'rows');

V_avg = permute(ap.groupfun(@nanmean,V_cat(:,:,use_align,use_v),[],[],[],V_animal_day_idx(use_v,2)),[1,2,4,3]);
px_avg = plab.wf.svd2px(U_master,V_avg);
ap.imscroll(px_avg);
axis image;
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));

use_frames = 15:28;
px_tavg = squeeze(max(px_avg(:,:,use_frames,:),[],3));
ap.imscroll(px_tavg,ld_x);
axis image;
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));


use_frames = 20:23;
px_tmax_cat = squeeze(mean(plab.wf.svd2px(U_master,V_cat(:,use_frames,use_align,:)),3));

ap.imscroll(px_tmax_cat);
axis image;
clim(max(abs(clim)).*[-1,1]);
ap.wf_draw('ccf','k');
colormap(ap.colormap('PWG',[],1.5));

% (draw_roi);
use_roi_trace = roi.trace(use_v);
use_V_animal_day_idx = V_animal_day_idx(use_v,:);

figure; hold on;
arrayfun(@(animal) ...
    plot(use_V_animal_day_idx(use_V_animal_day_idx(:,1) == animal,2), ...
    use_roi_trace(use_V_animal_day_idx(:,1) == animal),'color',[0.5,0.5,0.5]), ...
    1:max(use_V_animal_day_idx(:,1)));

m = ap.groupfun(@nanmean,use_roi_trace,[],use_V_animal_day_idx(:,2));
e = ap.groupfun(@AP_sem,use_roi_trace,[],use_V_animal_day_idx(:,2));
errorbar(ld_x,m,e,'k','linewidth',2);


% K-means cortex maps
n_k = 4;
[kidx,ck] = kmeans(reshape(c2,[],size(c2,3))',n_k,'Distance','correlation');
ck = reshape(ck',[size(c2,[1,2]),n_k]);

figure;tiledlayout(1,n_k,'tilespacing','none','tileindexing','columnmajor');
for i = 1:n_k
    nexttile;
    imagesc(ck(:,:,i));
    axis image off;
    clim(max(abs(clim)).*[-1,1]);
    colormap(ap.colormap('PWG',[],1.5));
    ap.wf_draw('ccf','k');
end

% Plot MUA by k-means cluster
use_align = 3;
figure;h = tiledlayout(1,n_k);
for curr_kidx = 1:n_k

    curr_use_ctx = kidx == curr_kidx;
    [grp,~,grp_idx] = unique(animal_ld_idx(curr_use_ctx,1:2),'rows');

    x2 = ap.groupfun(@mean,c2(:,:,curr_use_ctx),[],[],grp_idx);
    x22 = ap.groupfun(@mean,x2,[],[],grp(:,2));
    figure; h2 = tiledlayout('flow','TileSpacing','none');
    c = max(abs(x22),[],'all').*[-1,1].*0.8;
    ld_x = unique(grp(:,2));
    for i = ld_x'
        nexttile(h2);
        imagesc(x22(:,:,ld_x==i));
        axis image off;
        clim(c);
        ap.wf_draw('ccf','k');
        title(i);
    end
    colormap(ap.colormap('PWG',[],1.5));


    x1 = ap.groupfun(@sum,r(curr_use_ctx,:,use_align),grp_idx,[]);
    x1_norm = (x1-nanmean(x1(:,1:500),2))./nanmean(x1(:,1:500),2);
    x11 = max(x1_norm(:,use_t(1):use_t(2)),[],2);

    nexttile(h); hold on;
    arrayfun(@(x) plot(grp(grp(:,1)==x,2), x11(grp(:,1)==x)),unique(grp(:,1)));
    xlabel('Learned day');
    ylabel('MUA')

    [m,e] = grpstats(x11,grp(:,2),{'mean','sem'});
    errorbar(unique(grp(:,2)),m,e,'k');
    drawnow;

    %%% TESTING (for stats)
    ld_x = unique(grp(:,2));
    act_grid = nan(max(grp(:,1)),length(ld_x));
    for curr_entry = 1:size(x11,1)
        act_grid(grp(curr_entry,1),ld_x==grp(curr_entry,2)) = ...
            x11(curr_entry);
    end
    p = signrank(act_grid(:,ld_x == -2),act_grid(:,ld_x == -1));
    %%%%%
    
end

% Plot responsive cell fraction by k-means cluster
figure;h = tiledlayout(1,n_k);

stim_responsive_cells_cat = cell2mat(vertcat(stim_responsive_cells{use_animals}));
for curr_kidx = 1:n_k

    curr_use_ctx = kidx == curr_kidx;
    [grp,~,grp_idx] = unique(animal_ld_idx(curr_use_ctx,1:2),'rows');

    stim_cells_combine = ap.groupfun(@sum,stim_responsive_cells_cat(curr_use_ctx,:),grp_idx,[]);
    stim_cells_frac  = stim_cells_combine(:,1)./stim_cells_combine(:,2);

    nexttile(h); hold on;
    arrayfun(@(x) plot(grp(grp(:,1)==x,2), stim_cells_frac(grp(:,1)==x)),unique(grp(:,1)));
    xlabel('Learned day');
    ylabel('Frac stim-responsive')

    [m,e] = grpstats(stim_cells_frac,grp(:,2),{'mean','sem'});
    errorbar(unique(grp(:,2)),m,e,'k');
    drawnow;

end


% Ephys-only mice, by depth
use_mice = 1:3;
% use_mice = setdiff(4:13,5);

use_align = 3;
max_depths = 5;

a = cellfun(@(x,ld) cell2mat(cellfun(@(x) sum(x(1:min(size(x,1),max_depths),:,use_align)), ...
    x(1:min(length(x),length(ld))),'uni',false)),day_mua_all(use_mice),ld(use_mice),'uni',false);
b = cell2mat(cellfun(@(x) (x-nanmean(x(:,1:500),2))./nanmean(x(:,1:500),2),a,'uni',false));

ld_sub = cell2mat(cellfun(@(x,y) x(1:min(size(y,1),length(x))),ld(use_mice),a,'uni',false));

x = ap.groupfun(@mean,b,ld_sub,[]);
ld_y = unique(ld_sub);
figure;imagesc([],ld_y,x)
colormap(AP_colormap('BWK'));
clim([-2,2]);

b2 = max(b(:,use_t(1):use_t(2)),[],2);
x2 = ap.groupfun(@mean,b2,ld_sub,[]);
figure;plot(ld_y,x2);



% Plot single animal over days
plot_animal = 3;
figure;
h = tiledlayout(1,length(day_mua_all{plot_animal}));
title(h,plot_animal);
for curr_day = 1:length(day_mua_all{plot_animal})
    curr_activity = day_mua_all{plot_animal}{curr_day};

    softnorm = 1;
    curr_activity_norm = (curr_activity-nanmean(curr_activity(:,200:400,:),2))./ ...
        (nanmean(curr_activity(:,200:400,:),2)+softnorm);

    nexttile; hold on;
    spacing = 3;
    col = lines(size(curr_activity,3));
    for curr_align = 1:size(curr_activity_norm,3)
        AP_stackplot(curr_activity_norm(:,:,curr_align)',[],spacing,[],col(curr_align,:));
    end
    title(ld{plot_animal}(curr_day));

end
linkaxes(h.Children);














