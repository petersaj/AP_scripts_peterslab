%%%% CK ephys data

%% Plot recording locations

animals = {'HA016','HA017','HA018','HA019','HA020','AP036','AP037'};

animal_col = ap.colormap('tube',length(animals));
ccf_draw = ap.ccf_draw;
ccf_draw.draw_name('Caudoputamen');

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    probe_color = animal_col(curr_animal,:);
    ccf_draw.draw_probes_histology(animal,probe_color);
    drawnow;
end

%% Get lick behavior before ephys

animals = {'HA016','HA017','HA018','HA019','HA020','AP036','AP037'};

anticipatory_licks = cell(length(animals),1);
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    recordings = plab.find_recordings(animal);
    ephys_start_idx = find([recordings.ephys],1);

    % % Use big task if present, small if not
    % task_workflow = 'visual_operant_lick_two_stim_static_big_stim';
    % task_recordings = plab.find_recordings(animal,[],task_workflow);
    % if isempty(task_recordings)
    %     task_workflow = 'visual_operant_lick_two_stim_static';
    %     task_recordings = plab.find_recordings(animal,[],task_workflow);
    % end
    
    % Use all static
    task_workflow = 'visual*static*';
    task_recordings = plab.find_recordings(animal,[],task_workflow);

    use_recordings = find(datetime({task_recordings.day}) < datetime({recordings(ephys_start_idx).day}));

    for curr_recording = 1:length(use_recordings)

        rec_day = task_recordings(use_recordings(curr_recording)).day;
        rec_time = task_recordings(use_recordings(curr_recording)).recording{end};

        load_parts.bhv = true;
        ap.load_recording;

        % Get trial parameters
        n_trials = sum(cellfun(@(x) length(x) == 2,{trial_events.timestamps.StimOn}));
        trial_stim_x = vertcat(trial_events.values(1:n_trials).TrialX);
        trial_static_stim_time = vertcat(trial_events.values(1:n_trials).TrialStimStaticTime);
        trial_quiescence_time = vertcat(trial_events.values(1:n_trials).TrialQuiescence);

        if contains(bonsai_workflow,'move')
            % Moving stim

            % Get stim on/move times
            % (photodiode CS- = on/off for CS+ = on,pulse on move,off)
            stim_pd_n = (trial_stim_x == -90)*1 + (trial_stim_x == 90)*2;
            stim_pd_on_grouped = mat2cell(photodiode_on_times(1:sum(stim_pd_n)),stim_pd_n);
            stim_pd_off_grouped = mat2cell(photodiode_off_times(1:sum(stim_pd_n)),stim_pd_n);

            stimOn_times = cellfun(@(x) x(1), stim_pd_on_grouped);
            stim_move_times = cellfun(@(x) x(1), stim_pd_off_grouped);

            % Time to first lick after stim in center (CS- would-be)
            reward_available_times = stim_move_times + trial_events.parameters.StimMoveTime;
            reward_available_firstlick_times = interp1(lick_times,lick_times,reward_available_times,'next');

            lick_align = {stimOn_times,stim_move_times,reward_available_firstlick_times};

        elseif contains(bonsai_workflow,'static')
            % Static stim

            % Get stim on/reward available times
            stimOn_times = photodiode_on_times(1:n_trials);

            % Time to first lick after stim in center (CS- would-be)
            reward_available_times = stimOn_times + trial_static_stim_time(1:n_trials);
            reward_available_firstlick_times = interp1(lick_times,lick_times,reward_available_times,'next');

            lick_align = {stimOn_times,reward_available_firstlick_times};
        end

        % Plot aligned licks
        if false
            figure('name',sprintf('%s %s',animal,rec_day));
            h = tiledlayout(4,length(lick_align),'TileIndexing','ColumnMajor');

            lick_window = [-7,10];
            lick_binsize = 0.01;
            for curr_align = 1:length(lick_align)
                plot_align = lick_align{curr_align};

                [lick_psth_r,lick_raster_r,lick_t] = ap.psth(lick_times,plot_align(trial_stim_x == 90 & ~isnan(plot_align)),...
                    'window',lick_window,'bin_size',lick_binsize,'smoothing',50);

                [lick_psth_l,lick_raster_l] = ap.psth(lick_times,plot_align(trial_stim_x == -90 & ~isnan(plot_align)),...
                    'window',lick_window,'bin_size',lick_binsize,'smoothing',50);

                nexttile(h,tilenum(h,1,curr_align)); hold on;
                plot(lick_t,lick_psth_r,'r');
                plot(lick_t,lick_psth_l,'b');
                xline(0,'color',[0.7,0.7,0]);

                nexttile(h,tilenum(h,2,curr_align),[3,1]); hold on;

                [lick_trial_r,lick_t_raster_r_idx] = find(lick_raster_r);
                [lick_trial_l,lick_t_raster_l_idx] = find(lick_raster_l);

                lick_t_raster_r = lick_t(lick_t_raster_r_idx);
                lick_t_raster_l = lick_t(lick_t_raster_l_idx);

                plot(lick_t_raster_r,lick_trial_r,'.r');
                plot(lick_t_raster_l,lick_trial_l+size(lick_raster_r,1),'.b');
                hold on; set(gca,'YDir','reverse');
                xline(0,'color',[0.7,0.7,0]);
                axis tight

            end
            linkaxes(h.Children,'x');
            drawnow;
        end

        % Get anticipatory licks (1s before reward available)
        [~,cs_plus_lick] = ap.psth(lick_times,reward_available_firstlick_times( ...
            trial_stim_x == 90 & ~isnan(reward_available_firstlick_times)),...
            'window',repelem(-0.5,2,1),'bin_size',1);

        [~,cs_minus_lick] = ap.psth(lick_times,reward_available_firstlick_times( ...
            trial_stim_x == -90 & ~isnan(reward_available_firstlick_times)),...
            'window',repelem(-0.5,2,1),'bin_size',1);

        % % Get stim onset licks
        % [~,cs_plus_lick] = ap.psth(lick_times,stimOn_times( ...
        %     trial_stim_x == 90 & ~isnan(stimOn_times)),...
        %     'window',repelem(0.5,2,1),'bin_size',1);
        % 
        % [~,cs_minus_lick] = ap.psth(lick_times,stimOn_times( ...
        %     trial_stim_x == -90 & ~isnan(stimOn_times)),...
        %     'window',repelem(0.5,2,1),'bin_size',1);

        anticipatory_licks{curr_animal}(curr_recording,:) = {cs_plus_lick,cs_minus_lick};
    end
end

% % (difference: concat)
% x = cellfun(@(x) mean(vertcat(x{:,1}))-mean(vertcat(x{:,2})),anticipatory_licks);

% (CS+ lick: concat)
x = cellfun(@(x) mean(vertcat(x{:,1})),anticipatory_licks);

% % (CS+ lick: avg)
% x = cellfun(@(x) mean(cellfun(@mean,x(:,1))),anticipatory_licks);

% % (discrimination)
% x = cellfun(@(x) (mean(vertcat(x{:,1}))-mean(vertcat(x{:,2})))./(mean(vertcat(x{:,1}))+mean(vertcat(x{:,2}))),anticipatory_licks);

% % (discrimination)
% anticipatory_licks_mean = cellfun(@(x) cellfun(@mean,x),anticipatory_licks,'uni',false);
% x = cellfun(@(x) mean(-diff(x,[],2)./sum(x,2),1),anticipatory_licks_mean);

% % (frac non-zero)
% x = cellfun(@(x) mean(cellfun(@(x) mean(x>0),x(:,1)),1),anticipatory_licks);

figure;plot(x,'.k','MarkerSize',20);
ylabel('Antipatory licks');


%% Striatal visual responses

animals = {'HA016','HA017','HA018','HA019','HA020','AP036','AP037'};

use_recs = cell(size(animals));
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    histology_filepattern = plab.locations.filename('server',animal,[],[], ...
        'histology','**','AP_histology_processing.mat');
    histology_dir = dir(histology_filepattern);
    if ~isempty(histology_dir)
        histology_filename = fullfile(histology_dir.folder,histology_dir.name);
    end

    % Get days with mapped posterior recordings
    load(histology_filename)
    day_pattern = digitsPattern(4)+'-'+digitsPattern(2)+'-'+digitsPattern(2);
    mapped_paths = {AP_histology_processing.annotation.ephys_path};
    ephys_mapped_days = unique(extract(mapped_paths(~cellfun(@isempty,mapped_paths)),day_pattern));

    use_recs{curr_animal} = ephys_mapped_days;
end

% Set times for PSTH
raster_window = [-0.2,0.8];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Get PSTHs and responsive units
unit_ccf_all = cell(length(animals),1);
unit_psth_all = cell(length(animals),1);
responsive_units_all = cell(length(animals),1);
for curr_animal = 1:length(animals)
    for curr_day = 1:length(use_recs{curr_animal})

        % Load data
        animal = animals{curr_animal};

        % % Use big task/passive if present, small if not
        % task_workflow = 'visual_operant_lick_two_stim_static_big_stim';
        % task_recordings = plab.find_recordings(animal,[],task_workflow);
        % if ~isempty(task_recordings)
        %     passive_workflow = 'lcr_passive_corner_CS\+_big_stim';
        % else
        %     passive_workflow = 'lcr_passive_corner_CS\+';
        % end
        passive_workflow = 'lcr_passive_corner_CS\+_big_stim';

        rec_day = use_recs{curr_animal}{curr_day};
        rec_time = plab.find_recordings(animal,rec_day,passive_workflow).recording{end};

        ap.load_recording;

        % Get quiescent trials
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');

        % Get stim times to use
        stim_x = vertcat(trial_events.values.TrialStimX);

        n_trials = min(length(stimOn_times),length(stim_x));
        align_times = cellfun(@(x) ...
            stimOn_times(stim_x(1:n_trials) == x & ...
            quiescent_trials(1:n_trials)),num2cell(unique(stim_x)),'uni',false);

        % Get PSTH by 2D histogram
        n_units = size(templates,1);
        unit_psth = nan(n_units,length(t_bins)-1,length(align_times));
        for curr_align = find(~cellfun(@isempty,align_times))'
            t_peri_event = align_times{curr_align} + t_bins;

            use_spikes = spike_times_timelite >= min(t_peri_event,[],'all') & ...
                spike_times_timelite <= max(t_peri_event,[],'all');

            spikes_binned_continuous = histcounts2(spike_times_timelite(use_spikes),spike_templates(use_spikes), ...
                reshape(t_peri_event',[],1),1:size(templates,1)+1)./psth_bin_size;

            use_continuous_bins = reshape(padarray(true(size(t_peri_event(:,1:end-1)')),[1,0],false,'post'),[],1);
            spikes_binned = permute(reshape(spikes_binned_continuous(use_continuous_bins,:), ...
                size(t_peri_event,2)-1,size(t_peri_event,1),size(templates,1)),[3,1,2]);

            unit_psth(:,:,curr_align) = nanmean(spikes_binned,3);
        end

        % Get responsive units
        response_align = align_times{unique(stim_x) == 90}; % response = R vis

        baseline_t = [-0.2,0];
        response_t = [0,0.2];

        baseline_bins = response_align + baseline_t;
        response_bins = response_align + response_t;

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

        responsive_units = event_response_p > 0.95;

        % Store unit PSTHs and responsive units (striatum units only)
        striatum_probe_area_idx = find(strcmp(probe_areas.acronym,'CP'));
        striatal_templates = any(cell2mat(arrayfun(@(x) ...
            isbetween(template_tipdist,probe_areas.tip_distance(x,2),probe_areas.tip_distance(x,1)) & ...
            template_shanks == probe_areas.probe_shank(x), ...
            striatum_probe_area_idx,'uni',false)'),2);

        dms_ccf_dv = 400;
        dms_templates = striatal_templates & (template_ccf(:,2) < dms_ccf_dv);

        % (don't use if only a few striatal units)
        if length(dms_templates) < 20
            continue
        end

        unit_ccf_all{curr_animal}{curr_day} = template_ccf(dms_templates,:);
        unit_psth_all{curr_animal}{curr_day} = unit_psth(dms_templates,:,unique(stim_x) == 90);
        responsive_units_all{curr_animal}{curr_day} = responsive_units(dms_templates);

    end
    ap.print_progress_fraction(curr_animal,length(animals));

end


% Smooth and normalize PSTH data
unit_psth_cat_smooth = smoothdata(cell2mat(cellfun(@(x) x(:,:,end),horzcat(unit_psth_all{:}),'uni',false)'),2,'gaussian',100);

softnorm = 1;
unit_psth_cat_norm = (unit_psth_cat_smooth-nanmean(unit_psth_cat_smooth(:,1:200,:),[2,3]))./ ...
    (nanmean(unit_psth_cat_smooth(:,1:200,:),[2,3])+softnorm);

responsive_units_cat = cell2mat(horzcat(responsive_units_all{:})');
responsive_units_ind = find(responsive_units_cat);

% % (sort amplitude)
% [~,sort_idx] = sort(nanmean(unit_psth_cat_norm(responsive_units_cat,200:400),2));
% (sort time)
[~,max_idx] = max(unit_psth_cat_norm(responsive_units_cat,200:400),[],2);
[~,sort_idx] = sort(max_idx);

figure;
imagesc(t_centers,[],unit_psth_cat_norm(responsive_units_ind(sort_idx),:));
clim([-10,10]);
title('Vis (V units)')
colormap(AP_colormap('BWR'));

% Plot responsive units in CCF
unit_ccf_cat = cell2mat(horzcat(unit_ccf_all{:})');

unit_color = padarray(+responsive_units_cat,[0,2],0,'pre');
unit_color(~any(unit_color,2),:) = 0.5;
unit_size = max(1,2*nanmean(unit_psth_cat_norm(:,200:400),2));

% Plot units on ccf (colored responsive, scaled by response, jittered x/y)
ap.ccf_outline_3d([],["brain","CP","GPe","RT"]);

xy_jitter = 20;
scatter3(unit_ccf_cat(:,1)+rand(size(unit_ccf_cat,1),1)*xy_jitter, ...
    unit_ccf_cat(:,3)+rand(size(unit_ccf_cat,1),1)*xy_jitter, ...
    unit_ccf_cat(:,2),unit_size,unit_color,'filled','MarkerFaceAlpha',0.5);


% Get fraction visually responsive units per animal
% x = cellfun(@(x) nanmean(vertcat(x{:})),responsive_units_all);
x = cellfun(@(x) mean(cellfun(@mean,x)),responsive_units_all);
figure;plot(x,'.k','MarkerSize',20);
ylabel('Frac vis. units');

% Get PSTH max
mua_psth = cell2mat(cellfun(@(x) nanmean(vertcat(x{:}),1),unit_psth_all,'uni',false));
mua_baseline = nanmean(mua_psth(:,isbetween(t_centers,-0.1,0)),2);
mua_psth_avg = (nanmean(mua_psth(:,isbetween(t_centers,0.05,0.15)),2)-mua_baseline)./mua_baseline;


%% widefield pre-ephys: task

animals = {'HA016','HA017','HA018','HA019','HA020','AP036','AP037'};

wf_kernel = cell(length(animals),1);
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    recordings = plab.find_recordings(animal);
    ephys_start_idx = find([recordings.ephys],1);

    % Use big task if present, small if not
    task_workflow = 'visual_operant_lick_two_stim_static_big_stim';
    task_recordings = plab.find_recordings(animal,[],task_workflow);
    if isempty(task_recordings)
        task_workflow = 'visual_operant_lick_two_stim_static';
        task_recordings = plab.find_recordings(animal,[],task_workflow);
    end

    use_days = find((datetime({task_recordings.day}) < datetime({recordings(ephys_start_idx).day})) & ...
        cellfun(@any,{task_recordings.widefield}));

    for curr_day = 1:length(use_days)

        % Set preload variables
        preload_vars = who;

        rec_day = task_recordings(use_days(curr_day)).day;
        rec_time = task_recordings(use_days(curr_day)).recording{end};

        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        n_trials = sum(cellfun(@(x) length(x) == 2,{trial_events.timestamps.StimOn}));

        % Set parameters for regression
        time_bins = [wf_t;wf_t(end)+1/wf_framerate];
        n_components = 200;
        frame_shifts = -5:40;
        lambda = 5;
        cv_fold = 5;

        skip_t = 60; % seconds start/end to skip for artifacts
        skip_frames = round(skip_t*wf_framerate);
     
        % Get stim on/move times  
        trial_stim_x = vertcat(trial_events.values(1:n_trials).TrialX);

        if contains(bonsai_workflow,'move')
            % Moving stim
            % (photodiode CS- = on/off for CS+ = on,pulse on move,off)
            stim_pd_n = (trial_stim_x == -90)*1 + (trial_stim_x == 90)*2;
            stim_pd_on_grouped = mat2cell(photodiode_on_times(1:sum(stim_pd_n)),stim_pd_n);
            stim_pd_off_grouped = mat2cell(photodiode_off_times(1:sum(stim_pd_n)),stim_pd_n);

            stimOn_times = cellfun(@(x) x(1), stim_pd_on_grouped);
            stim_move_times = cellfun(@(x) x(1), stim_pd_off_grouped);

        elseif contains(bonsai_workflow,'static')
            % Static stim
            stimOn_times = photodiode_on_times(1:n_trials);          
        end

        % Set regressors
        stim_regressors = zeros(0,length(time_bins)-1);
        stim_regressors(1,:) = histcounts(stimOn_times(trial_stim_x == 90),time_bins);
        stim_regressors(2,:) = histcounts(stimOn_times(trial_stim_x == -90),time_bins);
        if contains(bonsai_workflow,'move')
            stim_regressors(3,:) = histcounts(stim_move_times(trial_stim_x == 90),time_bins);
        end
      
        % Get stim kernel
        kernels = ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
            stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

        % Store
        wf_kernel{curr_animal}{curr_day} = kernels(:,:,1);

        % Clear load variables
        clearvars('-except',preload_vars{:});
        fprintf('%s %d/%d\n',animal,curr_day,length(use_days));
    end
end

n_components = 200;
wf_U = plab.wf.load_master_U(n_components);

data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\stim_lick\data';
data_filename = fullfile(data_path,'wf_kernel_task_ck');
save(data_filename,'wf_kernel');


px = plab.wf.svd2px(wf_U,cell2mat(permute(cellfun(@(x) ...
    mean(cat(3,x{:}),3),wf_kernel,'uni',false),[2,3,1])));

ap.imscroll(px);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG',[],1.5));
axis image


%% widefield pre-ephys: passive

animals = {'HA016','HA017','HA018','HA019','HA020','AP036','AP037'};

wf_kernel = cell(length(animals),1);
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    recordings = plab.find_recordings(animal);
    ephys_start_idx = find([recordings.ephys],1);

    % Use big task/passive if present, small if not
    % task_workflow = 'visual_operant_lick_two_stim_static_big_stim';
    % task_recordings = plab.find_recordings(animal,[],task_workflow);
    % if ~isempty(task_recordings)
    %     passive_workflow = 'lcr_passive_corner_CS\+_big_stim';
    % else
    %     task_workflow = 'visual_operant_lick_two_stim_static';
    %     task_recordings = plab.find_recordings(animal,[],task_workflow);
    %     passive_workflow = 'lcr_passive_corner_CS\+';
    % end
    % use_days = find((datetime({task_recordings.day}) < ...
    %     datetime({recordings(ephys_start_idx).day})) & ...
    %     cellfun(@any,{task_recordings.widefield}));

    % Use big passive
    task_workflow = 'visual*static*';
    task_recordings = plab.find_recordings(animal,[],task_workflow);

    passive_workflow = 'lcr_passive_corner_CS\+_big_stim';   
    passive_recordings = plab.find_recordings(animal,[],passive_workflow);
    use_days = {passive_recordings( ...
        datetime({passive_recordings.day}) >= datetime(task_recordings(1).day) & ...
        cellfun(@any,{passive_recordings.widefield})).day};

    for curr_day = 1:length(use_days)

        % Set preload variables
        preload_vars = who;

        rec_day = use_days{curr_day};        
        rec_time = plab.find_recordings(animal,rec_day,passive_workflow).recording{end};

        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % Set parameters for regression
        time_bins = [wf_t;wf_t(end)+1/wf_framerate];
        n_components = 200;
        frame_shifts = -5:40;
        lambda = 5;
        cv_fold = 5;

        skip_t = 60; % seconds start/end to skip for artifacts
        skip_frames = round(skip_t*wf_framerate);

        % Get stim kernels
        stim_x = vertcat(trial_events.values.TrialStimX);
        stim_x_unique = unique(stim_x);

        n_trials = min(length(stimOn_times),length(stim_x));

        stim_regressors = cell2mat(arrayfun(@(x) ...
            histcounts(stimOn_times(stim_x(1:n_trials) == x),time_bins), ...
            stim_x_unique,'uni',false));

        kernels = ...
            ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
            stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

        % Store (just 90deg stim)
        wf_kernel{curr_animal}{curr_day} = kernels(:,:,stim_x_unique == 90);

        % Clear load variables
        clearvars('-except',preload_vars{:});
        fprintf('%s %d/%d\n',animal,curr_day,length(use_days));
    end
end

n_components = 200;
wf_U = plab.wf.load_master_U(n_components);

data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\stim_lick\data';
data_filename = fullfile(data_path,'wf_kernel_passive_ck');
save(data_filename,'wf_kernel');

px = plab.wf.svd2px(wf_U,cell2mat(permute(cellfun(@(x) ...
    mean(cat(3,x{:}),3),wf_kernel,'uni',false),[2,3,1])));

ap.imscroll(px);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG',[],1.5));
axis image
set(gcf,'name','passive');


