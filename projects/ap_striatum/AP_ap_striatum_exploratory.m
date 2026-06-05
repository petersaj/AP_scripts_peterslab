%% A/P striatum
% some kind of iteration of comparing A/P striatum in naive, 2-stim
% learning, appetitive/aversive learning

%% DS AV posterior striatum: plot recording locations

% Animals that had posterior recordings + learned visual
animals = {'AP022','DS007','DS010','DS011', ... % V-A, performed both
    'DS005','DS000','DS013','DS014','DS015','DS016', ... % A-V, performed both
    'AP021','DS001', ... % V-A, didn't perform A mixed
    'DS003','DS004', ... % A-V, didn't perform A mixed
    'AP019', ... % V-A, unmixed tasks
    'AP020', ... % V-A, didn't learn A, unmixed tasks
    'AP018' ... % didn't learn A
    };

% not included:
% DS006 - didn't learn V = naive-ish?
% DS013 - didn't learn V, but started in mixed?

% Plot recording locations
animal_col = vertcat(ap.colormap('tube'),lines(7));
ccf_draw = ap.ccf_draw;
ccf_draw.draw_name('Caudoputamen');

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    probe_color = animal_col(curr_animal,:);
    ccf_draw.draw_probes_nte(animal,probe_color);
end

%% DS AV posterior striatum: plot recording locations

% Animals that had posterior recordings + learned visual
animals = {'AP022','DS007','DS010','DS011', ... % V-A, performed both
    'DS005','DS000','DS013','DS014','DS015','DS016', ... % A-V, performed both
    'AP021','DS001', ... % V-A, didn't perform A mixed
    'DS003','DS004', ... % A-V, didn't perform A mixed
    'AP019', ... % V-A, unmixed tasks
    'AP020', ... % V-A, didn't learn A, unmixed tasks
    'AP018' ... % didn't learn A
    };

% Find anterior/posterior striatum recordings (from NTE)
str_rec_days = cell(length(animals),2); % [ant,pos]
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};

    nte_filepattern = plab.locations.filename('server',animal,'*',[],'ephys','*probe_positions.mat');
    nte_fns = dir(nte_filepattern);

    % Get days for all NTE files
    day_pattern = digitsPattern(4) + '-' + digitsPattern(2) + '-' + digitsPattern(2);
    rec_days = extract({nte_fns.folder},day_pattern);

    % Load NTE files
    nte_all = cell(size(nte_fns));
    for curr_recording = 1:length(nte_fns)
        nte_all{curr_recording} = load(fullfile(nte_fns(curr_recording).folder, ...
            nte_fns(curr_recording).name));
    end

    % Get NTE files with A/P striatum
    ap_striatum_threshold = 550;
    anterior_recordings = cellfun(@(x) x.probe_positions_ccf{1}(1,2),nte_all) < ap_striatum_threshold;
    striatum_recordings = cellfun(@(x) any(strcmp(x.probe_areas{1}.name,'Caudoputamen')),nte_all);

    str_ant_rec = rec_days(striatum_recordings & anterior_recordings);
    str_pos_rec = rec_days(striatum_recordings & ~anterior_recordings);

    str_rec_days{curr_animal,1} = str_ant_rec;
    str_rec_days{curr_animal,2} = str_pos_rec;

end

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Get PSTHs and responsive units
unit_ccf_all = cell(length(animals),2);
unit_psth_all = cell(length(animals),2);
responsive_units_all = cell(length(animals),2);
for curr_animal = 1:length(animals)
    for curr_ap = 1:2
        for curr_day = 1:length(str_rec_days{curr_animal,curr_ap})
            for curr_modality = 1:2

                % Load data
                animal = animals{curr_animal};
                rec_day = str_rec_days{curr_animal,curr_ap}{curr_day};

                switch curr_modality
                    case 1
                        workflow = 'lcr_passive';
                    case 2
                        workflow = 'hml_passive_audio';
                end
                recordings = plab.find_recordings(animal,rec_day,workflow);
                rec_time = recordings.recording{end};

                ap.load_recording;

                % Get quiescent trials
                stim_window = [0,0.5];
                quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
                    timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
                    timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
                    (1:length(stimOn_times))');

                % Get stim times to use
                if isfield(trial_events.values,'TrialStimX')
                    align_category = vertcat(trial_events.values.TrialStimX);
                elseif isfield(trial_events.values,'StimFrequence')
                    align_category = vertcat(trial_events.values.StimFrequence);
                end
                align_times = cellfun(@(x) ...
                    stimOn_times(align_category(1:length(stimOn_times)) == x & ...
                    quiescent_trials),num2cell(unique(align_category)),'uni',false);

                % Get PSTH by 2D histogram
                n_units = size(templates,1);
                unit_psth = nan(n_units,length(t_bins)-1,length(align_times));
                for curr_align = 1:length(align_times)
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
                switch curr_modality
                    case 1
                        response_align = align_times{3}; % response = R vis
                    case 2
                        response_align = align_times{2}; % response = M aud
                end

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

                % Get striatum templates
                striatum_idx = strcmp(probe_areas{1}.name,'Caudoputamen');
                striatum_depth = prctile(probe_areas{1}.tip_distance(striatum_idx,:),[0,100],'all');
                striatum_templates = isbetween(template_tipdist,striatum_depth(1),striatum_depth(2));

                % Store unit PSTHs and responsive units (striatum units only)
                unit_ccf_all{curr_animal,curr_ap}{curr_day,curr_modality} = template_ccf(striatum_templates,:,:);
                unit_psth_all{curr_animal,curr_ap}{curr_day,curr_modality} = unit_psth(striatum_templates,:,:);
                responsive_units_all{curr_animal,curr_ap}{curr_day,curr_modality} = responsive_units(striatum_templates);

            end
        end
    end
    ap.print_progress_fraction(curr_animal,length(animals));
end

% Smooth and normalize PSTH data
unit_psth_all_smooth = cellfun(@(x) cellfun(@(x) ...
    smoothdata(x,2,'gaussian',100),x,'uni',false), ...
    unit_psth_all,'uni',false);

softnorm = 1;
unit_psth_all_norm = cellfun(@(x) cellfun(@(x) ...
    (x-nanmean(x(:,1:500,:,:),[2,3]))./(nanmean(x(:,1:500,:,:),[2,3])+softnorm),x,'uni',false), ...
    unit_psth_all_smooth,'uni',false);


% Plot PSTHs for V/A units during V/A stim
str_ap = 1;
use_animals = true(size(animals));

vis_psth_cat = cell2mat(cellfun(@(x) vertcat(x{:,1}),unit_psth_all_norm(use_animals,str_ap),'uni',false));
aud_psth_cat = cell2mat(cellfun(@(x) vertcat(x{:,2}),unit_psth_all_norm(use_animals,str_ap),'uni',false));

v_units_cat = find(cell2mat(cellfun(@(x) vertcat(x{:,1}),responsive_units_all(use_animals,str_ap),'uni',false)));
[~,v_sort_idx] = sort(nanmean(vis_psth_cat(v_units_cat,500:700,3),2));

a_units_cat = find(cell2mat(cellfun(@(x) vertcat(x{:,2}),responsive_units_all(use_animals,str_ap),'uni',false)));
[~,a_sort_idx] = sort(nanmean(aud_psth_cat(a_units_cat,500:700,2),2));

figure; colormap(AP_colormap('BWR'));
h = tiledlayout(2,2);

nexttile;
imagesc(vis_psth_cat(v_units_cat(v_sort_idx),:,3));
clim([-10,10]);
title('Vis (V units)')

nexttile;
imagesc(aud_psth_cat(v_units_cat(v_sort_idx),:,2));
clim([-10,10]);
title('Aud (V units)')

nexttile;
imagesc(vis_psth_cat(a_units_cat(a_sort_idx),:,3));
clim([-10,10]);
title('Vis (A units)')

nexttile;
imagesc(aud_psth_cat(a_units_cat(a_sort_idx),:,2));
clim([-10,10]);
title('Aud (A units)')


% Plot PSTHs for V/A units during V/A stim
use_animals = true(size(animals));

vis_psth_cat = cell2mat(cellfun(@(x) vertcat(x{:,1}),unit_psth_all_norm(use_animals,str_ap),'uni',false));
aud_psth_cat = cell2mat(cellfun(@(x) vertcat(x{:,2}),unit_psth_all_norm(use_animals,str_ap),'uni',false));

responsive_units_cat = find(cell2mat(cellfun(@(x) vertcat(x{:,1}),responsive_units_all(use_animals,str_ap),'uni',false)));
% % (sort amplitude)
% [~,sort_idx] = sort(nanmean(vis_psth_cat(responsive_units_cat,500:700,3),2));
% (sort timing)
[~,max_t] = max(vis_psth_cat(responsive_units_cat,:,3),[],2);
[~,sort_idx] = sort(max_t);


figure;
colormap(AP_colormap('BWR'));
h = tiledlayout(1,2);

nexttile;
imagesc(vis_psth_cat(responsive_units_cat(sort_idx),:,3));
clim([-10,10]);
title('Vis (V units)')

nexttile;
imagesc(aud_psth_cat(responsive_units_cat(sort_idx),:,2));
clim([-10,10]);
title('Aud (V units)')


% Plot responsive units in CCF
unit_ccf_cat = cell2mat(cellfun(@(x) vertcat(x{:,1}),unit_ccf_all(use_animals,str_ap),'uni',false));
unit_color = [ismember((1:size(unit_ccf_cat,1))',a_units_cat), ...
    zeros(size(unit_ccf_cat,1),1), ...
    ismember((1:size(unit_ccf_cat,1))',v_units_cat)];
unit_color(~any(unit_color,2),:) = 0.5;

unit_size = max(1,10*nanmean(vis_psth_cat(:,500:700,3),2));

% Plot units on ccf (colored responsive, scaled by response, jittered x/y)
ap.ccf_outline_3d([],[],["brain","caudoputamen"]);

xy_jitter = 20;
scatter3(unit_ccf_cat(:,1)+rand(size(unit_ccf_cat,1),1)*xy_jitter, ...
    unit_ccf_cat(:,3)+rand(size(unit_ccf_cat,1),1)*xy_jitter, ...
    unit_ccf_cat(:,2),unit_size,unit_color,'filled','MarkerFaceAlpha',0.5);


%% Spatial with PG data (prelim test)

% Animals that had posterior recordings + learned visual
animals = {'PG001','PG003','PG006','PG007','PG008','PG009'};

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

    % Load data
    animal = animals{curr_animal};
    % workflow = 'lcr_passive_corner*';
    workflow = 'lcr_passive_size60';
    recordings = plab.find_recordings(animal,[],workflow);
    recordings = recordings([recordings.ephys]>0);

    if isempty(recordings)
        continue
    end

    for curr_day = 1:length(recordings)

        rec_day = recordings(curr_day).day;
        rec_time = recordings(curr_day).recording{end};
        ap.load_recording;

        % Get quiescent trials
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');

        % Get stim times to use
        align_category = vertcat(trial_events.values.TrialStimX);

        n_trials = min(length(stimOn_times),length(align_category));
        align_times = cellfun(@(x) ...
            stimOn_times(align_category(1:n_trials) == x & ...
            quiescent_trials(1:n_trials)),num2cell(unique(align_category)),'uni',false);

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
        response_align = align_times{unique(align_category) == 90}; % response = R vis

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
        unit_ccf_all{curr_animal}{curr_day} = template_ccf;
        unit_psth_all{curr_animal}{curr_day} = unit_psth;
        responsive_units_all{curr_animal}{curr_day} = responsive_units;

    end
    ap.print_progress_fraction(curr_animal,length(animals));

end


% Smooth and normalize PSTH data
unit_psth_cat_smooth = smoothdata(cell2mat(horzcat(unit_psth_all{:})'),2,'gaussian',100);

softnorm = 1;
unit_psth_cat_norm = unit_psth_cat_smooth-nanmean(unit_psth_cat_smooth(:,1:200,:),[2,3])./ ...
    (nanmean(unit_psth_cat_smooth(:,1:200,:),[2,3])+softnorm);

responsive_units_cat = cell2mat(horzcat(responsive_units_all{:})');
responsive_units_ind = find(responsive_units_cat);

[~,sort_idx] = sort(nanmean(unit_psth_cat_norm(responsive_units_cat,200:400,3),2));

figure;
imagesc(unit_psth_cat_norm(responsive_units_ind(sort_idx),:,3));
clim([-10,10]);
title('Vis (V units)')
colormap(AP_colormap('BWR'));

% Plot responsive units in CCF
unit_ccf_cat = cell2mat(horzcat(unit_ccf_all{:})');

unit_color = padarray(+responsive_units_cat,[0,2],0,'pre');
unit_color(~any(unit_color,2),:) = 0.5;
unit_size = max(1,2*nanmean(unit_psth_cat_norm(:,200:400,3),2));

% Plot units on ccf (colored responsive, scaled by response, jittered x/y)
ap.ccf_outline_3d([],[],["brain","CP","GPe"]);

xy_jitter = 20;
scatter3(unit_ccf_cat(:,1)+rand(size(unit_ccf_cat,1),1)*xy_jitter, ...
    unit_ccf_cat(:,3)+rand(size(unit_ccf_cat,1),1)*xy_jitter, ...
    unit_ccf_cat(:,2),unit_size,unit_color,'filled','MarkerFaceAlpha',0.5);


