%% Da's secondary-ish project with multiple visual stimuli

%% TESTING: all animals that have learned the stim task

% Find all animals from task
% task_dir = dir(fullfile(plab.locations.server_data_path, ...
%     '**','bonsai','**','stim_wheel_right_stage2.bonsai'));

task_dir_topfolders = unique(extractBetween({task_dir.folder}, ...
    [plab.locations.server_data_path,filesep], ...
    filesep));

animals = task_dir_topfolders(matches(task_dir_topfolders,lettersPattern(2)+digitsPattern(3)));

% Get learned days
learned_days = struct('animal',animals,'days',cell(size(animals)),'learned',cell(size(animals)));
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    task_recordings = plab.find_recordings(animal,[],{'stim_wheel_right_stage1','stim_wheel_right_stage2'});

    for curr_recording = 1:length(task_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = task_recordings(curr_recording).day;
        rec_time = task_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % (only use if there were at least 10 trials)
        if n_trials < 10
            continue
        end

        % Get reaction stat
        use_stat = 'mean';
        rxn_stat_p = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);
        learned_day = rxn_stat_p < 0.05;

        % Package data
        learned_days(curr_animal).days{curr_recording} = rec_day;
        learned_days(curr_animal).learned(curr_recording) = learned_day;

        ap.print_progress_fraction(curr_recording,length(task_recordings));
    end
end

% Save
save_path = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data';
save_fn = fullfile(save_path,'position_learned_days');
save(save_fn,'learned_days');

%% ~~~~~~~~~~~~~ All animals, all save

%% Grab behavior info for all animals

animals = { ...
    'AP019','AP021','AP022','DS001','DS007','DS010','DS011', ... % V > A
    'DS000','DS003','DS004','DS014','DS015','DS016','DS005', ... % A > V
    'DS019','DS020','DS021','AP027','AP028','AP029', ... % V-alt > V > Afreq
    'AP020','AP018', ... % V > A-nonlearn
    'DS006','DS013', ... % A > V-nonlearn
    'HA003','HA004', ... % V-alt > V
    'HA000','HA001','HA002' ... % V-alt > V
    };

bhv_fieldnames = ...
    {'animal','day','workflow','rxn_stat_p','rxn_stat','rxn_stat_null'};
bhv = cell2struct(cell(length(bhv_fieldnames),0),bhv_fieldnames);

for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    % Find learned days
    task_recordings = plab.find_recordings(animal,[],'stim_wheel*');
    use_recordings = task_recordings(cellfun(@any,{task_recordings.widefield}));

    bhv(curr_animal).animal = animal;
    bhv(curr_animal).day = cell(size(use_recordings));
    bhv(curr_animal).workflow = cell(size(use_recordings));
    bhv(curr_animal).rxn_stat_p = nan(size(use_recordings));
    bhv(curr_animal).rxn_stat = nan(size(use_recordings));
    bhv(curr_animal).rxn_stat_null = nan(size(use_recordings));

    for curr_recording = 1:length(use_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = use_recordings(curr_recording).day;
        rec_time = use_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        if n_trials < 10
            continue
        end

        % Store workflow info
        bhv(curr_animal).day{curr_recording} = rec_day;
        bhv(curr_animal).workflow{curr_recording} = bonsai_workflow;

        % Get reaction stat
        use_stat = 'mean';
        [bhv(curr_animal).rxn_stat_p(curr_recording), ...
            bhv(curr_animal).rxn_stat(curr_recording), ...
            bhv(curr_animal).rxn_stat_null(curr_recording)] = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

        ap.print_progress_fraction(curr_recording,length(use_recordings));

    end

end

% Save
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data';
save_fn = fullfile(data_path,'bhv');
save(save_fn,'bhv')
fprintf('Saved %s\n',save_fn);

%% (testing: grab matching days from bhv)

bhv_filename = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\bhv.mat';
load(bhv_filename);

task_workflow = 'stim_wheel_right_stage\d$';
task_workflow = 'stim_wheel_right_stage\d_audio*';

workflow_match = cellfun(@(bhv) cellfun(@(x) ~isempty(regexp(x,task_workflow)), ...
    bhv,'ErrorHandler',@(varargin) false),{bhv.workflow},'uni',false);
workflow_learned = cellfun(@(x) x < 0.05,{bhv.rxn_stat_p},'uni',false);

use_days = cellfun(@(x,y) find(x(y)),workflow_match,workflow_learned,'uni',false);


%% Stim average all animals

bhv_filename = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\bhv.mat';
load(bhv_filename);

% Set workflows
vis_workflow = 'lcr_passive';
aud_workflow = 'hml_passive_audio';
task_workflow = 'stim_wheel*';

% Set times for PSTH        
raster_window = [-0.5,1];
psth_bin_size = 0.03;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

stim_V = cell(length(bhv),1);

for curr_animal = 18:length(bhv)
   
    animal = bhv(curr_animal).animal;
    fprintf([animal,'\n']);

    rec_days = bhv(curr_animal).day;

    for curr_day = find(~cellfun(@isempty,rec_days))
        for curr_modality = 1:2
            switch curr_modality
                case 1
                    curr_workflow = vis_workflow;
                case 2
                    curr_workflow = aud_workflow;
            end

            % Grab pre-load vars
            preload_vars = who;

            % Load data
            curr_recording = plab.find_recordings(animal,rec_days{curr_day},curr_workflow);
            if isempty(curr_recording)
                continue
            end

            rec_day = curr_recording.day;
            rec_time = curr_recording.recording{end};

            load_parts = struct;
            load_parts.widefield = true;
            load_parts.widefield_master = true;
            
            % (Try/catch: rare catastrophic dropped frames)
            try
                ap.load_recording;
            catch me
                continue
            end

            % If widefield isn't aligned, skip
            if ~load_parts.widefield_align
                continue
            end

            % Stim times (quiescent trials only)
            switch curr_modality
                case 1
                    stim_type = vertcat(trial_events.values.TrialStimX);             
                case 2
                    stim_type = vertcat(trial_events.values.StimFrequence);
            end
            
            stim_window = [0,0.5];
            quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
                timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
                timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
                1:length(stimOn_times))';
      
            % sometimes not all trials shown?
            stim_type = stim_type(1:length(stimOn_times));
            align_times = cellfun(@(x) stimOn_times(stim_type == x & quiescent_trials), ...
                num2cell(unique(stim_type)),'uni',false);

            % Align V to align times
            aligned_V_mean = nan(size(wf_V,1),length(t_centers),length(align_times));
            for curr_align = 1:length(align_times)

                % (skip if <5 usable trials)
                if length(align_times{curr_align}) < 5
                    continue
                end

                peri_event_t = align_times{curr_align} + t_centers;

                aligned_V = interp1(wf_t,wf_V',peri_event_t);

                aligned_V_baselinesub = aligned_V - ...
                    mean(aligned_V(:,t_centers < 0,:),2);

                aligned_V_mean(:,:,curr_align) = ...
                    permute(nanmean(aligned_V_baselinesub,1),[3,2,1]);
            end

            % Store mean aligned ROI trace
            stim_V{curr_animal}{curr_day,curr_modality} = aligned_V_mean;

            % Prep for next loop
            clearvars('-except',preload_vars{:});

        end
        ap.print_progress_fraction(curr_day,length(rec_days));
    end
    disp(['Done ' animal]);
end

% Save
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data';
save_fn = fullfile(data_path,'passive_avg');
save(save_fn,'stim_V','-v7.3')
fprintf('Saved %s\n',save_fn);

%% Stim kernel all animals

bhv_filename = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\bhv.mat';
load(bhv_filename);

% Set workflows
vis_workflow = 'lcr_passive';
aud_workflow = 'hml_passive_audio';
task_workflow = 'stim_wheel*';

stim_V = cell(length(bhv),1);

for curr_animal = 1:length(bhv)
   
    animal = bhv(curr_animal).animal;
    fprintf([animal,'\n']);

    rec_days = bhv(curr_animal).day;

    for curr_day = find(~cellfun(@isempty,rec_days))
        for curr_modality = 1:2
            switch curr_modality
                case 1
                    curr_workflow = vis_workflow;
                case 2
                    curr_workflow = aud_workflow;
            end

            % Grab pre-load vars
            preload_vars = who;

            % Load data
            curr_recording = plab.find_recordings(animal,rec_days{curr_day},curr_workflow);
            if isempty(curr_recording)
                continue
            end

            rec_day = curr_recording.day;
            rec_time = curr_recording.recording{end};

            load_parts = struct;
            load_parts.widefield = true;
            load_parts.widefield_master = true;
            
            % (Try/catch: rare catastrophic dropped frames)
            try
                ap.load_recording;
            catch me
                continue
            end

            % If widefield isn't aligned, skip
            if ~load_parts.widefield_align
                continue
            end

            % Stim times (quiescent trials only)
            switch curr_modality
                case 1
                    stim_type = vertcat(trial_events.values.TrialStimX);             
                case 2
                    stim_type = vertcat(trial_events.values.StimFrequence);
            end

            time_bins = [wf_t;wf_t(end)+1/wf_framerate];

            stim_regressors = cell2mat(arrayfun(@(x) ...
                histcounts(stimOn_times(stim_type == x),time_bins), ...
                unique(stim_type),'uni',false));

            n_components = 500;

            frame_shifts = 0:30;
            lambda = 20;

            skip_t = 3; % seconds start/end to skip for artifacts
            skip_frames = round(skip_t*wf_framerate);
            [kernels,predicted_signals,explained_var] = ...
                ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
                stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda);

            % Store mean aligned ROI trace
            stim_V{curr_animal}{curr_day,curr_modality} = kernels;

            % Prep for next loop
            clearvars('-except',preload_vars{:});

        end
        ap.print_progress_fraction(curr_day,length(rec_days));
    end
    disp(['Done ' animal]);
end

% Save
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data';
save_fn = fullfile(data_path,'passive_kernel');
save(save_fn,'stim_V','-v7.3')
fprintf('Saved %s\n',save_fn);

%% Task kernel all animals

bhv_filename = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\bhv.mat';
load(bhv_filename);

% Set workflow
task_workflow = 'stim_wheel*';

stim_V = cell(length(bhv),1);

for curr_animal = 1:length(bhv)

    animal = bhv(curr_animal).animal;
    fprintf([animal,'\n']);

    rec_days = bhv(curr_animal).day;

    for curr_day = find(~cellfun(@isempty,rec_days))

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        curr_recording = plab.find_recordings(animal,rec_days{curr_day},task_workflow);
        if isempty(curr_recording)
            continue
        end

        rec_day = curr_recording.day;
        rec_time = curr_recording.recording{end};

        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;

        % (Try/catch: rare catastrophic dropped frames)
        try
            ap.load_recording;
        catch me
            continue
        end

        % If widefield isn't aligned, skip
        if ~load_parts.widefield_align
            continue
        end

        time_bins = [wf_t;wf_t(end)+1/wf_framerate];

        stim_regressors = histcounts(stimOn_times,time_bins);

        n_components = 200;

        frame_shifts = 0:20;
        lambda = 20;
        cv_fold = 3;

        skip_t = 60; % seconds start/end to skip for artifacts
        skip_frames = round(skip_t*wf_framerate);
        [kernels,predicted_signals,explained_var] = ...
            ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
            stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

        stim_V{curr_animal}{curr_day} = kernels;

        % Prep for next loop
        clearvars('-except',preload_vars{:});

        ap.print_progress_fraction(curr_day,length(rec_days));
    end
    disp(['Done ' animal]);
end

% Save
data_path = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data';
save_fn = fullfile(data_path,'task_kernel');
save(save_fn,'stim_V','-v7.3')
fprintf('Saved %s\n',save_fn);


%% >>>>>> Batch analysis

%% Behavior

bhv_filename = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\bhv.mat';
load(bhv_filename);

% Concatenate 
rxn_stat_cat = horzcat(bhv.rxn_stat)';
rxn_stat_null_cat = horzcat(bhv.rxn_stat_null)';
rxn_stat_idx_cat = (rxn_stat_null_cat - rxn_stat_cat)./(rxn_stat_null_cat + rxn_stat_cat);

rxn_stat_p_cat = horzcat(bhv.rxn_stat_p)';

% Get modality order
workflow_animal = {bhv.workflow}';
vis_workflow = 'stim_wheel_right_stage\d$';
aud_workflow = 'stim_wheel_right_stage\d_audio*';
va_animals = cellfun(@(x) max([...
    find(cellfun(@any,regexp(x,vis_workflow)),1) < ...
    find(cellfun(@any,regexp(x,aud_workflow)),1), ...
    false]),workflow_animal);
av_animals = cellfun(@(x) max([...
    find(cellfun(@any,regexp(x,vis_workflow)),1) > ...
    find(cellfun(@any,regexp(x,aud_workflow)),1), ...
    false]),workflow_animal);

% (ALT: define manually)
va_animals = ismember({bhv.animal},{'AP019','AP021','AP022','DS001','DS007','DS010','DS011'})';
av_animals = ismember({bhv.animal},{'DS000','DS004','DS014','DS015','DS016','DS005'})';

workflow_cat = horzcat(bhv.workflow)';
va_animals_cat = cell2mat(cellfun(@(x,y) repmat(x,length(y),1),num2cell(va_animals),workflow_animal,'uni',false));
av_animals_cat = cell2mat(cellfun(@(x,y) repmat(x,length(y),1),num2cell(av_animals),workflow_animal,'uni',false));

% Get days from V/A task start
v_day = cell2mat(cellfun(@(x) ...
    (cumsum(cellfun(@any,regexp(x,vis_workflow))).*AP_nanout(~cellfun(@any,regexp(x,vis_workflow))))', ...
    workflow_animal,'uni',false));

a_day = cell2mat(cellfun(@(x) ...
    (cumsum(cellfun(@any,regexp(x,aud_workflow))).*AP_nanout(~cellfun(@any,regexp(x,aud_workflow))))', ...
    workflow_animal,'uni',false));

% Get and plot number of days to learning
v_learned_day = cellfun(@(x,p,stat) ...
     sum(cellfun(@any,regexp( ...
     x(1:find(cellfun(@any,regexp(x,vis_workflow)) & p < 0.05 & stat < 1,1)),vis_workflow))), ...
     workflow_animal,{bhv.rxn_stat_p}',{bhv.rxn_stat}');
 
a_learned_day = cellfun(@(x,p,stat) ...
     sum(cellfun(@any,regexp( ...
     x(1:find(cellfun(@any,regexp(x,aud_workflow)) & p < 0.05 & stat < 1,1)),aud_workflow))), ...
     workflow_animal,{bhv.rxn_stat_p}',{bhv.rxn_stat}');

% Plot days to learn modality 1/2 by group
order_grp = sum([va_animals,av_animals].*[1:2],2).*AP_nanout(sum([va_animals,av_animals],2) == 0);
va_ld_avg = [ap.groupfun(@mean,v_learned_day,order_grp), ...
    ap.groupfun(@median,a_learned_day,order_grp)];
va_ld_sem = [ap.groupfun(@AP_sem,v_learned_day,order_grp), ...
    ap.groupfun(@AP_sem,a_learned_day,order_grp)];

% (plot modality by group)
figure; h = tiledlayout(1,2);
nexttile;
cats = categorical({'Vis','Aud'});
cats = reordercats(cats,cellstr(cats));
errorbar(cats,va_ld_avg',va_ld_sem','linewidth',2);
legend({'VA','AV'});
ylabel('Med. days to learn');

% (plot modality order by group)
nexttile;
cats = categorical({'Modality 1','Modality 2'});
cats = reordercats(cats,cellstr(cats));
errorbar(cats,[va_ld_avg(1,:)',fliplr(va_ld_avg(2,:))'], ...
    [va_ld_sem(1,:)',fliplr(va_ld_sem(2,:))'],'linewidth',2);
legend({'VA','AV'});
ylabel('Mean days to learn');

% Plot reaction time stat
figure; h = tiledlayout(1,2);
plot_bhv = rxn_stat_idx_cat;

nexttile; hold on; title('Visual learning');
[x,g] = ap.groupfun(@nanmean,plot_bhv(va_animals_cat),v_day(va_animals_cat));
e = ap.groupfun(@AP_sem,plot_bhv(va_animals_cat),v_day(va_animals_cat));
ap.errorfill(g(g<=7),x(g<=7),e(g<=7),'r');

[x,g] = ap.groupfun(@nanmean,plot_bhv(av_animals_cat),v_day(av_animals_cat));
e = ap.groupfun(@AP_sem,plot_bhv(av_animals_cat),v_day(av_animals_cat));
ap.errorfill(g(g<=7),x(g<=7),e(g<=7),'b');

nexttile; hold on; title('Auditory learning');
[x,g] = ap.groupfun(@nanmean,plot_bhv(va_animals_cat),a_day(va_animals_cat));
e = ap.groupfun(@AP_sem,plot_bhv(va_animals_cat),a_day(va_animals_cat));
ap.errorfill(g(g<=7),x(g<=7),e(g<=7),'r');

[x,g] = ap.groupfun(@nanmean,plot_bhv(av_animals_cat),a_day(av_animals_cat));
e = ap.groupfun(@AP_sem,plot_bhv(av_animals_cat),a_day(av_animals_cat));
ap.errorfill(g(g<=7),x(g<=7),e(g<=7),'b');
legend({'VA','AV'})
linkaxes(h.Children);

%% Passive maps

% Kernel
data_fn = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\passive_kernel.mat';
load(data_fn);
t = (0:30)/35;
n_components = size(stim_V{1}{1},1);

% % Average
% data_fn = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\passive_avg.mat';
% load(data_fn);
% t = conv(-0.5:0.03:1,[1,1]/2,'valid');
% n_components = size(stim_V{1}{1},1);

bhv_filename = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\bhv.mat';
load(bhv_filename);

U_master = plab.wf.load_master_U;

%%% PROBLEM! bhv doesn't match kernels (sometime extra bhv day), don't know
%%% why at the moment, maybe just truncated if error on last day, but
%%% didn't see error on example one

stim_V_cat = vertcat(stim_V{:});

% just truncate bhv values for now
rxn_stat_p_animal = {bhv.rxn_stat_p}';
rxn_stat_p_animal_trunc = cellfun(@(x,y) x(1:size(y,1))',rxn_stat_p_animal,stim_V,'uni',false);
rxn_stat_p_cat = vertcat(rxn_stat_p_animal_trunc{:});

workflow_animal = {bhv.workflow}';
workflow_animal_trunc = cellfun(@(x,y) x(1:size(y,1))',workflow_animal,stim_V,'uni',false);
workflow_cat = vertcat(workflow_animal_trunc{:});

vis_workflow = 'stim_wheel_right_stage\d$';
aud_workflow = 'stim_wheel_right_stage\d_audio*';
va_animals = cellfun(@(x) max([...
    find(cellfun(@any,regexp(x,vis_workflow)),1) < ...
    find(cellfun(@any,regexp(x,aud_workflow)),1), ...
    false]),workflow_animal);
av_animals = cellfun(@(x) max([...
    find(cellfun(@any,regexp(x,vis_workflow)),1) > ...
    find(cellfun(@any,regexp(x,aud_workflow)),1), ...
    false]),workflow_animal);

% (ALT: define manually)
% va_animals = ismember({bhv.animal},{'AP019','AP021','AP022','DS001','DS007','DS010','DS011'})';
% va_animals = ismember({bhv.animal},{'AP019','AP021','AP022','DS001','DS007','DS010','DS011','DS019','DS020','DS021'})';
va_animals = ismember({bhv.animal},{'AP019','AP021','AP022','DS001','DS007','DS010','DS011','DS019','DS020','DS021','AP027','AP029','AP029'})';
av_animals = ismember({bhv.animal},{'DS000','DS004','DS014','DS015','DS016','DS005'})';

va_animals_cat = cell2mat(cellfun(@(x,y) repmat(x,length(y),1),num2cell(va_animals),workflow_animal_trunc,'uni',false));
av_animals_cat = cell2mat(cellfun(@(x,y) repmat(x,length(y),1),num2cell(av_animals),workflow_animal_trunc,'uni',false));

task_workflow = 'stim_wheel_right_stage\d$';
use_rec = av_animals_cat & cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat < 0.05;
v = plab.wf.svd2px(U_master(:,:,1:n_components),nanmean(cat(4,stim_V_cat{use_rec,1}),4));

task_workflow = 'stim_wheel_right_stage\d_audio*';
use_rec = av_animals_cat & cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat < 0.05;
a = plab.wf.svd2px(U_master(:,:,1:n_components),nanmean(cat(4,stim_V_cat{use_rec,2}),4));

ap.imscroll([v,a],t);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));

% Plot max for all stim
stim_t = t > 0 & t < 0.2;
v_max = squeeze(max(v(:,:,stim_t,:),[],3));
a_max = squeeze(max(a(:,:,stim_t,:),[],3));

figure;
h = tiledlayout(2,3);
for i = 1:3
    nexttile;
    imagesc(v_max(:,:,i));
    axis image off;
    clim(col_lim)
    colormap(gca,ap.colormap('WR'));
    ap.wf_draw('ccf','k');
end
for i = 1:3
    nexttile;
    imagesc(a_max(:,:,i));
    axis image off;
    clim(col_lim)
    colormap(gca,ap.colormap('WB'));
    ap.wf_draw('ccf','k');
end

% Colored learned map
v_stim = 3;
a_stim = 2;

task_workflow = 'stim_wheel_right_stage\d$';
use_rec = cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat < 0.05;
% use_rec = va_animals_cat & cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat > 0.05;
v_learned = plab.wf.svd2px(U_master(:,:,1:n_components),nanmean(cat(4,stim_V_cat{use_rec,1}),4));
v_learned = imgaussfilt(v_learned,3)-imgaussfilt(v_learned,20);

task_workflow = 'stim_wheel_right_stage\d_audio*';
use_rec =  av_animals_cat & cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat < 0.05;
% use_rec = av_animals_cat & cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat > 0.05;
a_learned = plab.wf.svd2px(U_master(:,:,1:n_components),nanmean(cat(4,stim_V_cat{use_rec,2}),4));
a_learned = imgaussfilt(a_learned,3)-imgaussfilt(a_learned,20);

v_col = ap.colormap('WR');
a_col = ap.colormap('WB');
col_lim = [0,1.4e-4]; % kernel
% col_lim = [0,6e-3]; % avg

v_learned_gray = 1+round(mat2gray(v_learned(:,:,:,v_stim),col_lim).*(size(v_col,1)-1));
v_learned_col = permute(reshape(v_col(v_learned_gray,:),[size(v_learned_gray),3]),[1,2,4,3]);

a_learned_gray = 1+round(mat2gray(a_learned(:,:,:,a_stim),col_lim).*(size(a_col,1)-1));
a_learned_col = permute(reshape(a_col(a_learned_gray,:),[size(a_learned_gray),3]),[1,2,4,3]);

va_learned_col = min(v_learned_col,a_learned_col);

ap.imscroll(va_learned_col,[],true);
axis image;
ap.wf_draw('ccf','k');

plot_frames = find(t>=0 & t<= 0.2);
figure; h = tiledlayout(3,length(plot_frames),'TileSpacing','none','TileIndexing','ColumnMajor');
for curr_frame = plot_frames
    nexttile;
    image(v_learned_col(:,:,:,curr_frame));
    ap.wf_draw('ccf','k');
    axis image off;
    title(sprintf('%.2f',(curr_frame-1)/35));

    nexttile;
    image(a_learned_col(:,:,:,curr_frame));
    ap.wf_draw('ccf','k');
    axis image off;

    nexttile;
    image(va_learned_col(:,:,:,curr_frame));
    ap.wf_draw('ccf','k');
    axis image off;
end

% Plot max V/A/overlay
stim_t = t > 0 & t < 0.2;
v_learned_max = max(v_learned(:,:,stim_t,v_stim),[],3);
a_learned_max = max(a_learned(:,:,stim_t,a_stim),[],3);

figure;
h = tiledlayout(1,3,'TileSpacing','none');

nexttile; imagesc(v_learned_max); 
clim(col_lim); axis image off;
colormap(gca,ap.colormap('WR'));
ap.wf_draw('ccf','k');

nexttile; imagesc(a_learned_max); 
clim(col_lim); axis image off;
colormap(gca,ap.colormap('WB'));
ap.wf_draw('ccf','k');

v_learned_max_gray = 1+round(mat2gray(v_learned_max,col_lim).*(size(v_col,1)-1));
v_learned_max_col = reshape(v_col(v_learned_max_gray,:),size(v_learned_max,1),size(v_learned_max,2),3,[]);

a_learned_max_gray = 1+round(mat2gray(a_learned_max,col_lim).*(size(a_col,1)-1));
a_learned_max_col = reshape(a_col(a_learned_max_gray,:),size(a_learned_max,1),size(a_learned_max,2),3,[]);

va_learned_max_col = min(v_learned_max_col,a_learned_max_col);
nexttile; image(va_learned_max_col);
axis image off;
ap.wf_draw('ccf','k');





%% ROI

% Passive kernel
data_fn = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\passive_kernel.mat';
load(data_fn);
t = (0:30)/35;
n_components = size(stim_V{1}{1});

% % Passive average
% data_fn = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\passive_avg.mat';
% load(data_fn);
% t = conv(-0.5:0.03:1,[1,1]/2,'valid');
% n_components = size(stim_V{1}{1});

% % Task kernel
% data_fn = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\task_kernel.mat';
% load(data_fn);
% t = (0:20)/35;
% n_components = size(stim_V{1}{1});

bhv_filename = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\bhv.mat';
load(bhv_filename);

U_master = plab.wf.load_master_U;

%%% PROBLEM! bhv doesn't match kernels (sometime extra bhv day), don't know
%%% why at the moment, maybe just truncated if error on last day, but
%%% didn't see error on example one

% Task/passive params
if contains(data_fn,'task')
    stim_V_cat = horzcat(stim_V{:})';
    [v_stim,v_col,a_stim,a_col] = deal(1);
elseif contains(data_fn,'passive')
    stim_V_cat = vertcat(stim_V{:});
    v_stim = 3;
    v_col = 1;

    a_stim = 2;
    a_col = 2;
end

% just truncate bhv values for now
rxn_stat_p_animal = {bhv.rxn_stat_p}';
rxn_stat_p_animal_trunc = cellfun(@(x,y) x(1:length(y))',rxn_stat_p_animal,stim_V,'uni',false);
rxn_stat_p_cat = vertcat(rxn_stat_p_animal_trunc{:});

workflow_animal = {bhv.workflow}';
workflow_animal_trunc = cellfun(@(x,y) x(1:length(y))',workflow_animal,stim_V,'uni',false);
workflow_cat = vertcat(workflow_animal_trunc{:});

vis_workflow = 'stim_wheel_right_stage\d$';
aud_workflow = 'stim_wheel_right_stage\d_audio*';
va_animals = cellfun(@(x) max([...
    find(cellfun(@any,regexp(x,vis_workflow)),1) < ...
    find(cellfun(@any,regexp(x,aud_workflow)),1), ...
    false]),workflow_animal);
av_animals = cellfun(@(x) max([...
    find(cellfun(@any,regexp(x,vis_workflow)),1) > ...
    find(cellfun(@any,regexp(x,aud_workflow)),1), ...
    false]),workflow_animal);

% (ALT: define manually)
va_animals = ismember({bhv.animal},{'AP019','AP021','AP022','DS001','DS007','DS010','DS011'})';
av_animals = ismember({bhv.animal},{'DS000','DS004','DS014','DS015','DS016','DS005'})';

va_animals_cat = cell2mat(cellfun(@(x,y) repmat(x,length(y),1),num2cell(va_animals),workflow_animal_trunc,'uni',false));
av_animals_cat = cell2mat(cellfun(@(x,y) repmat(x,length(y),1),num2cell(av_animals),workflow_animal_trunc,'uni',false));

% Get learned day
v_learned_day = cellfun(@(x,p,stat) ...
     sum(cellfun(@any,regexp( ...
     x(1:find(cellfun(@any,regexp(x,vis_workflow)) & p < 0.05 & stat < 1,1)),vis_workflow))), ...
     workflow_animal,{bhv.rxn_stat_p}',{bhv.rxn_stat}');
 
a_learned_day = cellfun(@(x,p,stat) ...
     sum(cellfun(@any,regexp( ...
     x(1:find(cellfun(@any,regexp(x,aud_workflow)) & p < 0.05 & stat < 1,1)),aud_workflow))), ...
     workflow_animal,{bhv.rxn_stat_p}',{bhv.rxn_stat}');

v_ld = cell2mat(cellfun(@(x,ld) (1:length(x))'-ld,stim_V,num2cell(v_learned_day),'uni',false));
a_ld = cell2mat(cellfun(@(x,ld) (1:length(x))'-ld,stim_V,num2cell(a_learned_day),'uni',false));

% Movie V/A
task_workflow = 'stim_wheel_right_stage\d$';
use_rec = cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat < 0.05;
v = plab.wf.svd2px(U_master(:,:,1:n_components),nanmean(cat(4,stim_V_cat{use_rec,v_col}),4));

task_workflow = 'stim_wheel_right_stage\d_audio*';
use_rec = cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat < 0.05;
a = plab.wf.svd2px(U_master(:,:,1:n_components),nanmean(cat(4,stim_V_cat{use_rec,a_col}),4));

ap.imscroll(v);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));

ap.imscroll(a);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));



% DRAW ROI
plot_days = -3:5;

col = ap.colormap('BKR',max(abs(plot_days))*2+1);
col_days = -max(abs(plot_days)):max(abs(plot_days));

task_workflow = 'stim_wheel_right_stage\d$';
use_rec = va_animals_cat & cellfun(@any,regexp(workflow_cat,task_workflow)) & ismember(v_ld,plot_days);
v_roi = cell2mat(cellfun(@(x) ap.wf_roi(U_master(:,:,1:n_components), ...
    x(:,:,v_stim),[],[],roi.mask),stim_V_cat(use_rec,v_col), ...
    'uni',false,'ErrorHandler',@(varargin) nan(1,length(t),'single')));

[v_roi_mean,v_roi_group] = ap.groupfun(@nanmean,v_roi,v_ld(use_rec));
v_roi_sem = ap.groupfun(@AP_sem,v_roi,v_ld(use_rec));

figure; h = tiledlayout(1,length(plot_days));
for curr_idx = 1:size(v_roi_mean,1)
    nexttile;
    ap.errorfill(t,v_roi_mean(curr_idx,:),v_roi_sem(curr_idx,:), ...
        col(ismember(col_days,v_roi_group(curr_idx)),:));
end
linkaxes(h.Children);


task_workflow = 'stim_wheel_right_stage\d_audio*';
use_rec = va_animals_cat & cellfun(@any,regexp(workflow_cat,task_workflow)) & ismember(a_ld,plot_days);
a_roi = cell2mat(cellfun(@(x) ap.wf_roi(U_master(:,:,1:n_components), ...
    x(:,:,a_stim),[],[],roi.mask),stim_V_cat(use_rec,a_col), ...
    'uni',false,'ErrorHandler',@(varargin) nan(1,length(t),'single')));

[a_roi_mean,a_roi_group] = ap.groupfun(@mean,a_roi,a_ld(use_rec));
a_roi_sem = ap.groupfun(@AP_sem,a_roi,a_ld(use_rec));

figure; h = tiledlayout(1,length(plot_days));
for curr_idx = 1:size(a_roi_mean,1)
    nexttile;
    ap.errorfill(t,a_roi_mean(curr_idx,:),a_roi_sem(curr_idx,:), ...
        col(ismember(col_days,a_roi_group(curr_idx)),:));
end
linkaxes(h.Children);




%% Task maps 

kernel_fn = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\task_kernel.mat';
load(kernel_fn);
t = (0:20)/35;
n_components = size(stim_V{1}{1},1);

bhv_filename = 'C:\Users\petersa\Documents\PetersLab\analysis\av_association\data\bhv.mat';
load(bhv_filename);

U_master = plab.wf.load_master_U;

%%% PROBLEM! bhv doesn't match kernels (sometime extra bhv day), don't know
%%% why at the moment, maybe just truncated if error on last day, but
%%% didn't see error on example one

stim_V_cat = horzcat(stim_V{:})';

% just truncate bhv values for now
rxn_stat_p_animal = {bhv.rxn_stat_p}';
rxn_stat_p_animal_trunc = cellfun(@(x,y) x(1:size(y,2))',rxn_stat_p_animal,stim_V,'uni',false);
rxn_stat_p_cat = vertcat(rxn_stat_p_animal_trunc{:});

workflow_animal = {bhv.workflow}';
workflow_animal_trunc = cellfun(@(x,y) x(1:size(y,2))',workflow_animal,stim_V,'uni',false);
workflow_cat = vertcat(workflow_animal_trunc{:});

vis_workflow = 'stim_wheel_right_stage\d$';
aud_workflow = 'stim_wheel_right_stage\d_audio*';
va_animals = cellfun(@(x) max([...
    find(cellfun(@any,regexp(x,vis_workflow)),1) < ...
    find(cellfun(@any,regexp(x,aud_workflow)),1), ...
    false]),workflow_animal);
av_animals = cellfun(@(x) max([...
    find(cellfun(@any,regexp(x,vis_workflow)),1) > ...
    find(cellfun(@any,regexp(x,aud_workflow)),1), ...
    false]),workflow_animal);

% (ALT: define manually)
% va_animals = ismember({bhv.animal},{'AP019','AP021','AP022','DS001','DS007','DS010','DS011'})';
va_animals = ismember({bhv.animal},{'AP019','AP021','AP022','DS001','DS007','DS010','DS011','DS019','DS020','DS021'})';
av_animals = ismember({bhv.animal},{'DS000','DS004','DS014','DS015','DS016','DS005'})';

va_animals_cat = cell2mat(cellfun(@(x,y) repmat(x,length(y),1),num2cell(va_animals),workflow_animal_trunc,'uni',false));
av_animals_cat = cell2mat(cellfun(@(x,y) repmat(x,length(y),1),num2cell(av_animals),workflow_animal_trunc,'uni',false));

% Plot
task_workflow = 'stim_wheel_right_stage\d$';
task_workflow = 'stim_wheel_right_stage\d_audio';
% use_rec = cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat < 0.05;
use_rec = va_animals_cat & cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat < 0.05;
v_learned = plab.wf.svd2px(U_master(:,:,1:n_components),nanmean(cat(4,stim_V_cat{use_rec}),4));
% v_learned = imgaussfilt(v_learned,3)-imgaussfilt(v_learned,20);

task_workflow = 'stim_wheel_right_stage\d_audio';
% use_rec =  cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat < 0.05;
use_rec = av_animals_cat & cellfun(@any,regexp(workflow_cat,task_workflow)) & rxn_stat_p_cat < 0.05;
a_learned = plab.wf.svd2px(U_master(:,:,1:n_components),nanmean(cat(4,stim_V_cat{use_rec}),4));
% a_learned = imgaussfilt(a_learned,3)-imgaussfilt(a_learned,20);

v_col = ap.colormap('WG');
a_col = ap.colormap('WP');
col_lim = [0,3e-4]/2; % kernel

v_learned_gray = 1+round(mat2gray(v_learned,col_lim).*(size(v_col,1)-1));
v_learned_col = permute(reshape(v_col(v_learned_gray,:),[size(v_learned_gray),3]),[1,2,4,3]);

a_learned_gray = 1+round(mat2gray(a_learned,col_lim).*(size(a_col,1)-1));
a_learned_col = permute(reshape(a_col(a_learned_gray,:),[size(a_learned_gray),3]),[1,2,4,3]);

va_learned_col = min(v_learned_col,a_learned_col);

ap.imscroll(va_learned_col,t,true);
axis image;
ap.wf_draw('ccf','k');

plot_frames = find(t>=0 & t<= 0.2);
figure; h = tiledlayout(3,length(plot_frames),'TileSpacing','none','TileIndexing','ColumnMajor');
for curr_frame = plot_frames
    nexttile;
    image(v_learned_col(:,:,:,curr_frame));
    ap.wf_draw('ccf','k');
    axis image off;
    title(sprintf('%.2f',(curr_frame-1)/35));

    nexttile;
    image(a_learned_col(:,:,:,curr_frame));
    ap.wf_draw('ccf','k');
    axis image off;

    nexttile;
    image(va_learned_col(:,:,:,curr_frame));
    ap.wf_draw('ccf','k');
    axis image off;
end


% Plot max V/A/overlay
stim_t = t > 0 & t < 0.2;
v_learned_max = max(v_learned(:,:,stim_t),[],3);
a_learned_max = max(a_learned(:,:,stim_t),[],3);

figure;
h = tiledlayout(1,3,'TileSpacing','none');

nexttile; imagesc(v_learned_max); 
clim(col_lim); axis image off;
colormap(gca,v_col);
ap.wf_draw('ccf','k');

nexttile; imagesc(a_learned_max); 
clim(col_lim); axis image off;
colormap(gca,a_col);
ap.wf_draw('ccf','k');

v_learned_max_gray = 1+round(mat2gray(v_learned_max,col_lim).*(size(v_col,1)-1));
v_learned_max_col = reshape(v_col(v_learned_max_gray,:),size(v_learned_max,1),size(v_learned_max,2),3,[]);

a_learned_max_gray = 1+round(mat2gray(a_learned_max,col_lim).*(size(a_col,1)-1));
a_learned_max_col = reshape(a_col(a_learned_max_gray,:),size(a_learned_max,1),size(a_learned_max,2),3,[]);

va_learned_max_col = min(v_learned_max_col,a_learned_max_col);
nexttile; image(va_learned_max_col);
axis image off;
ap.wf_draw('ccf','k');


%% ~~~~~~~~~~~~~ Specific animals, grand average

%% Task-specific: AVERAGE passive stim

% Set task workflow
% task_workflow = 'stim_wheel_right_stage\d_size_up';
% task_workflow = 'stim_wheel_right_stage\d_angle_size60';
% task_workflow = 'stim_wheel_right_stage\d_audio_volume';
task_workflow = 'stim_wheel_right_stage\d_audio*';
task_workflow = 'stim_wheel_right_stage\d';

% Set passive workflow
% passive_workflow = 'lcr_passive';
passive_workflow = 'hml_passive_audio';

% animals = {'AP021','AP022','DS001','DS007','DS010','DS011'}; % V>A
% animals = {'DS000','DS003','DS004','DS013','DS014','DS015','DS016'}; % A>V
% animals = {'DS006','DS013'}; % A>no-learn V

animals = { ...
    'AP019','AP021','AP022','DS001','DS007','DS010','DS011', ... % V > A
    'DS000','DS003','DS004','DS014','DS015','DS016','DS005', ... % A > V
    'DS019','DS020','DS021','AP027','AP028','AP029', ... % V-alt > V > Afreq
    'AP020','AP018', ... % V > A-nonlearn
    'DS006','DS013', ... % A > V-nonlearn
    'HA003','HA004', ... % V-alt > V
    'HA000','HA001','HA002' ... % V-alt > V
    };

% Grab V's from passive stimuli on post-learning days
stim_v = cell(size(animals));
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    % Find learned days
    task_recordings = plab.find_recordings(animal,[],task_workflow);
    curr_rxn_stat_p = nan(size(task_recordings));

    for curr_recording = 1:length(task_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = task_recordings(curr_recording).day;
        rec_time = task_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        if n_trials < 10
            continue
        end

        % Get reaction stat
        use_stat = 'mean';
        curr_rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    % Grab passive LCR response from learned days
    learned_passive_recordings = plab.find_recordings(animal, ...
        {task_recordings(curr_rxn_stat_p < 0.05).day},passive_workflow);

    for curr_recording = 1:length(learned_passive_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = learned_passive_recordings(curr_recording).day;
        rec_time = learned_passive_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        % Stim-align data (quiescent trials)
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        trial_t_window = [-0.3,1];
        trial_t_rate = 35;
        trial_t = trial_t_window(1):1/trial_t_rate:trial_t_window(2);
        interp_t = reshape(stimOn_times(quiescent_trials),[],1) + trial_t;

        n_wf_components = 200;
        wf_V_stimalign = interp1(wf_t,wf_V(1:n_wf_components,:)',interp_t,'previous');

        switch passive_workflow
            case 'lcr_passive'
                stim_type = vertcat(trial_events.values.TrialStimX);
            case 'hml_passive_audio'
                stim_type = vertcat(trial_events.values.StimFrequence);
        end

        wf_V_stimalign_avg = nan(n_wf_components,length(trial_t),length(unique(stim_type)),class(wf_V_stimalign));
        stim_types_used_idx = ismember(unique(stim_type),unique(stim_type(quiescent_trials)));

        wf_V_stimalign_avg(:,:,stim_types_used_idx) = ...
            permute(ap.groupfun(@mean,wf_V_stimalign,stim_type(quiescent_trials),[],[]),[3,2,1]);

        stim_v{curr_animal}{curr_recording} = wf_V_stimalign_avg;

        ap.print_progress_fraction(curr_recording,length(learned_passive_recordings));
    end

end
disp('Done.');

U_master = plab.wf.load_master_U;

stim_v_avg = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(cat(4,x{:}),4),stim_v(~cellfun(@isempty,stim_v)),'uni',false),[1,3,4,2])),4);

stim_v_avg = stim_v_avg - nanmean(stim_v_avg(:,trial_t<0,:),2);

stim_v_px = plab.wf.svd2px(U_master(:,:,1:n_wf_components),stim_v_avg);
ap.imscroll(stim_v_px,trial_t);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(AP_colormap('PWG'));
ap.wf_draw('ccf','k');


%% Task-specific: KERNEL passive stim (linear)

% task_workflow = 'stim_wheel_right_stage\d';
% task_workflow = 'stim_wheel_right_stage\d_size_up';
% task_workflow = 'stim_wheel_right_stage\d_angle_size60';
task_workflow = 'stim_wheel_right_stage\d_audio_volume';
% task_workflow = 'stim_wheel_right_stage\d_audio_frequency';

% passive_workflow = 'lcr_passive';
passive_workflow = 'hml_passive_audio';

% animals = {'AP021','AP022','DS001','DS007','DS010','DS011'}; % V>A
% animals = {'DS000','DS003','DS004','DS013','DS014','DS015','DS016'}; % A>V
% animals = {'DS019','DS020','DS021','AP027','AP028','AP029'}; % V>A freq
% animals = {'DS006','DS013'}; % A>no-learn V
animals = {'AP018','AP020'}; % V>no-learn A

% Grab V's from passive stimuli on post-learning days
stim_v = cell(size(animals));
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    % Find learned days
    task_recordings = plab.find_recordings(animal,[],task_workflow);
    curr_rxn_stat_p = nan(size(task_recordings));

    for curr_recording = 1:length(task_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = task_recordings(curr_recording).day;
        rec_time = task_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get reaction stat
        use_stat = 'mean';
        curr_rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    % Grab passive LCR response from learned days
    learned_passive_recordings = plab.find_recordings(animal, ...
        {task_recordings(curr_rxn_stat_p < 0.05).day},passive_workflow);

    for curr_recording = 1:length(learned_passive_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = learned_passive_recordings(curr_recording).day;
        rec_time = learned_passive_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        switch bonsai_workflow
            case 'lcr_passive'
                stim_type = vertcat(trial_events.values.TrialStimX);
            case 'hml_passive_audio'
                stim_type = vertcat(trial_events.values.StimFrequence);
        end

        time_bins = [wf_t;wf_t(end)+1/wf_framerate];

        stim_regressors = cell2mat(arrayfun(@(x) ...
            histcounts(stimOn_times(stim_type == x),time_bins), ...
            unique(stim_type),'uni',false));

        n_components = 500;

        frame_shifts = 0:30;
        lambda = 20;

        skip_t = 3; % seconds start/end to skip for artifacts
        skip_frames = round(skip_t*wf_framerate);
        [kernels,predicted_signals,explained_var] = ...
            ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
            stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda);

        stim_v{curr_animal}{curr_recording} = kernels;

        ap.print_progress_fraction(curr_recording,length(learned_passive_recordings));
    end

end
disp('Done.');

U_master = plab.wf.load_master_U;

stim_v_avg = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(cat(4,x{:}),4),stim_v,'uni',false),[1,3,4,2])),4);

stim_v_px = plab.wf.svd2px(U_master(:,:,1:size(stim_v_avg,1)),stim_v_avg);
ap.imscroll(stim_v_px);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(AP_colormap('PWG'));
ap.wf_draw('ccf','k');



% x = cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:200),x),cat(2,stim_v{:}),'uni',false);
% x2 = nanmedian(cell2mat(permute(x,[1,3,4,5,2])),5);
% ap.imscroll(x2);
% axis image;
% clim(max(abs(clim)).*[-1,1]);
% colormap(AP_colormap('PWG'));
% ap.wf_draw('ccf','k');

%% Task-specific: KERNEL task stim (linear)

% task_workflow = 'stim_wheel_right_stage\d_size_up';
% task_workflow = 'stim_wheel_right_stage\d_angle_size60';
% task_workflow = 'stim_wheel_right_stage\d_audio_volume';
task_workflow = 'stim_wheel_right_stage\d_audio*';
% task_workflow = 'stim_wheel_right_stage\d';
% task_workflow = '*audio_frequency*';

% animals = {'DS019','DS020','DS021','HA003','HA004'};

% animals = {'AP021','AP022','DS001','DS007','DS010','DS011'}; % V>A
% animals = {'DS000','DS003','DS004','DS013','DS014','DS015','DS016'}; % A>V
% animals = {'DS019','DS020','DS021','AP027','AP028','AP029'}; % V > A freq
% animals = {'DS006','DS013'}; % A>no-learn V

animals = { ...
    'AP019','AP021','AP022','DS001','DS007','DS010','DS011', ...
    'DS000','DS003','DS004','DS014','DS015','DS016','DS005', ...
    'DS019','DS020','DS021','AP027','AP028','AP029', ...
    'AP020','AP018', ...
    'DS006','DS013', ...
    'HA003','HA004', ...
    'HA000','HA001','HA002'};

% Grab V's from passive stimuli on post-learning days
stim_v = cell(size(animals));
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    % Find learned days
    task_recordings = plab.find_recordings(animal,[],task_workflow);
    curr_rxn_stat_p = nan(size(task_recordings));

    for curr_recording = 1:length(task_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = task_recordings(curr_recording).day;
        rec_time = task_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        if n_trials < 10
            continue
        end

        % Get reaction stat
        use_stat = 'mean';
        curr_rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    % Grab passive LCR response from learned days
    learned_recordings = plab.find_recordings(animal, ...
        {task_recordings(curr_rxn_stat_p < 0.05).day},task_workflow);

    for curr_recording = 1:length(learned_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = learned_recordings(curr_recording).day;
        rec_time = learned_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        time_bins = [wf_t;wf_t(end)+1/wf_framerate];

        stim_regressors = histcounts(stimOn_times,time_bins);

        n_components = 200;

        frame_shifts = 0:20;
        lambda = 20;
        cv_fold = 3;

        skip_t = 60; % seconds start/end to skip for artifacts
        skip_frames = round(skip_t*wf_framerate);
        [kernels,predicted_signals,explained_var] = ...
            ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
            stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

        stim_v{curr_animal}{curr_recording} = kernels;

        ap.print_progress_fraction(curr_recording,length(learned_recordings));
    end

end
disp('Done.');

U_master = plab.wf.load_master_U;

stim_v_avg = nanmean(cell2mat(permute(horzcat(stim_v{:}),[1,3,2])),3);

stim_v_px = plab.wf.svd2px(U_master(:,:,1:size(stim_v_avg,1)),stim_v_avg);
ap.imscroll(stim_v_px);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(AP_colormap('PWG'));
ap.wf_draw('ccf','k');


% x = cellfun(@(x) plab.wf.svd2px(U_master(:,:,1:200),x),cat(2,stim_v{:}),'uni',false);
% x2 = nanmedian(cell2mat(permute(x,[1,3,4,5,2])),5);
% ap.imscroll(x2);
% axis image;
% clim(max(abs(clim)).*[-1,1]);
% colormap(AP_colormap('PWG'));
% ap.wf_draw('ccf','k');


%% Task-specific: KERNEL passive stim (logistic)

% task_workflow = 'stim_wheel_right_stage\d_size_up';
% task_workflow = 'stim_wheel_right_stage\d_angle_size60';
task_workflow = 'stim_wheel_right_stage\d_audio_volume';
% task_workflow = 'stim_wheel_right_stage\d';

% passive_workflow = 'lcr_passive';
passive_workflow = 'hml_passive_audio';

% animals = {'DS019','DS020','DS021','HA003','HA004'};

animals = {'AP021','AP022','DS001','DS007','DS010','DS011'}; % V>A
% animals = {'DS000','DS003','DS004','DS013','DS014','DS015','DS016'}; % A>V


% Grab V's from passive stimuli on post-learning days
stim_v = cell(size(animals));
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    % Find learned days
    task_recordings = plab.find_recordings(animal,[],task_workflow);
    curr_rxn_stat_p = nan(size(task_recordings));

    for curr_recording = 1:length(task_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = task_recordings(curr_recording).day;
        rec_time = task_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % (temp: skip if <10 trials - usually this is because short bad one
        % and long good one, not always last - find the good one)
        if n_trials < 10
            continue
        end

        % Get reaction stat
        use_stat = 'mean';
        curr_rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    % Grab passive LCR response from learned days
    learned_passive_recordings = plab.find_recordings(animal, ...
        {task_recordings(curr_rxn_stat_p < 0.05).day},passive_workflow);

    for curr_recording = 1:length(learned_passive_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = learned_passive_recordings(curr_recording).day;
        rec_time = learned_passive_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        switch bonsai_workflow
            case 'lcr_passive'
                stim_type = vertcat(trial_events.values.TrialStimX);
            case 'hml_passive_audio'
                stim_type = vertcat(trial_events.values.StimFrequence);
        end


        switch bonsai_workflow
            case {'lcr_passive','lcr_passive_size60'}
                stim_type = vertcat(trial_events.values.TrialStimX);
                use_stim = 3;
            case 'hml_passive_audio'
                stim_type = vertcat(trial_events.values.StimFrequence);
                use_stim = 2;
        end

        time_bins = [wf_t;wf_t(end)+1/wf_framerate];

        stim_regressors = cell2mat(arrayfun(@(x) ...
            histcounts(stimOn_times(stim_type == x),time_bins), ...
            unique(stim_type),'uni',false));

        n_components = 100;

        timesteps = -10:30;
        regressors = repmat(wf_V(1:n_components,:)',1,1,length(timesteps));
        for i = 1:length(timesteps)
            regressors(:,:,i) = circshift(regressors(:,:,i),-timesteps(i),1);
        end

        regressors_flat = reshape(regressors,[],prod(size(regressors,[2,3])));

        % Model - elastic net (can't get pure ridge alpha=0)
        [B, FitInfo] = lassoglm(regressors_flat, stim_regressors(use_stim,:)','binomial','Alpha',1e-5,'Lambda',1e-2);

        kernels = reshape(B,n_components,[]);

        stim_v{curr_animal}{curr_recording} = kernels;

        ap.print_progress_fraction(curr_recording,length(learned_passive_recordings));
    end

end
disp('Done.');

U_master = plab.wf.load_master_U;

stim_v_avg = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(cat(4,x{:}),4),stim_v,'uni',false),[1,3,4,2])),4);

stim_v_px = plab.wf.svd2px(U_master(:,:,1:size(stim_v_avg,1)),stim_v_avg);
ap.imscroll(stim_v_px);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(AP_colormap('PWG'));
ap.wf_draw('ccf','k');


%% Task-specific: KERNEL passive stim (logistic classification)

% task_workflow = 'stim_wheel_right_stage\d_size_up';
% task_workflow = 'stim_wheel_right_stage\d_angle_size60';
task_workflow = 'stim_wheel_right_stage\d_audio_volume';
% task_workflow = 'stim_wheel_right_stage\d';

% passive_workflow = 'lcr_passive';
passive_workflow = 'hml_passive_audio';

% animals = {'DS019','DS020','DS021','HA003','HA004'};

% animals = {'AP021','AP022','DS001','DS007','DS010','DS011'}; % V>A
animals = {'DS000','DS003','DS004','DS013','DS014','DS015','DS016'}; % A>V


% Grab V's from passive stimuli on post-learning days
stim_v = cell(size(animals));
for curr_animal = 1:length(animals)

    animal = animals{curr_animal};
    fprintf([animal,'\n']);

    % Find learned days
    task_recordings = plab.find_recordings(animal,[],task_workflow);
    curr_rxn_stat_p = nan(size(task_recordings));

    for curr_recording = 1:length(task_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = task_recordings(curr_recording).day;
        rec_time = task_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % (temp: skip if <10 trials - usually this is because short bad one
        % and long good one, not always last - find the good one)
        if n_trials < 10
            continue
        end

        % Get reaction stat
        use_stat = 'mean';
        curr_rxn_stat_p(curr_recording) = AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

    end

    % Grab passive LCR response from learned days
    learned_passive_recordings = plab.find_recordings(animal, ...
        {task_recordings(curr_rxn_stat_p < 0.05).day},passive_workflow);

    for curr_recording = 1:length(learned_passive_recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = learned_passive_recordings(curr_recording).day;
        rec_time = learned_passive_recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.widefield = true;
        load_parts.widefield_master = true;
        ap.load_recording;

        switch bonsai_workflow
            case 'lcr_passive'
                stim_type = vertcat(trial_events.values.TrialStimX);
            case 'hml_passive_audio'
                stim_type = vertcat(trial_events.values.StimFrequence);
        end


        switch bonsai_workflow
            case {'lcr_passive','lcr_passive_size60'}
                stim_type = vertcat(trial_events.values.TrialStimX);
                use_stim = 3;
            case 'hml_passive_audio'
                stim_type = vertcat(trial_events.values.StimFrequence);
                use_stim = 2;
        end


        % max snapshot, lassoglm (doesn't look good)
        [~,~,stim_type_unique] = unique(stim_type);

        stim_window = stimOn_times(stim_type_unique == use_stim) + [0,0.2];
        stim_max_px = cell2mat(arrayfun(@(x) max(plab.wf.svd2px(wf_U, ...
            wf_V(:,wf_t >= stim_window(x,1) ...
            & wf_t <= stim_window(x,2))),[],3), ...
            permute(1:size(stim_window,1),[1,3,2]),'uni',false));

        baseline_window = stimOn_times + [-0.3,0];
        baseline_max_px = cell2mat(arrayfun(@(x) max(plab.wf.svd2px(wf_U, ...
            wf_V(:,wf_t >= baseline_window(x,1) ...
            & wf_t <= baseline_window(x,2))),[],3), ...
            permute(1:size(baseline_window,1),[1,3,2]),'uni',false));

        wf_snapshots = reshape(wf_U,[],size(wf_U,3))'* ...
            cat(2,reshape(stim_max_px,[],size(stim_max_px,3)), ...
            reshape(baseline_max_px,[],size(baseline_max_px,3)));
        stim_shapshot = [true(size(stim_window,1),1);false(size(baseline_window,1),1)];


        n_components = 2000;

        model = fitclinear(wf_snapshots(1:n_components,:)', stim_shapshot, ...
            'Learner', 'logistic', 'Regularization', 'ridge', 'Lambda', 0);

        kernels = model.Beta;

        stim_v{curr_animal}{curr_recording} = kernels;

        ap.print_progress_fraction(curr_recording,length(learned_passive_recordings));
    end

end
disp('Done.');

U_master = plab.wf.load_master_U;

stim_v_avg = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(cat(4,x{:}),4),stim_v,'uni',false),[1,3,4,2])),4);

stim_v_px = plab.wf.svd2px(U_master(:,:,1:size(stim_v_avg,1)),stim_v_avg);
ap.imscroll(stim_v_px);
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(AP_colormap('PWG'));
ap.wf_draw('ccf','k');


%% Testing kernel parameters

% % VA V-learn
% animal = 'DS011';
% rec_day = '2024-07-20';
% rec_time = '1755';
% verbose = true;
% ap.load_recording;

% VA A-learn
animal = 'DS011';
rec_day = '2024-07-28';
rec_time = '1042';
verbose = true;
ap.load_recording;

%% (testing linear)

switch bonsai_workflow
    case {'lcr_passive','lcr_passive_size60'}
        stim_type = vertcat(trial_events.values.TrialStimX);
    case 'hml_passive_audio'
        stim_type = vertcat(trial_events.values.StimFrequence);
end


sample_rate = 100;
sample_t = (wf_t(1):1/sample_rate:wf_t(end))';

wf_V_resample = interp1(wf_t,wf_V',sample_t,'previous')';

time_bins = [sample_t;sample_t(end)+1/sample_rate];

stim_regressors = cell2mat(arrayfun(@(x) ...
    histcounts(stimOn_times(stim_type == x),time_bins), ...
    unique(stim_type),'uni',false));

% Regular linear regression (to binary variable, so bad practice)
n_components = 200;

t_shifts = [0,0.5];
frame_shifts = round(t_shifts(1)*sample_rate):round(t_shifts(2)*sample_rate);
lambda = 20;

[kernels,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V_resample(1:n_components,:),stim_regressors,-frame_shifts,lambda);

ap.imscroll(plab.wf.svd2px(wf_U(:,:,1:size(kernels,1)),kernels));
axis image;
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));



%% (testing logistic)

n_components = 100;

timesteps = -10:30;
regressors = repmat(wf_V(1:n_components,:)',1,1,length(timesteps));
for i = 1:length(timesteps)
    regressors(:,:,i) = circshift(regressors(:,:,i),-timesteps(i),1);
end

regressors_flat = reshape(regressors,[],prod(size(regressors,[2,3])));

use_stim = 2;

% GLM - not regularized
model = fitglm(regressors_flat, stim_regressors(use_stim,:)','linear','Distribution','binomial','Link','logit');

px = plab.wf.svd2px(wf_U(:,:,1:n_components),reshape(model.Coefficients.Estimate(2:end),n_components,[]));
ap.imscroll(px)


% Model - elastic net (can't get pure ridge alpha=0)
[B, FitInfo] = lassoglm(regressors_flat, stim_regressors(use_stim,:)','binomial','Alpha',1e-5,'Lambda',1e-2);

px = plab.wf.svd2px(wf_U(:,:,1:n_components),reshape(B,n_components,[]));
ap.imscroll(px); axis image;
colormap(ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);


% Model - elastic net (can't get pure ridge alpha=0) (do in chunks?)
n_split = 4;
split_idx = ceil((1:length(wf_t))/(length(wf_t)/n_split));

B = nan(size(regressors_flat,2),n_split,'single');
for curr_split = 1:n_split
    B(:,curr_split) = lassoglm(regressors_flat(split_idx==curr_split,:), ...
        stim_regressors(use_stim,split_idx==curr_split)','binomial','Alpha',1e-5,'Lambda',1e-2);
    ap.print_progress_fraction(curr_split,n_split)
end

px = plab.wf.svd2px(wf_U(:,:,1:n_components),reshape(B,n_components,length(timesteps),[]));
ap.imscroll(median(px,4)); axis image;
colormap(ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);


% Exclude other stim times (didn't change anything)
alt_stim_times = stimOn_times(~ismember(stim_type,90)) + [0,1];
use_t = mod(discretize(wf_t,reshape(alt_stim_times',[],1)),2) ~= 1;

[B, FitInfo] = lassoglm(regressors_flat(use_t,:), stim_regressors(use_stim,use_t)','binomial','Alpha',1e-7,'Lambda',1e-1,'CV',5);

px = plab.wf.svd2px(wf_U(:,:,1:n_components),reshape(B,n_components,[]));
ap.imscroll(px); axis image;


% fitclinear (this looks just exactly like the average?)
model = fitclinear(regressors_flat, stim_regressors(use_stim,:)', 'Learner', 'logistic', 'Regularization', 'ridge', 'Lambda', 1);
px = plab.wf.svd2px(wf_U(:,:,1:n_components),reshape(model.Beta,n_components,[]));
ap.imscroll(px); axis image;
colormap(ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);


% max snapshot, lassoglm (doesn't look good)
[~,~,stim_type_unique] = unique(stim_type);

stim_window = stimOn_times(stim_type_unique == use_stim) + [0,0.2];
stim_max_px = cell2mat(arrayfun(@(x) max(plab.wf.svd2px(wf_U, ...
    wf_V(:,wf_t >= stim_window(x,1) ...
    & wf_t <= stim_window(x,2))),[],3), ...
    permute(1:size(stim_window,1),[1,3,2]),'uni',false));

baseline_window = stimOn_times + [-0.3,0];
baseline_max_px = cell2mat(arrayfun(@(x) max(plab.wf.svd2px(wf_U, ...
    wf_V(:,wf_t >= baseline_window(x,1) ...
    & wf_t <= baseline_window(x,2))),[],3), ...
    permute(1:size(baseline_window,1),[1,3,2]),'uni',false));

wf_snapshots = reshape(wf_U,[],size(wf_U,3))'* ...
    cat(2,reshape(stim_max_px,[],size(stim_max_px,3)), ...
    reshape(baseline_max_px,[],size(baseline_max_px,3)));
stim_shapshot = [true(size(stim_window,1),1);false(size(baseline_window,1),1)];

n_components = 100;
[B, FitInfo] = lassoglm(wf_snapshots(1:n_components,:)', stim_shapshot,'binomial','Alpha',0.1,'Lambda',1);
px = plab.wf.svd2px(wf_U(:,:,1:n_components),B);

figure;
imagesc(px); 
axis image;
colormap(ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);


% max snapshot, fitclinear (doesn't look good) 
n_components = 2000;

model = fitclinear(wf_snapshots(1:n_components,:)', stim_shapshot, ...
    'Learner', 'logistic', 'Regularization', 'ridge', 'Lambda', 0);
px = plab.wf.svd2px(wf_U(:,:,1:n_components),model.Beta);
figure;imagesc(px); axis image;
colormap(ap.colormap('PWG'));
clim(max(abs(clim)).*[-1,1]);













