%% Analyses for reviewers (cleaned up from _1 scratch)

%% R1 p1: Correct performance rate by day

%%% Load non-activity data
load_dataset = 'noact';
Marica_2025.figures.load_data;
%%%

% Plot reaction time and association index, split within day
n_daysplit = 3;

outcome_mean_daysplit = cell2mat(cellfun(@(x,idx) ap.groupfun(@mean,+x, ...
    ap.quantile_bin(length(x),n_daysplit)),bhv.trial_outcome,'uni',false)')';

[outcome_daysplit_mean,outcome_group_x] = ...
    ap.groupfun(@mean,outcome_mean_daysplit,bhv.days_from_learning);
outcome_daysplit_sem = ap.groupfun(@AP_sem,outcome_mean_daysplit, ...
    bhv.days_from_learning);

plot_days = -3:2;
plot_day_idx = ismember(outcome_group_x,plot_days);

figure;
outcome_group_x_daysplit = outcome_group_x+(0:n_daysplit)./n_daysplit;
errorbar(reshape(outcome_group_x_daysplit(plot_day_idx,:)',[],1), ...
    reshape(padarray(outcome_daysplit_mean(plot_day_idx,:),[0,1],nan,'post')',[],1), ...
    reshape(padarray(outcome_daysplit_sem(plot_day_idx,:),[0,1],nan,'post')',[],1),'k','linewidth',2);
xline(0,'r');
ylabel('Fraction correct');
xlabel('Day from learning');
ap.prettyfig;


%% R2 p4: CCF-aligned probe histology

[av,tv,st] = ap_histology.load_ccf;

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

% Set atlas bins to plot through striatum
striatum_ccf_id = find(contains(lower(st.safe_name),'caudoputamen'));
striatum_ccf_ap = prctile(find(max(av == striatum_ccf_id,[],[2,3])),[0,100]);

n_atlas_bins = 10;
atlas_bins = round(linspace(striatum_ccf_ap(1),striatum_ccf_ap(2),n_atlas_bins+1));

histology_atlas_bin_max = cell(length(animals),n_atlas_bins);

for curr_animal_idx = 1:length(animals)

    histology_path = plab.locations.filename('server',animals{curr_animal_idx},[],[],'histology');
    load(fullfile(histology_path,'AP_histology_processing.mat'));

    % Load images
    image_path = histology_path;
    image_dir = dir(fullfile(image_path,'*.tif'));
    image_filenames = cellfun(@(path,name) fullfile(path,name), ...
        {image_dir.folder},{image_dir.name},'uni',false);
    [~,sort_idx] = ap_histology.natsortfiles(image_filenames);

    images = cell(size(image_dir));
    for curr_im = 1:length(sort_idx)
        images{curr_im} = tiffreadVolume( ...
            image_filenames{sort_idx(curr_im)});
    end

    % Grab atlas images
    n_slices = length(images);
    slice_atlas = struct('tv',cell(n_slices,1), 'av',cell(n_slices,1));
    slice_atlas_ccf = struct('ap',cell(n_slices,1),'ml',cell(n_slices,1),'dv',cell(n_slices,1));
    for curr_slice = 1:length(images)
        [slice_atlas(curr_slice),slice_atlas_ccf(curr_slice)] = ...
            ap_histology.grab_atlas_slice(av,tv, ...
            AP_histology_processing.histology_ccf.slice_vector, ...
            AP_histology_processing.histology_ccf.slice_points(curr_slice,:), 1);
    end

    % Build volume of histology images
    histology_volume = zeros(size(tv),'single');
    probe_channel = 2;
    for curr_im_idx = 1:length(images)

        % Rigid transform
        im_rigid_transformed = ap_histology.rigid_transform( ...
            images{curr_im_idx}(:,:,probe_channel),curr_im_idx,AP_histology_processing);

        % Affine/nonlin transform
        if isfield(AP_histology_processing.histology_ccf,'control_points') && ...
                (size(AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx},1) == ...
                size(AP_histology_processing.histology_ccf.control_points.atlas{curr_im_idx},1)) && ...
                size(AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx},1) >= 3
            % Manual alignment (if >3 matched points)
            histology2atlas_tform = fitgeotform2d( ...
                AP_histology_processing.histology_ccf.control_points.histology{curr_im_idx}, ...
                AP_histology_processing.histology_ccf.control_points.atlas{curr_im_idx},'pwl');
        elseif isfield(AP_histology_processing.histology_ccf,'atlas2histology_tform')
            % Automatic alignment
            histology2atlas_tform = invert(AP_histology_processing.histology_ccf.atlas2histology_tform{curr_im_idx});
        end

        atlas_slice_aligned = imwarp(im_rigid_transformed, ...
            histology2atlas_tform,'nearest','OutputView', ...
            imref2d(size(slice_atlas(curr_im_idx).av)));

        % % Check match
        % figure; imshowpair(slice_atlas(curr_im_idx).av,atlas_slice_aligned);

        % Add points to volume in CCF space
        curr_ccf_idx = sub2ind(size(tv), ...
            round(slice_atlas_ccf(curr_im_idx).ap(:)), ...
            round(slice_atlas_ccf(curr_im_idx).dv(:)), ...
            round(slice_atlas_ccf(curr_im_idx).ml(:)));

        histology_volume(curr_ccf_idx) = histology_volume(curr_ccf_idx) + ...
            single(atlas_slice_aligned(:));

    end

    % Get max of histology volume in atlas bins
    for curr_atlas_bin = 1:n_atlas_bins
        plot_atlas_ap = round(mean(atlas_bins(curr_atlas_bin+[0,1])));
        histology_atlas_bin_max{curr_animal_idx,curr_atlas_bin} = ...
            permute(max(histology_volume(atlas_bins(curr_atlas_bin): ...
            atlas_bins(curr_atlas_bin+1),:,:),[],1),[2,3,1]);
    end
    
    ap.print_progress_fraction(curr_animal_idx,length(animals));
end

% Plot probe channel overlay by colored animals
animal_colors = ap.colormap('tube',14);
overlay_dilation = 1;
histology_clim = repelem({[200,500]},length(animals),1);

figure; tiledlayout('TileSpacing','none');
for curr_atlas_bin = 1:n_atlas_bins

    % Plot CCF borders from the middle of the bin
    plot_atlas_ap = round(mean(atlas_bins(curr_atlas_bin+[0,1])));
    curr_ccf_borders = imdilate(boundarymask(permute(av(plot_atlas_ap,:,:),[2,3,1])),ones(overlay_dilation));

    % Max, color, and flip contrast
    curr_histology_volume_max_gray = cellfun(@(x,c) ...
        mat2gray(x,c),histology_atlas_bin_max(:,curr_atlas_bin),histology_clim,'uni',false);
        
    curr_histology_volume_max_colored = ...
        cat(4,curr_histology_volume_max_gray{:}).*permute(animal_colors,[3,4,2,1]);
    
    curr_histology_volume_max_colored_white = ...
        curr_histology_volume_max_colored+(1-cat(4,curr_histology_volume_max_gray{:}));

    curr_histology_combined = min(curr_histology_volume_max_colored_white,[],4);

    % Plot CCF over combined colored image
    curr_overlay = imoverlay(curr_histology_combined,curr_ccf_borders,'k');
    nexttile;imagesc(curr_overlay);axis image off;
    drawnow;
end

%% R3 M3: Fraction of quiescent passive trials

%%% Load non-activity data
load_dataset = 'noact';
Marica_2025.figures.load_data;
%%%

% Get which trials in passive were quiescent
animals = unique(bhv.animal,'stable');

passive_quiescent_trials = struct( ...
    'stim_x',cell(length(animals),1), ...
    'quiescent',cell(length(animals),1));

for animal_idx=1:length(animals)
   
    animal = animals{animal_idx};

    % Find passive recording days that also have task
    workflow_passive = {'lcr_passive'};
    recordings_passive = plab.find_recordings(animal, [], workflow_passive);
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);
    training_days = ismember({recordings_passive.day}, {recordings_task.day});
    train_rec_passive = recordings_passive(training_days);
    bhv_days = {train_rec_passive.day};
    ephys_days =  bhv_days([train_rec_passive.ephys]);    

    for use_rec = 1:length(ephys_days)

        rec_day = train_rec_passive(use_rec).day;
        rec_time = train_rec_passive(use_rec).recording{end};
        load_parts.behavior = true;
        ap.load_recording
   
        % Trial stim values
        trial_stim_values = vertcat(trial_events.values.TrialStimX);
        trial_stim_values = trial_stim_values(1:length(stimOn_times));

        % Quiescent trials
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');

        passive_quiescent_trials(animal_idx).stim_x{use_rec,1} = trial_stim_values;
        passive_quiescent_trials(animal_idx).quiescent{use_rec,1} = quiescent_trials;

    end
    ap.print_progress_fraction(animal_idx,length(animals));
end

stim_unique = unique(passive_quiescent_trials(1).stim_x{1});
stim_color = ap.colormap('BKR',length(stim_unique));

% Plot quiescent trials pre/post learning
passive_quiescent_stim = cell2mat(cellfun(@(stim,quiescent) ...
    ap.groupfun(@mean,+quiescent,stim)', ...
    cat(1,passive_quiescent_trials.stim_x), ...
    cat(1,passive_quiescent_trials.quiescent),'uni',false));

[passive_quiescent_stim_ld,passive_quiescent_stim_ld_grp] = ...
    ap.nestgroupfun({@mean,@mean},passive_quiescent_stim, ...
    grp2idx(bhv.animal),bhv.days_from_learning >= 0);
passive_quiescent_stim_ld_sem = ...
    ap.nestgroupfun({@mean,@ap.sem},passive_quiescent_stim, ...
    grp2idx(bhv.animal),bhv.days_from_learning >= 0);

figure;
x_labels = ["Pre-learn","Post-learn"];
errorbar(reordercats(categorical(x_labels),x_labels), ...
    passive_quiescent_stim_ld,passive_quiescent_stim_ld_sem, ...
    'linewidth',2,'capsize',0)
axis padded;
ylabel('Frac. quiescent trials');
set(gca,'ColorOrder',stim_color)
ap.prettyfig;

% Plot quiescent trials binned by trial number
n_trialsplit = length(stim_unique)*5; % n presentations of each stim
passive_trialsplit_idx = cell2mat(cellfun(@(stim) ...
    (floor((0:length(stim)-1)/n_trialsplit)+1)', ...
    cat(1,passive_quiescent_trials.stim_x),'uni',false));

passive_ld_idx = cell2mat(cellfun(@(rec,stim) repelem(rec,length(stim))', ...
    num2cell(bhv.days_from_learning), ...
    cat(1,passive_quiescent_trials.stim_x),'uni',false));

[passive_quiescent_stim_trialsplit,passive_quiescent_stim_trialsplit_grp] = ...
    ap.groupfun(@mean,+cell2mat(cat(1,passive_quiescent_trials.quiescent)), ...
    [passive_ld_idx,passive_trialsplit_idx,cell2mat(cat(1,passive_quiescent_trials.stim_x))]);

passive_quiescent_stim_trialsplit_sem = ...
    ap.groupfun(@ap.sem,+cell2mat(cat(1,passive_quiescent_trials.quiescent)), ...
    [passive_ld_idx,passive_trialsplit_idx,cell2mat(cat(1,passive_quiescent_trials.stim_x))]);

passive_quiescent_stim_trialsplit_grid = ...
    accumarray([grp2idx(passive_quiescent_stim_trialsplit_grp(:,1)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,2)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,3))], ...
    passive_quiescent_stim_trialsplit,[],[],nan);

passive_quiescent_stim_trialsplit_sem_grid = ...
    accumarray([grp2idx(passive_quiescent_stim_trialsplit_grp(:,1)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,2)), ...
    grp2idx(passive_quiescent_stim_trialsplit_grp(:,3))], ...
    passive_quiescent_stim_trialsplit_sem,[],[],nan);

passive_quiescent_stim_daysplit_x = reshape((unique(passive_quiescent_stim_trialsplit_grp(:,1)) + ...
    linspace(0,1,max(passive_quiescent_stim_trialsplit_grp(:,2)+1)))',[],1);

figure;
plot_x_idx = isbetween(passive_quiescent_stim_daysplit_x,-3,3);

passive_quiescent_stim_trialsplit_grid_reshape = ...
    reshape(permute(padarray(passive_quiescent_stim_trialsplit_grid,[0,1],nan,'post'),[2,1,3]),[],3);
passive_quiescent_stim_trialsplit_sem_grid_reshape = ...
    reshape(permute(padarray(passive_quiescent_stim_trialsplit_sem_grid,[0,1],nan,'post'),[2,1,3]),[],3);

ap.errorfill(passive_quiescent_stim_daysplit_x(plot_x_idx), ...
    passive_quiescent_stim_trialsplit_grid_reshape(plot_x_idx,:), ...
    passive_quiescent_stim_trialsplit_sem_grid_reshape(plot_x_idx,:),stim_color);
xline(0,'k');
xlabel('Days from learning');
ylabel('Frac. quiescent trials');
axis padded;
ap.prettyfig;

% ~~~ STATS ~~~
sig_flag = @(p) discretize(p < 0.05,[0,1,Inf],["","*"]);

[passive_quiescent_stim_ld_animal,passive_quiescent_stim_ld_animal_grp] = ...
    ap.groupfun(@mean,passive_quiescent_stim, ...
    [grp2idx(bhv.animal),bhv.days_from_learning >= 0]);

for curr_stim = 1:length(stim_unique)
    stat_p = ranksum( ...
        passive_quiescent_stim_ld_animal(passive_quiescent_stim_ld_animal_grp(:,2) == 0,curr_stim), ...
        passive_quiescent_stim_ld_animal(passive_quiescent_stim_ld_animal_grp(:,2) == 1,curr_stim));

    print_stat('Ranksum pre/post stim %d: p = %.2g%s\n', ...
        stim_unique(curr_stim),stat_p,sig_flag(stat_p));
end

