%% Testing encoding/decoding with toy data

% Load label images
labels_path = fullfile(fileparts(which('AP_av_association_revisions')),'label_images');

% (overlapping)
label_stim = ~imresize(imbinarize(max(imread(fullfile(labels_path,"label_stim.BMP")),[],3)),1/10);
label_move = ~imresize(imbinarize(max(imread(fullfile(labels_path,"label_move.BMP")),[],3)),1/10);

% (overlapping xx)
label_stim = ~imresize(imbinarize(max(imread(fullfile(labels_path,"label_stxx.BMP")),[],3)),1/10);
label_move = ~imresize(imbinarize(max(imread(fullfile(labels_path,"label_moxx.BMP")),[],3)),1/10);

% (stacked)
label_stim_split = ~imresize(imbinarize(max(imread(fullfile(labels_path,"label_stim.BMP")),[],3)),1/20);
label_move_split = ~imresize(imbinarize(max(imread(fullfile(labels_path,"label_move.BMP")),[],3)),1/20);
label_stim = vertcat(label_stim_split,zeros(size(label_move_split)));
label_move = vertcat(zeros(size(label_stim_split)),label_move_split);

% Load example day
animal = 'DS010';
recordings = plab.find_recordings(animal,[],'stim_wheel_right_stage\d');
use_recording = find([recordings.widefield],1,'last');
rec_time = recordings(use_recording).recording{end};
rec_day = recordings(use_recording).day;
verbose = true;
load_parts.widefield = false;
ap.load_recording


% Create toy data
wf_t = widefield_expose_times(1:2:end);
wf_framerate = 1/mean(diff(wf_t));
time_bin_centers = wf_t;
time_bins = [wf_t;wf_t(end)+1/wf_framerate];

stim_event_trace = histcounts(stimOn_times,time_bins);

wheel_starts = timelite.timestamps(diff([0;wheel_move]) == 1);
wheel_event_trace = histcounts(wheel_starts,time_bins);


% (first pass)
event_t = 0.2;
event_samples = round(event_t*wf_framerate);
event_kernel = ones(1,event_samples);
toy_data = reshape(convn( ...
    label_stim.*reshape(stim_event_trace,1,1,[]) + ...
    label_move.*reshape(wheel_event_trace,1,1,[]), ...
    reshape(event_kernel,1,1,[]),'same'),[],length(time_bin_centers));

% (add random noise)
toy_data = toy_data + randn(size(toy_data));

% Decode
t_shifts = 0:event_samples;
lambda = 10;
k_stim = ap.regresskernel(toy_data,stim_event_trace,t_shifts,lambda);
k_move = ap.regresskernel(toy_data,wheel_event_trace,t_shifts,lambda);

ap.imscroll(reshape(k_stim,[size(label_stim),size(k_stim,2)]));axis image;
ap.imscroll(reshape(k_move,[size(label_stim),size(k_stim,2)]));axis image;


% Encode
t_shifts = -event_samples:0;
k_act = ap.regresskernel(cat(1,stim_event_trace,wheel_event_trace),toy_data,t_shifts,0);
ap.imscroll(reshape(permute(k_act(1,:,:),[3,2,1]),[size(label_stim),size(k_stim,2)]));axis image;
ap.imscroll(reshape(permute(k_act(2,:,:),[3,2,1]),[size(label_stim),size(k_stim,2)]));axis image;


% TRYING MOVE BEING VARIABLE

% (make toy data with a few different movements)

% (add stim: constant)
event_t = 0.5;
event_samples = round(event_t*wf_framerate);
event_kernel = ones(1,1,event_samples);
toy_data = convn(label_stim.*reshape(stim_event_trace,1,1,[]),event_kernel,'same');

% (add stim: randomly scaled)
event_t = 0.2;
event_samples = round(event_t*wf_framerate);
event_kernel = ones(1,1,event_samples);
halfnorm_dist = makedist('HalfNormal','mu',0,'sigma',2);
rand_mult_values = random(halfnorm_dist,1,length(time_bin_centers));
toy_data = convn(label_stim.*reshape(stim_event_trace.*rand_mult_values,1,1,[]),event_kernel,'same');

% (add moves: N different types of move)
move_scale = 5;
n_moves = 5;
label_move_different = move_scale*label_move+1*randn([size(label_move),n_moves]);
wheel_starts_grp = ap.shake(ap.quantile_bin(length(wheel_starts),n_moves));
for curr_move = 1:n_moves
    event_samples = curr_move*5;
    event_kernel = ones(1,1,event_samples);

    curr_wheel_event_trace = histcounts(wheel_starts(wheel_starts_grp==curr_move),time_bins);
    toy_data = toy_data + convn(label_move_different(:,:,curr_move).* ...
        reshape(curr_wheel_event_trace,1,1,[]),event_kernel,'same');
end

% (add independent noise)
toy_data = toy_data + 3*randn(size(toy_data));

% (add correlated noise)
toy_data = toy_data + reshape(smooth(randn(1,size(toy_data,3)),100)*10,1,1,[]);

% (reshape)
toy_data = reshape(toy_data,[],length(time_bin_centers));


event_samples = 7;
t_shifts = -round(event_samples*0.5):round(event_samples*1.5);
lambda = 1;

k_stim = ap.regresskernel(toy_data,stim_event_trace,t_shifts,lambda);
ap.imscroll(reshape(k_stim,[size(label_stim),size(k_stim,2)]));axis image;

k_move = ap.regresskernel(toy_data,wheel_event_trace,t_shifts,lambda);
ap.imscroll(reshape(k_move,[size(label_stim),size(k_move,2)]));axis image;



t_shifts = -event_samples:event_samples;
k_act = ap.regresskernel(cat(1,stim_event_trace,wheel_event_trace),toy_data,t_shifts,0);
ap.imscroll(reshape(permute(k_act(1,:,:),[3,2,1]),[size(label_stim),size(k_act,2)]));axis image;
ap.imscroll(reshape(permute(k_act(2,:,:),[3,2,1]),[size(label_stim),size(k_act,2)]));axis image;


% Stim-triggered average
stim_frames = find(stim_event_trace);
stim_frames_pull = stim_frames' + [-10:10];
stim_avg = reshape(permute(nanmean(interp1( ...
    (1:length(time_bin_centers))', ...
    toy_data',stim_frames_pull),1),[3,2,1]), ...
    size(label_stim,1),size(label_stim,2),[]);
ap.imscroll(stim_avg);axis image


% So far it looks like encoding super cleanly pulls responses apart? Maybe
% the problem with real data is that there are a lot of autocorrelations
% etc, while this toy example is random?


%% Comparing encoding/decoding

% Load example day
animal = 'DS010';
recordings = plab.find_recordings(animal,[],'stim_wheel_right_stage\d');
use_recording = find([recordings.widefield],1,'last');
rec_time = recordings(use_recording).recording{end};
rec_day = recordings(use_recording).day;
verbose = true;
ap.load_recording

%%%%%%%%%%%%%% PUPIL PROCESS (adapted from PG code)

% Load pupil SLEAP tracking
pupil_sleap_dir = dir(fullfile(fileparts(mousecam_fn),'**','*.h5'));
pupil_sleap_fn = fullfile(pupil_sleap_dir.folder,pupil_sleap_dir.name);

tracks = h5read(pupil_sleap_fn,'/tracks');       % frames x nodes x 2 (x/y-position)
% (scores aren't used at the moment - PG didn't need previously)
% pointScores = h5read(pupil_sleap_fn,'/point_scores')'; % frames x nodes
% instanceScores = h5read(pupil_sleap_fn, '/instance_scores')'; % transpose to 1 x frames

% Fit circle (solve for a,b,c in x^2 + y^2 + a*x + b*y + c = 0)
% (use only non-NaN vertices, remove points with too few vertices)
circle_fit_fun = @(x,y) [x(~any(isnan([x,y]),2)) y(~any(isnan([x,y]),2)) ...
    ones(sum(~any(isnan([x,y]),2)),1)]\ ...
    -(x(~any(isnan([x,y]),2)).^2+y(~any(isnan([x,y]),2)).^2);

min_pupil_points = 3;
pupil_valid_frames = sum(~any(isnan(tracks),3),2) >= min_pupil_points;

pupil_circle_fit = nan(3,length(mousecam_times));
pupil_circle_fit(:,pupil_valid_frames) = cell2mat(arrayfun(@(frame) ...
    circle_fit_fun(tracks(frame,:,1)',tracks(frame,:,2)'), ...
    find(pupil_valid_frames)','uni',false));

pupil = struct( ...
    'x',-pupil_circle_fit(1,:)./2, ...
    'y',-pupil_circle_fit(2,:)./2, ...
    'diameter',2*sqrt(sum(pupil_circle_fit(1:2,:).^2,1)/4-pupil_circle_fit(3,:)));

% (PG pipeline: z-scored, lowpass, sgolay filt
% pupil_diameterPx_filt = sgolayfilt(lowpass(fillmissing(pupil_diameterZ',"linear"),4,30),3,15);

%%%%%%%%%%%%%%


% Widefield time bins for regressors
time_bins = [wf_t;wf_t(end)+1/wf_framerate];

% Stim onset regressors
stim_regressors = histcounts(stimOn_times,time_bins);

% Move regressors
wheel_starts = timelite.timestamps(diff([0;wheel_move]) == 1);
move_onset_regressors = histcounts(wheel_starts,time_bins);

wheel_velocity_resample = interp1(timelite.timestamps,wheel_velocity,wf_t);
pupil_diameter_resample = interp1(mousecam_times,pupil_radius,wf_t);
pupil_velocity_resample = interp1(sqrt(diff(pupil.x).^2+diff(pupil.y).^2),wf_t);

nanzscore = @(x) (x-nanmean(x))./nanstd(x);

move_regressors = vertcat( ...
    move_onset_regressors,nanzscore(wheel_velocity_resample'), ...
    nanzscore(pupil_diameter_resample'),nanzscore(pupil_velocity_resample'));

% Decoding
n_components = 200;
frame_shifts = -10:30;
lambda = 15;
cv_fold = 5;

skip_t = 60; % seconds start/end to skip for artifacts
skip_frames = round(skip_t*wf_framerate);
[kernels,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
    stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

kernels_px = plab.wf.svd2px(wf_U(:,:,1:size(kernels,1)),kernels);
ap.imscroll(kernels_px,frame_shifts);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image
set(gcf,'name','Decoding');

% Encoding
skip_t = 60; % seconds start/end to skip for artifacts
skip_frames = round(skip_t*wf_framerate);
[kernels,predicted_signals,explained_var] = ...
    ap.regresskernel( ...
    vertcat(stim_regressors(:,skip_frames:end-skip_frames), ...
    move_regressors(:,skip_frames:end-skip_frames)), ...
    wf_V(1:n_components,skip_frames:end-skip_frames), ...
    frame_shifts,lambda,[],cv_fold);

kernels_px = plab.wf.svd2px(wf_U(:,:,1:size(kernels,3)),permute(kernels,[3,2,1]));
ap.imscroll(kernels_px,frame_shifts);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image
set(gcf,'name','Encoding');

% REVIEWER SUGGESTION: regress out movement, decoding on residual

skip_t = 60; % seconds start/end to skip for artifacts
skip_frames = round(skip_t*wf_framerate);
[kernels,predicted_signals_move] = ...
    ap.regresskernel( ...
    move_regressors(:,skip_frames:end-skip_frames), ...
    wf_V(1:n_components,skip_frames:end-skip_frames), ...
    frame_shifts,lambda,[],cv_fold);

[kernels,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames) - ...
    predicted_signals_move, ...
    stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

kernels_px = plab.wf.svd2px(wf_U(:,:,1:size(kernels,1)),kernels);
ap.imscroll(kernels_px,frame_shifts);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image
set(gcf,'name','Decoding (move residual)');

% (alternate - remove stim and decode move)
[move_decode_kernels,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
    move_onset_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

[kernels,predicted_signals_stim] = ...
    ap.regresskernel( ...
    stim_regressors(:,skip_frames:end-skip_frames), ...
    wf_V(1:n_components,skip_frames:end-skip_frames), ...
    frame_shifts,lambda,[],cv_fold);

[move_decode_kernels_nomove,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames) - ...
    predicted_signals_stim, ...
    move_onset_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

kernels_px = cat(2, ...
    plab.wf.svd2px(wf_U(:,:,1:n_components),move_decode_kernels), ...
    plab.wf.svd2px(wf_U(:,:,1:n_components),move_decode_kernels_nomove));
ap.imscroll(kernels_px,frame_shifts);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image
set(gcf,'name','Decoding (move residual)');

% (alternate - remove move and decode move)
[move_decode_kernels,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
    move_onset_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

[kernels,predicted_signals_stim] = ...
    ap.regresskernel( ...
    move_onset_regressors(:,skip_frames:end-skip_frames), ...
    wf_V(1:n_components,skip_frames:end-skip_frames), ...
    frame_shifts,lambda,[],cv_fold);

[move_decode_kernels_nomove,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames) - ...
    predicted_signals_stim, ...
    move_onset_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

kernels_px = cat(2, ...
    plab.wf.svd2px(wf_U(:,:,1:n_components),move_decode_kernels), ...
    plab.wf.svd2px(wf_U(:,:,1:n_components),move_decode_kernels_nomove));
ap.imscroll(kernels_px,frame_shifts);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image
set(gcf,'name','Decoding (move residual)');


% (sanity check: remove stim, decode stim)
[stim_decode_kernels,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
    stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

[kernels,predicted_signals_stim] = ...
    ap.regresskernel( ...
    stim_regressors(:,skip_frames:end-skip_frames), ...
    wf_V(1:n_components,skip_frames:end-skip_frames), ...
    frame_shifts,lambda,[],cv_fold);

[stim_decode_kernels_nostim,predicted_signals,explained_var] = ...
    ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames) - ...
    predicted_signals_stim, ...
    stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

kernels_px = cat(2, ...
    plab.wf.svd2px(wf_U(:,:,1:n_components),stim_decode_kernels), ...
    plab.wf.svd2px(wf_U(:,:,1:n_components),stim_decode_kernels_nostim));
ap.imscroll(kernels_px,frame_shifts);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image
set(gcf,'name','Decoding (move residual)');


% -> this all seems to work, do in batch for all post-learn?


%% Decoding after encoding-regressing out movement

% Load behavior
% (individual animals packaged here:)
% https://github.com/Soda1212-0808/DS_PetersLab/blob/main/Data_exploration/projects/cross_modal_PFC_striatum_2025/ds_behavior_single_mouse.m
% (across-animals packaged here:)
% Song_2025.package.save_beahvior
data_path = fullfile(plab.locations.server_path,'Lab','Papers','Song_2025','data');
load(fullfile(data_path,'behavior.mat'));

% Get animals from behavior structure
animals = cell(2,1);
animals{1} = behavior_aligned{1}.reaction_time.name; % VA
animals{2} = behavior_aligned{2}.reaction_time.name; % AV

% Loop through modalities/animal groups/animals, get kernels
animal_kernels = struct;
for curr_animal_group = 1:length(animals)
    for curr_animal_idx = 1:length(animals{curr_animal_group})

        animal = animals{curr_animal_group}{curr_animal_idx};

        % Get corresponding learned days from behavior
        % (copied code from Song_2025.package.save_beahvior to get day match)
        mouse_data_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Da_Song\Data_back_up\project_cross_model\wf_data';
        raw_data_behavior = load(fullfile(mouse_data_path,'behavior',[animal,'_behavior.mat']));
        raw_data_task = load(fullfile(mouse_data_path,'task',[animal,'_task.mat']));
        if strcmp(animal,'AP019')
            use_rec_idx =  (1:15)';
        else
            [~,use_rec_idx] = ismember(raw_data_task.workflow_day,raw_data_behavior.workflow_day);
        end

        rec_days = raw_data_behavior.workflow_day(use_rec_idx);
        rec_days_unimodal_idx = find(~contains(behavior_each_mice{curr_animal_group}.workflow_name{curr_animal_idx},'mixed'));

        % Loop though animals, run encoding/decoding
        for curr_recording = rec_days_unimodal_idx'

            % All in a try/catch loop: skip on errors
            % (at least one day is too long = regression memory error)
            try

            % Get day and recording
            rec_day = rec_days{curr_recording};
            recordings = plab.find_recordings(animal,rec_day,'stim_wheel*');
            % (copied from ds_behavior_single_mouse: use longest rec)
            if length(recordings.index)>1
                for mm=1:length(recordings.index)
                    rec_time = recordings.recording{mm};
                    timelite_fn = plab.locations.filename('server',animal,rec_day,rec_time,'timelite.mat');
                    timelite = load(timelite_fn);
                    time(mm)=length(timelite.timestamps);
                end
                [~,index_real]=max(time);
            else
                index_real=1;
            end

            rec_time = recordings.recording{index_real};

            preload_vars = who;
            load_parts.widefield = true;
            load_parts.widefield_master = true;
            load_parts.mousecam = true;
            ap.load_recording;

            % Load and prep pupil SLEAP tracking
            pupil_sleap_dir = dir(fullfile(fileparts(mousecam_fn),'**','*.h5'));
            if isempty(pupil_sleap_dir)
                continue
            end
            pupil_sleap_fn = fullfile(pupil_sleap_dir.folder,pupil_sleap_dir.name);

            tracks = h5read(pupil_sleap_fn,'/tracks');       % frames x nodes x 2 (x/y-position)
            % (scores aren't used at the moment - PG didn't need previously)
            % pointScores = h5read(pupil_sleap_fn,'/point_scores')'; % frames x nodes
            % instanceScores = h5read(pupil_sleap_fn, '/instance_scores')'; % transpose to 1 x frames

            % Fit circle (solve for a,b,c in x^2 + y^2 + a*x + b*y + c = 0)
            % (use only non-NaN vertices, remove points with too few vertices)
            circle_fit_fun = @(x,y) [x(~any(isnan([x,y]),2)) y(~any(isnan([x,y]),2)) ...
                ones(sum(~any(isnan([x,y]),2)),1)]\ ...
                -(x(~any(isnan([x,y]),2)).^2+y(~any(isnan([x,y]),2)).^2);

            min_pupil_points = 3;
            pupil_valid_frames = sum(~any(isnan(tracks),3),2) >= min_pupil_points;

            pupil_circle_fit = nan(3,length(mousecam_times));
            pupil_circle_fit(:,pupil_valid_frames) = cell2mat(arrayfun(@(frame) ...
                circle_fit_fun(tracks(frame,:,1)',tracks(frame,:,2)'), ...
                find(pupil_valid_frames)','uni',false));

            pupil = struct( ...
                'x',-pupil_circle_fit(1,:)./2, ...
                'y',-pupil_circle_fit(2,:)./2, ...
                'diameter',2*sqrt(sum(pupil_circle_fit(1:2,:).^2,1)/4-pupil_circle_fit(3,:)));

            % Regression parameters
            time_bins = [wf_t;wf_t(end)+1/wf_framerate];
            n_components = 200;
            frame_shifts = -10:30;
            lamda_encode = 0;
            lambda_decode = 15;
            cv_fold = 5;
            skip_t = 60; % seconds start/end to skip for artifacts
            skip_frames = round(skip_t*wf_framerate);

            % Stim onset regressors
            stim_regressors = histcounts(stimOn_times,time_bins);

            % Move regressors
            % (move onset, wheel velocity, pupil diameter, pupil velocity)
            wheel_starts = timelite.timestamps(diff([0;wheel_move]) == 1);
            move_onset_regressors = histcounts(wheel_starts,time_bins);

            wheel_velocity_resample = interp1(timelite.timestamps,wheel_velocity,wf_t);
            pupil_diameter_resample = interp1(mousecam_times,pupil.diameter,wf_t);
            pupil_velocity_resample = interp1(sqrt(diff(pupil.x).^2+diff(pupil.y).^2),wf_t);

            move_regressors = vertcat( ...
                move_onset_regressors,wheel_velocity_resample', ...
                pupil_diameter_resample',pupil_velocity_resample');

            % Encoding: stim kernel
            kernels_stimmove_encode = ...
                ap.regresskernel( ...
                vertcat(stim_regressors(:,skip_frames:end-skip_frames), ...
                move_regressors(:,skip_frames:end-skip_frames)), ...
                wf_V(1:n_components,skip_frames:end-skip_frames), ...
                frame_shifts,lamda_encode,[],cv_fold);

            % Encoding: stim and move separately (to use residuals)
            [kernels_stim_encode,predicted_signals_stim] = ...
                ap.regresskernel(stim_regressors(:,skip_frames:end-skip_frames), ...
                wf_V(1:n_components,skip_frames:end-skip_frames), ...
                frame_shifts,lamda_encode,[],cv_fold);

            [kernels_move_encode,predicted_signals_move] = ...
                ap.regresskernel(move_regressors(:,skip_frames:end-skip_frames), ...
                wf_V(1:n_components,skip_frames:end-skip_frames), ...
                frame_shifts,lamda_encode,[],cv_fold);

            % Decoding: regress to stim and move (on full and residuals)
            kernels_stim_decode_full = ...
                ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames), ...
                stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda_decode,[],cv_fold);

            kernels_stim_decode_moveresiduals = ...
                ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames) - ...
                predicted_signals_move, ...
                stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda_decode,[],cv_fold);

            kernels_stim_decode_stimresiduals = ...
                ap.regresskernel(wf_V(1:n_components,skip_frames:end-skip_frames) - ...
                predicted_signals_stim, ...
                stim_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda_decode,[],cv_fold);

            % Store variables
            animal_kernels(curr_animal_group).stimmove_encode{curr_animal_idx}{curr_recording} = permute(kernels_stimmove_encode,[3,2,1]);
            animal_kernels(curr_animal_group).stim_encode{curr_animal_idx}{curr_recording} = permute(kernels_stim_encode,[3,2,1]);
            animal_kernels(curr_animal_group).move_encode{curr_animal_idx}{curr_recording} = permute(kernels_move_encode,[3,2,1]);
            animal_kernels(curr_animal_group).stim_decode_full{curr_animal_idx}{curr_recording} = kernels_stim_decode_full;
            animal_kernels(curr_animal_group).stim_decode_moveresiduals{curr_animal_idx}{curr_recording} =  kernels_stim_decode_moveresiduals;
            animal_kernels(curr_animal_group).stim_decode_stimresiduals{curr_animal_idx}{curr_recording} = kernels_stim_decode_stimresiduals;

            % Clear recording variables and print progress
            clearvars('-except',preload_vars{:});
            fprintf('Finished %s day %d/%d\n',animal,curr_recording,length(rec_days))
            reset(gpuDevice());

            catch me
                fprintf('Error: "%s", skipping: %s day %d/%d\n',me.message,animal,curr_recording,length(rec_days))
                warning(me.message);
                clearvars('-except',preload_vars{:});
                continue
            end
        end
    end
end

save_path = "C:\Users\petersa\Documents\PetersLab\papers\2025_Song\data_test";
save_fn = fullfile(save_path,'encoding_decoding_kernels');
save(save_fn,'animal_kernels');
fprintf('\n---\nSaved %s\n---\n',save_fn)

% (for test plotting)
wf_U = plab.wf.load_master_U;
n_components = 200;

kernels_px = plab.wf.svd2px(wf_U(:,:,1:n_components),x);
ap.imscroll(kernels_px);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image
set(gcf,'name','Encoding');

% %%%%%% UNUSED AT THE MOMENT - USE ON FIG GENERATION?
% recordings_learned = behavior_each_mice{curr_animal_group}.learned{curr_animal_idx}(bhv_curr_task_idx);












