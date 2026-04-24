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


time_bins = [wf_t;wf_t(end)+1/wf_framerate];
stim_regressors = histcounts(stimOn_times,time_bins);

wheel_starts = timelite.timestamps(diff([0;wheel_move]) == 1);
move_regressors = histcounts(wheel_starts,time_bins);


% Decoding
n_components = 200;

frame_shifts = -10:30;
lambda = 20;
cv_fold = 1;

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

% (just checking: try alternate - remove stim and decode move?
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
    move_regressors(:,skip_frames:end-skip_frames),-frame_shifts,lambda,[],cv_fold);

kernels_px = plab.wf.svd2px(wf_U(:,:,1:size(kernels,1)),kernels);
ap.imscroll(kernels_px,frame_shifts);
clim(max(abs(clim)).*[-1,1]);
colormap(ap.colormap('PWG'));
axis image
set(gcf,'name','Decoding (move residual)');



