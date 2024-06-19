% %% Save data (just for AP to package - don't use)
% 
% animal = 'AP014';
% rec_day = '2024-01-19';
% load_parts.widefield = true;
% load_parts.widefield_align = false;
% ap.load_recording;
% 
% data_dir = fullfile('D:\retinotopy_for_Jennifer',sprintf('%s_%s',animal,rec_day));
% 
% save(fullfile(data_dir,'noise_locations'),'noise_locations');
% save(fullfile(data_dir,'stim_times'),'stim_times');
% 
% widefield_spatial_components = wf_U_raw{1};
% save(fullfile(data_dir,'widefield_spatial_components'),'widefield_spatial_components');
% 
% widefield_temporal_components = wf_V_raw{1};
% save(fullfile(data_dir,'widefield_temporal_components'),'widefield_temporal_components');
% 
% widefield_times = wf_t_all{1};
% save(fullfile(data_dir,'widefield_times'),'widefield_times');
% 
% widefield_average = wf_avg;
% save(fullfile(data_dir,'widefield_average'),'widefield_average');


%% Get visual field sign retinotopy

% Dataset 1: 
% data_dir = 'D:\retinotopy_for_Jennifer\AP010_2023-08-09';

% Dataset 2:
data_dir = 'D:\retinotopy_for_Jennifer\AP014_2024-01-19';

% Load data
% - sparse noise display (vertical location x horizontal location x stim number)
load(fullfile(data_dir,'noise_locations.mat'));
% - time of each stimulus in noise_locations
load(fullfile(data_dir,'stim_times.mat'));
% - widefield spatial components (U from SVD)
load(fullfile(data_dir,'widefield_spatial_components.mat'));
% - widefield temporal components (V from SVD)
load(fullfile(data_dir,'widefield_temporal_components.mat'));
% - times of each widefield frame
load(fullfile(data_dir,'widefield_times.mat'));
% - average raw widefield image
load(fullfile(data_dir,'widefield_average.mat'));

% Define window around stimulus to get activity
surround_window = [0.3,0.5]; % seconds from stimulus onset
framerate = 1./nanmean(diff(widefield_times));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);

% Get average activity for each square on the screen
[n_y_squares,n_x_squares] = size(noise_locations,[1,2]);

response_n = nan(n_y_squares,n_x_squares);
response_grid = cell(n_y_squares,n_x_squares);
for px_y = 1:n_y_squares
    for px_x = 1:n_x_squares

        % (use both gray-to-white and gray-to-black flips)
        align_times = ...
            stim_times(find( ...
            (noise_locations(px_y,px_x,1:end-1) == 128 & ...
            noise_locations(px_y,px_x,2:end) == 255) | ...
            (noise_locations(px_y,px_x,1:end-1) == 128 & ...
            noise_locations(px_y,px_x,2:end) == 0))+1);

        response_n(px_y,px_x) = length(align_times);

        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < widefield_times(2) | ...
            align_times + surround_time(2) > widefield_times(end)) = [];

        % Get frames around each stim presentation
        align_surround_times = align_times + surround_time;
        frame_edges = [widefield_times;widefield_times(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);

        % Get stim-aligned baseline (= stim onset frame, before activity)
        align_baseline_times = align_times;
        align_frames_baseline = discretize(align_baseline_times,frame_edges);

        % Grab temporal components (V's) around stim, minus baseline
        peri_stim_v = ...
            reshape(widefield_temporal_components(:,align_frames)',size(align_frames,1),size(align_frames,2),[]) - ...
            nanmean(reshape(widefield_temporal_components(:,align_frames_baseline)',size(align_frames_baseline,1),size(align_frames_baseline,2),[]),2);

        % Get time-average temporal component for each stim
        mean_peri_stim_v = permute(mean(peri_stim_v,2),[3,1,2]);

        % Store average temporal components
        response_grid{px_y,px_x} = mean_peri_stim_v;

    end
end

% Get bootstrapped average response for each pixel
% (bootstrapping helps reduce noise)
n_boot = 1;
response_mean_bootstrap = cellfun(@(x) bootstrp(n_boot,@mean,x')',response_grid,'uni',false);

% Loop through bootstrapped averages, get visual field sign map
filter_sigma = 2; % (gaussian smooth standard deviation size in pixels)
vfs_boot = nan(size(widefield_spatial_components,1),size(widefield_spatial_components,2),n_boot);
for curr_boot = 1:n_boot

    % Get current bootstrap mean
    response_mean = cell2mat(cellfun(@(x) x(:,curr_boot),response_mean_bootstrap(:),'uni',false)');

    % Convert SVD components to pixels (pixels = U*V)
    im_px = reshape( ...
        reshape(widefield_spatial_components,[],size(widefield_spatial_components,3)) * ...
        reshape(response_mean,size(response_mean,1),[]), ...
        [size(widefield_spatial_components,[1:2]),size(response_mean,2)]);

    % Reshape to be vertical square x horizontal square x pixel
    stim_im_px = reshape(permute(im_px,[3,1,2]),n_y_squares,n_x_squares,[]);
    gauss_filt = fspecial('gaussian',[n_y_squares,n_x_squares],filter_sigma);
    stim_im_smoothed = imfilter(stim_im_px,gauss_filt);

    % Get center-of-mass screen response for each widefield pixel
    [yy,xx] = ndgrid(1:size(stim_im_smoothed,1),1:size(stim_im_smoothed,2));
    m_xr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,xx),1),2)./ ...
        sum(sum(stim_im_smoothed.^2,1),2),size(widefield_spatial_components,1),size(widefield_spatial_components,2));
    m_yr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,yy),1),2)./ ...
        sum(sum(stim_im_smoothed.^2,1),2),size(widefield_spatial_components,1),size(widefield_spatial_components,2));

    % Calculate and plot visual field sign sign

    % 1) Get gradient direction
    [~,Vdir] = imgradient(imgaussfilt(m_yr,1));
    [~,Hdir] = imgradient(imgaussfilt(m_xr,1));

    % 3) Get sin(difference in direction) to measure gradient
    % orthogonality. (This gives degree of retinotopy: horizontal and
    % vertical are orthogonal in screen space, so if they're also
    % orthogonal on the brain, that region is retinotopic. The range is
    % [-1:1], with 0 = not retinotopic, and sign indicating whether an area
    % is the same or mirrored orientation relative to the screen).
    angle_diff = sind(Vdir-Hdir);
    angle_diff(isnan(angle_diff)) = 0;

    % Store the visual field sign for this bootstrap
    vfs_boot(:,:,curr_boot) = angle_diff;
end

% Get the average visual field sign across bootstraps, and gaussian smooth
vfs = imgaussfilt(nanmean(vfs_boot,3),2);

% Plot visual field sign by itself
figure;
imagesc(vfs)
axis image off;
colormap(parula);

% Plot visual field sign overlaid on average image
figure;
ax_avg = axes;
ax_vfs = axes;
h1 = imagesc(ax_avg,widefield_average);
colormap(ax_avg,gray);
clim(ax_avg,prctile(widefield_average(:),[20,99]));
h2 = imagesc(ax_vfs,vfs);
colormap(ax_vfs,parula);
clim(ax_vfs,[-1,1]);
set(ax_avg,'Visible','off');
axes(ax_avg); axis image off;
set(ax_vfs,'Visible','off');
axes(ax_vfs); axis image off;
set(h2,'AlphaData',mat2gray(abs(vfs))*0.5);
colormap(ax_avg,gray);
linkaxes([ax_avg,ax_vfs]);






