% wf_retinotopy: run after loading sparse noise experiment with
% ap.load_experiment

% Use raw blue signal (best SNR)
surround_window = [0.3,0.5]; % 6s = [0.3,0.5], deconv = [0.05,0.15]
framerate = 1./nanmean(diff(wf_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
response_n = nan(n_y_squares,n_x_squares);
response_grid = cell(n_y_squares,n_x_squares);
for px_y = 1:n_y_squares
    for px_x = 1:n_x_squares

        align_times = ...
            stim_times(find( ...
            (noise_locations(px_y,px_x,1:end-1) == 128 & ...
            noise_locations(px_y,px_x,2:end) == 255) | ...
            (noise_locations(px_y,px_x,1:end-1) == 128 & ...
            noise_locations(px_y,px_x,2:end) == 0))+1);

        response_n(px_y,px_x) = length(align_times);

        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < wf_t(2) | ...
            align_times + surround_time(2) > wf_t(end)) = [];

        % Get stim-aligned responses, 2 choices:

        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]);

        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = align_times + surround_time;
        frame_edges = [wf_t;wf_t(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);

        % Get stim-aligned baseline (at stim onset)
        align_baseline_times = align_times;
        align_frames_baseline = discretize(align_baseline_times,frame_edges);

        % Don't use NaN frames (delete, dirty)
        nan_stim = any(isnan(align_frames),2) | isnan(align_frames_baseline);
        align_frames(nan_stim,:) = [];
        align_frames_baseline(nan_stim,:) = [];

        % Define the peri-stim V's as subtracting first frame (baseline)
        peri_stim_v = ...
            reshape(wf_V_raw{1}(:,align_frames)',size(align_frames,1),size(align_frames,2),[]) - ...
            nanmean(reshape(wf_V_raw{1}(:,align_frames_baseline)',size(align_frames_baseline,1),size(align_frames_baseline,2),[]),2);

        mean_peri_stim_v = permute(mean(peri_stim_v,2),[3,1,2]);

        % Store V's
        response_grid{px_y,px_x} = mean_peri_stim_v;

    end
end

% Get position preference for every pixel
U_downsample_factor = 1; %2 if max method
screen_resize_scale = 1; %3 if max method
filter_sigma = (screen_resize_scale*2);

% Downsample U
[Uy,Ux,nSV] = size(wf_U);
Ud = imresize(wf_U,1/U_downsample_factor,'bilinear');

% Convert V responses to pixel responses
use_svs = 1:min(2000,size(Ud,3)); % can de-noise if reduced
n_boot = 10;

response_mean_bootstrap = cellfun(@(x) bootstrp(n_boot,@mean,x')',response_grid,'uni',false);

%%%%%%% Get retinotopy (for each bootstrap)
vfs_boot = nan(size(Ud,1),size(Ud,2),n_boot);
for curr_boot = 1:n_boot

    response_mean = cell2mat(cellfun(@(x) x(:,curr_boot),response_mean_bootstrap(:),'uni',false)');
    stim_im_px = reshape(permute(plab.wf.svd2px(Ud(:,:,use_svs), ...
        response_mean(use_svs,:)),[3,1,2]),n_y_squares,n_x_squares,[]);
    gauss_filt = fspecial('gaussian',[n_y_squares,n_x_squares],filter_sigma);
    stim_im_smoothed = imfilter(imresize(stim_im_px,screen_resize_scale,'bilinear'),gauss_filt);

    %     % (for troubleshooting: scroll through px stim responses)
    %     AP_imscroll(reshape(permute(stim_im_px,[3,1,2]),size(Ud,1),size(Ud,2),[]));
    %     axis image;clim(max(abs(clim)).*[-1,1]);colormap(AP_colormap('BWR'));

    % Get center-of-mass screen response for each widefield pixel
    [yy,xx] = ndgrid(1:size(stim_im_smoothed,1),1:size(stim_im_smoothed,2));
    m_xr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,xx),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));
    m_yr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,yy),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));

    % Calculate and plot sign map (dot product between horz & vert gradient)

    % 1) get gradient direction
    [~,Vdir] = imgradient(imgaussfilt(m_yr,1));
    [~,Hdir] = imgradient(imgaussfilt(m_xr,1));

    % 3) get sin(difference in direction) if retinotopic, H/V should be
    % orthogonal, so the closer the orthogonal the better (and get sign)
    angle_diff = sind(Vdir-Hdir);
    angle_diff(isnan(angle_diff)) = 0;

    vfs_boot(:,:,curr_boot) = angle_diff;
end

vfs = imgaussfilt(nanmean(vfs_boot,3),2);

figure;
imagesc(vfs)
axis image off;
colormap(AP_colormap('BWR'));
title(sprintf('%s, %s',animal,rec_day));


%% Align retinotopy to master retinotopy
% (TURNING THIS OFF)
% % Plot CCF areas and coordinates aligned to master retinotopy
% 
% % Align retinotopy to CCF 
% [vfs_ccf_aligned,im_tform] = ap.wf_align(vfs,[],[],'new_animal');
% close(gcf);
% 
% % Apply alignment to average image
% wf_avg_aligned = imwarp(wf_avg,im_tform,'Outputview',imref2d(size(vfs_ccf_aligned)));
% 
% figure('color','w');
% tiledlayout('flow','TileSpacing','tight');
% 
% % Plot retinotopy with CCF overlay
% h = nexttile;
% imagesc(vfs_ccf_aligned);
% axis image off;
% colormap(h,AP_colormap('BWR'));
% ap.wf_draw('ccf','k');
% 
% % Plot average image with retinotopy/CCF/grid overlay
% h = nexttile;
% imagesc(wf_avg_aligned);
% axis image off;
% colormap(h,'gray');
% ap.wf_draw('ccf','b');
% % ap.wf_draw('grid_aligned','y');


%% Plot average/VFS overlay
% (TURNING THIS OFF)
% figure;
% ax_avg = axes;
% ax_vfs = axes;
% h1 = imagesc(ax_avg,wf_avg_aligned);
% colormap(ax_avg,gray);
% clim(ax_avg,prctile(wf_avg(:),[20,99]));
% h2 = imagesc(ax_vfs,vfs_ccf_aligned);
% colormap(ax_vfs,AP_colormap('BWR'));
% clim(ax_vfs,[-1,1]);
% set(ax_avg,'Visible','off');
% axes(ax_avg); axis image off;
% set(ax_vfs,'Visible','off');
% axes(ax_vfs); axis image off;
% set(h2,'AlphaData',mat2gray(abs(vfs_ccf_aligned))*0.5);
% colormap(ax_avg,gray);
% linkaxes([ax_avg,ax_vfs]);


