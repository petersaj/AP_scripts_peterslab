%% Exploratory ephys analysis


%% Ephys raster

if contains(bonsai_workflow,'lcr')
    % (L/C/R passive)
    align_times_all = stimOn_times;
    align_category_all = vertcat(trial_events.values.TrialStimX);
    % (get only quiescent trials)
    framerate = 30;
    wheel_window = [0,0.5];
    wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
    wheel_window_t_peri_event = align_times_all + wheel_window_t;
    event_aligned_move = interp1(timelite.timestamps, ...
        +wheel_move,wheel_window_t_peri_event,'previous');
    quiescent_trials = ~any(event_aligned_move,2);
    align_times = align_times_all(quiescent_trials);
    align_category = align_category_all(quiescent_trials);

    AP_cellraster(align_times,align_category);

elseif contains(bonsai_workflow,'stim_wheel')
    % (task)
    stimOn_times = photodiode_times(1:2:end);
    wheel_starts = timelite.timestamps(diff([0;wheel_move]) == 1);

    % (get wheel starts when no stim on screen: not sure this works yet)
    wheel_starts_iti = ...
        wheel_starts(interp1(photodiode_times, ...
        photodiode_values,wheel_starts,'previous') == 0);

    [~,rxn_sort_idx] = sort(stim_to_move);
    AP_cellraster({stimOn_times,wheel_move_time,wheel_starts_iti,reward_times}, ...
        {rxn_sort_idx,rxn_sort_idx,1:length(wheel_starts_iti),1:length(reward_times)});
end


%% Ephys sparse noise? 

stim_screen = discretize(noise_locations,[0,128,129,255],-1:1);
% Get stim times vector (x,y)
nY = size(stim_screen,1);
nX = size(stim_screen,2);
stim_times_grid = cell(nY,nX);
for x = 1:nX
    for y = 1:nY      
        align_stims = (stim_screen(y,x,1:end-1) == 0) & ...
            (stim_screen(y,x,2:end) == 0);
        align_times = stim_times(find(align_stims)+1);
        stim_times_grid{y,x} = align_times;
    end
end

% Vectorize stim times by y,x position
[stim_x,stim_y] = meshgrid(1:nX,1:nY);
stim_positions = cellfun(@(x,y,times) repmat([y,x],length(times),1), ...
    num2cell(stim_x),num2cell(stim_y),stim_times_grid,'uni',false);

% Pick unit
% use_spikes = spike_times_timeline(spike_templates == 276);
use_spikes = spike_times_timeline(spike_templates == 229);
% use_spikes = spike_times_timeline(ismember(spike_templates, ...
%     find(template_depths > 1700 & template_depths < 2300)));

% Get stim times vector (x,y)
stim_aligned_avg = cell(nY,nX);
raster_window = [-0.05,0.2];
smooth_size = 2;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
for x = 1:nX
    for y = 1:nY      
        
        psth_bin_size = 0.01;
        [psth,bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
            use_spikes,stim_times_grid{y,x}, ...
            raster_window, psth_bin_size);
        
        psth_smooth = conv2(psth,smWin,'same');
        stim_aligned_avg{y,x} = psth_smooth;
        
    end
end
stim_aligned_avg_cat = cell2mat(cellfun(@(x) permute(x,[1,3,2]),stim_aligned_avg,'uni',false));
AP_imscroll(stim_aligned_avg_cat,bins);
axis image;


%% Widefield/ephys regression maps

% Set upsample value for regression
upsample_factor = 1;
sample_rate = (1/mean(diff(wf_times)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = wf_times(find(wf_times > skip_seconds,1)):1/sample_rate:wf_times(find(wf_times-wf_times(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% (to group multiunit by evenly spaced depths)
n_depths = 6;
depth_group_edges = round(linspace(0,4000,n_depths+1));
[depth_group_n,~,depth_group] = histcounts(spike_depths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

% % (for manual depth)
% depth_group_edges = [1000,4000];
% n_depths = length(depth_group_edges) - 1;
% [depth_group_n,depth_group] = histcounts(spike_depths,depth_group_edges);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
%     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
%         ismember(spike_templates,find(msn)));
    
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
    
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

use_svs = 1:200;
kernel_t = [-0.5,0.5];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 10;
zs = [false,false];
cvfold = 5;
return_constant = false;
use_constant = true;
    
fVdf_deconv_resample = interp1(wf_times,wf_V(use_svs,:)',time_bin_centers)';
        
% TO USE DECONV
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample, ...
    binned_spikes_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);

% Convert kernel to pixel space
r_px = plab.wf.svd2px(wf_U(:,:,use_svs),k);

AP_imscroll(r_px,kernel_frames/wf_framerate);
caxis([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
colormap(AP_colormap('BWR'));
axis image;

% Get center of mass for each pixel
% (get max r for each pixel, filter out big ones)
r_px_max = squeeze(max(r_px,[],3));
r_px_max(isnan(r_px_max)) = 0;
% for i = 1:n_depths
%     r_px_max(:,:,i) = medfilt2(r_px_max(:,:,i),[10,10]);
% end
r_px_max_norm = bsxfun(@rdivide,r_px_max, ...
    permute(max(reshape(r_px_max,[],n_depths),[],1),[1,3,2]));
r_px_max_norm(isnan(r_px_max_norm)) = 0;
r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depths,[1,3,2])),3)./sum(r_px_max_norm,3);

% Plot map of cortical pixel by preferred depth of probe
r_px_com_col = ind2rgb(round(mat2gray(r_px_com,[1,n_depths])*255),jet(255));
figure;
a1 = axes('YDir','reverse');
imagesc(wf_avg); colormap(gray); caxis([0,prctile(wf_avg(:),99.7)]);
axis off; axis image;
a2 = axes('Visible','off'); 
p = imagesc(r_px_com_col);
axis off; axis image;
set(p,'AlphaData',mat2gray(max(r_px_max_norm,[],3), ...
     [0,double(prctile(reshape(max(r_px_max_norm,[],3),[],1),95))]));
set(gcf,'color','w');

c1 = colorbar('peer',a1,'Visible','off');
c2 = colorbar('peer',a2);
ylabel(c2,'Depth (\mum)');
colormap(c2,jet);
set(c2,'YDir','reverse');
set(c2,'YTick',linspace(0,1,n_depths));
set(c2,'YTickLabel',round(linspace(depth_group_edges(1),depth_group_edges(end),n_depths)));












