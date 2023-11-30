%% Exploratory ephys analysis


%% Cell raster

if contains(bonsai_workflow,'lcr')
    % (L/C/R passive)

    align_times_all = stimOn_times;
    align_category_all = vertcat(trial_events.values.TrialStimX);
    % (get only quiescent trials)
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        1:length(stimOn_times));
    align_times = align_times_all(quiescent_trials);
    align_category = align_category_all(quiescent_trials);

    AP_cellraster(align_times,align_category);

elseif contains(bonsai_workflow,'stim_wheel')
    % (task)
    wheel_starts = timelite.timestamps(diff([0;wheel_move]) == 1);

    % (get wheel starts when no stim on screen: not sure this works yet)
    wheel_starts_iti = ...
        wheel_starts(interp1(photodiode_times, ...
        photodiode_values,wheel_starts,'previous') == 0);

    [~,rxn_sort_idx] = sort(stim_to_move);
    AP_cellraster({stimOn_times(1:n_trials),stim_move_time,wheel_starts_iti,reward_times}, ...
        {rxn_sort_idx,rxn_sort_idx,1:length(wheel_starts_iti),1:length(reward_times)});
end

%% MUA PSTH by depth 

% % (evenly spaced depths)
% n_depths = 20;
% depth_group_edges = round(linspace(0,3840,n_depths+1))';
% [depth_group_n,~,depth_group] = histcounts(spike_depths,depth_group_edges);
% depth_groups_used = unique(depth_group);
% depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

% (for clickable manual depths)
figure('units','normalized','position',[0.05,0.2,0.1,0.7]); 
unit_axes = axes('YDir','reverse'); hold on; xlabel('Norm spike rate');ylabel('Depth');

if exist('probe_areas','var')
    probe_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
        probe_areas{1}.color_hex_triplet,'uni',false)),[1,3,2]);

    probe_areas_boundaries = probe_areas{1}.probe_depth;
    probe_areas_centers = mean(probe_areas_boundaries,2);

    probe_areas_image_depth = 0:1:max(probe_areas_boundaries,[],'all');
    probe_areas_image_idx = interp1(probe_areas_boundaries(:,1), ...
        1:height(probe_areas{1}),probe_areas_image_depth, ...
        'previous','extrap');
    probe_areas_image = probe_areas_rgb(probe_areas_image_idx,:,:);

    image(unit_axes,[0,1],probe_areas_image_depth,probe_areas_image);
    yline(unique(probe_areas_boundaries(:)),'color','k','linewidth',1);
    set(unit_axes,'YTick',probe_areas_centers,'YTickLabels',probe_areas{1}.acronym);
end

norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
unit_dots = scatter3( ...
    norm_spike_n,template_depths(unique(spike_templates)), ...
    unique(spike_templates),20,'k','filled','ButtonDownFcn',@unit_click);
multiunit_lines = arrayfun(@(x) line(xlim,[0,0],'linewidth',2,'visible','off'),1:2);
xlim(unit_axes,[-0.1,1]);
ylim([-50, max(channel_positions(:,2))+50]);
ylabel('Depth (\mum)')
xlabel('Normalized log rate')

title('Click MUA borders');
user_click_coords = ginput;
depth_group_edges = user_click_coords(:,2);
yline(depth_group_edges,'linewidth',2,'color','r');
depth_group_centers = movmean(depth_group_edges,2,'endpoints','discard');
text(zeros(length(depth_group_centers),1),depth_group_centers, ...
    num2cell(1:length(depth_group_centers)),'FontSize',20','color','r');
drawnow;

% Make depth groups
n_depths = length(depth_group_edges) - 1;
depth_group = discretize(spike_depths,depth_group_edges);

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% PSTH for all conditions
if contains(bonsai_workflow,'lcr')
    % (L/C/R passive)
    stim_x = vertcat(trial_events.values.TrialStimX);
    align_times = cellfun(@(x) stimOn_times(stim_x == x),num2cell(unique(stim_x)),'uni',false);
elseif contains(bonsai_workflow,'stim_wheel')
    % (task)
    align_times = {stimOn_times,stim_move_time,reward_times};
end

depth_psth = nan(n_depths,length(t_bins)-1,2);
for curr_align = 1:length(align_times)
    use_align = align_times{curr_align};
    t_peri_event = bsxfun(@plus,use_align,t_bins);
    for curr_depth = 1:n_depths        
        curr_spikes = spike_times_timeline(depth_group == curr_depth);
        
        curr_spikes_binned = cell2mat(arrayfun(@(x) ...
            histcounts(curr_spikes,t_peri_event(x,:)), ...
            [1:size(t_peri_event,1)]','uni',false))./psth_bin_size;
        curr_mean_psth = mean(curr_spikes_binned,1);
        
        depth_psth(curr_depth,:,curr_align) = curr_mean_psth;
    end
end

smooth_size = 100;
depth_psth_smooth = smoothdata(depth_psth,2,'gaussian',smooth_size);

figure; h = tiledlayout(1,length(align_times)+1);

% Plot units and depths
nexttile; set(gca,'YDir','reverse');hold on;
if exist('probe_areas','var')
    probe_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
        probe_areas{1}.color_hex_triplet,'uni',false)),[1,3,2]);

    probe_areas_boundaries = probe_areas{1}.probe_depth;
    probe_areas_centers = mean(probe_areas_boundaries,2);

    probe_areas_image_depth = 0:1:max(probe_areas_boundaries,[],'all');
    probe_areas_image_idx = interp1(probe_areas_boundaries(:,1), ...
        1:height(probe_areas{1}),probe_areas_image_depth, ...
        'previous','extrap');
    probe_areas_image = probe_areas_rgb(probe_areas_image_idx,:,:);

    image([0,1],probe_areas_image_depth,probe_areas_image);
    yline(unique(probe_areas_boundaries(:)),'color','k','linewidth',1);
    set(gca,'YTick',probe_areas_centers,'YTickLabels',probe_areas{1}.acronym);
end
norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
scatter(norm_spike_n,template_depths(unique(spike_templates)),20,'k','filled');
yline(depth_group_edges,'linewidth',2,'color','r');
xlim([0,1]);

% Plot MUA
depth_psth_smooth_norm = reshape(normalize(reshape(depth_psth_smooth, ...
    size(depth_psth_smooth,1),[]),2,'range'),size(depth_psth_smooth));
for curr_align = 1:length(align_times)
    nexttile; set(gca,'YDir','reverse'); hold on;

    curr_stackplot = -200*depth_psth_smooth_norm(:,:,curr_align) + ...
        depth_group_centers;
    plot(t_centers,curr_stackplot','k','linewidth',2);

    xline(0,'r');
end

linkaxes(h.Children,'y');
ylim([0,3840]);


%% Sparse noise? 

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


%% Widefield/ephys regression maps (MUA)

% Set upsample value for regression
upsample_factor = 1;
sample_rate = (1/mean(diff(wf_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

mua_method = 'click'; % depth, click

switch mua_method

    case 'depth'
        % (to group multiunit by evenly spaced depths)
        n_depths = 8;
        depth_group_edges = round(linspace(0,4000,n_depths+1));
        [depth_group_n,~,depth_group] = histcounts(spike_depths,depth_group_edges);
        depth_groups_used = unique(depth_group);
        depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

    case 'click'
        % (for clickable manual depths)
        h = figure('units','normalized','position',[0.05,0.2,0.1,0.7]);
        unit_axes = axes('YDir','reverse'); hold on; xlabel('Norm spike rate');ylabel('Depth');

        if exist('probe_areas','var')
            probe_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
                probe_areas{1}.color_hex_triplet,'uni',false)),[1,3,2]);

            probe_areas_boundaries = probe_areas{1}.probe_depth;
            probe_areas_centers = mean(probe_areas_boundaries,2);

            probe_areas_image_depth = 0:1:max(probe_areas_boundaries,[],'all');
            probe_areas_image_idx = interp1(probe_areas_boundaries(:,1), ...
                1:height(probe_areas{1}),probe_areas_image_depth, ...
                'previous','extrap');
            probe_areas_image = probe_areas_rgb(probe_areas_image_idx,:,:);

            image(unit_axes,[0,1],probe_areas_image_depth,probe_areas_image);
            yline(unique(probe_areas_boundaries(:)),'color','k','linewidth',1);
            set(unit_axes,'YTick',probe_areas_centers,'YTickLabels',probe_areas{1}.acronym);
        end

        norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
        unit_dots = scatter3( ...
            norm_spike_n,template_depths(unique(spike_templates)), ...
            unique(spike_templates),20,'k','filled','ButtonDownFcn',@unit_click);
        multiunit_lines = arrayfun(@(x) line(xlim,[0,0],'linewidth',2,'visible','off'),1:2);
        xlim(unit_axes,[-0.1,1]);
        ylim([-50, max(channel_positions(:,2))+50]);
        ylabel('Depth (\mum)')
        xlabel('Normalized log rate')

        title('Click MUA borders');
        user_click_coords = ginput;
        close(h);
        depth_group_edges = user_click_coords(:,2);
        yline(depth_group_edges,'linewidth',2,'color','r');
        depth_group_centers = movmean(depth_group_edges,2,'endpoints','discard');
        text(zeros(length(depth_group_centers),1),depth_group_centers, ...
            num2cell(1:length(depth_group_centers)),'FontSize',20','color','r');
end

% Draw units and borders
figure('units','normalized','position',[0.05,0.2,0.1,0.7]); 
unit_axes = axes('YDir','reverse'); hold on; xlabel('Norm spike rate');ylabel('Depth');

if exist('probe_areas','var')
    probe_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
        probe_areas{1}.color_hex_triplet,'uni',false)),[1,3,2]);

    probe_areas_boundaries = probe_areas{1}.probe_depth;
    probe_areas_centers = mean(probe_areas_boundaries,2);

    probe_areas_image_depth = 0:1:max(probe_areas_boundaries,[],'all');
    probe_areas_image_idx = interp1(probe_areas_boundaries(:,1), ...
        1:height(probe_areas{1}),probe_areas_image_depth, ...
        'previous','extrap');
    probe_areas_image = probe_areas_rgb(probe_areas_image_idx,:,:);

    image(unit_axes,[0,1],probe_areas_image_depth,probe_areas_image);
    yline(unique(probe_areas_boundaries(:)),'color','k','linewidth',1);
    set(unit_axes,'YTick',probe_areas_centers,'YTickLabels',probe_areas{1}.acronym);
end

norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
unit_dots = scatter3( ...
    norm_spike_n,template_depths(unique(spike_templates)), ...
    unique(spike_templates),20,'k','filled','ButtonDownFcn',@unit_click);
multiunit_lines = arrayfun(@(x) line(xlim,[0,0],'linewidth',2,'visible','off'),1:2);
xlim(unit_axes,[-0.1,1]);
ylim([-50, max(channel_positions(:,2))+50]);
ylabel('Depth (\mum)')
xlabel('Normalized log rate')

yline(depth_group_edges,'linewidth',2,'color','r');
depth_group_centers = movmean(depth_group_edges,2,'endpoints','discard');
text(zeros(length(depth_group_centers),1),depth_group_centers, ...
    num2cell(1:length(depth_group_centers)),'FontSize',20','color','r');
drawnow;

% Bin spikes
n_depths = length(depth_group_edges) - 1;
depth_group = discretize(spike_depths,depth_group_edges);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths   
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);    
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins); 
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

use_svs = 1:200;
kernel_t = [-0.5,0.5];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 20;
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

% Convert kernel to pixel space
r_px = plab.wf.svd2px(wf_U(:,:,use_svs),k);

AP_imscroll(r_px,kernel_frames/wf_framerate);
clim([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
colormap(AP_colormap('BWR'));
axis image;

% Get center of mass for each pixel
% (get max r for each pixel, filter out big ones)
r_px_max = squeeze(max(r_px,[],3));
r_px_max(isnan(r_px_max)) = 0;
r_px_max_norm = bsxfun(@rdivide,r_px_max, ...
    permute(max(reshape(r_px_max,[],n_depths),[],1),[1,3,2]));
r_px_max_norm(isnan(r_px_max_norm)) = 0;
r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depths,[1,3,2])),3)./sum(r_px_max_norm,3);

% Plot map of cortical pixel by preferred depth of probe
r_px_com_col = ind2rgb(round(mat2gray(r_px_com,[1,n_depths])*255),jet(255));
figure;
a1 = axes('YDir','reverse');
imagesc(wf_avg); colormap(gray); clim([0,prctile(wf_avg(:),99.7)]);
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

%% Widefield/ephys regression maps (single units)

% Set templates
use_templates = unique(spike_templates);

% Set upsample value for regression
upsample_factor = 1;
sample_rate = (1/mean(diff(wf_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Bin spikes
binned_spikes = zeros(length(use_templates),length(time_bins)-1);
for curr_unit_idx = 1:length(use_templates)
    curr_spike_times = spike_times_timeline(spike_templates == use_templates(curr_unit_idx));    
    binned_spikes(curr_unit_idx,:) = histcounts(curr_spike_times,time_bins); 
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

use_svs = 1:200;
kernel_t = [-0.2,0.2];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 20;
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

% Convert kernel to pixel space
r_px = plab.wf.svd2px(wf_U(:,:,use_svs),k);

AP_imscroll(r_px,kernel_frames/wf_framerate);
clim([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
colormap(AP_colormap('BWR'));
axis image;

% Combine time, sort by explained variance, plot
% r_px_timepoint = sqrt(squeeze(sum(r_px.^2,3)));
r_px_timepoint = squeeze(max(r_px,[],3));

[~,sort_idx] = sort(template_depths);
% [~,sort_idx] = sort(explained_var.total,'descend');
plot_templates = explained_var.total > 0;
plot_templates_sort = sort_idx(plot_templates(sort_idx));

AP_imscroll(r_px_timepoint(:,:,plot_templates_sort),plot_templates_sort);
axis image;
clim([0,prctile(r_px_timepoint(:),99.9)]);




