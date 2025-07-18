%% Exploratory ephys analysis

%% Bin and plot recording by depth

% Bin ephys into MUA by depth
depth_corr_window = 50; % MUA window in microns
depth_corr_window_spacing = 50; % MUA window spacing in microns

max_depths = 3840;

depth_corr_bins = [0:depth_corr_window_spacing:(max_depths-depth_corr_window); ...
    (0:depth_corr_window_spacing:(max_depths-depth_corr_window))+depth_corr_window];
depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

spike_binning_t = 0.005; % seconds
spike_binning_t_edges = min(timelite.timestamps):spike_binning_t:max(timelite.timestamps);
spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;

binned_spikes_depth = zeros(size(depth_corr_bins,2),length(spike_binning_t_edges)-1);
for curr_depth = 1:size(depth_corr_bins,2)
    curr_depth_templates_idx = ...
        find(template_depths >= depth_corr_bins(1,curr_depth) & ...
        template_depths < depth_corr_bins(2,curr_depth));

    % (to plot only responsive cells)
    % curr_depth_templates_idx = ...
    %     intersect(find(template_depths >= depth_corr_bins(1,curr_depth) & ...
    %     template_depths < depth_corr_bins(2,curr_depth)),responsive_units);

    binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timelite( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

% Smooth spikes
smooth_size = 10;
binned_spikes_depth_smooth = smoothdata(binned_spikes_depth,2,'gaussian',smooth_size);

% Plot
figure;
h = tiledlayout(5,1,'TileSpacing','none');
nexttile;
plot(timelite.timestamps,wheel_velocity,'k');
nexttile([4,1]);
imagesc(spike_binning_t_centers,[],binned_spikes_depth_smooth)
axis tight;
colormap(AP_colormap('WK'));
linkaxes(get(h,'Children'),'x');
clim([0,max(binned_spikes_depth_smooth,[],'all')*0.5])

% % Draw stim times (if stim presented)
% if isfield(trial_events.values,'TrialStimX')
%     stim_azimuth = vertcat(trial_events.values.TrialStimX);
%     xline(stimOn_times(stim_azimuth == -90),'b');
%     xline(stimOn_times(stim_azimuth == 0),'k');
%     xline(stimOn_times(stim_azimuth == 90),'r');
% end

%% MUA correlelogram

% Get correlation of MUA in sliding windows
depth_corr_window = 50; % MUA window in microns
depth_corr_window_spacing = 20; % MUA window spacing in microns

max_depths = 3840; % (hardcode, sometimes kilosort drops channels)

depth_corr_bins = [0:depth_corr_window_spacing:(max_depths-depth_corr_window); ...
    (0:depth_corr_window_spacing:(max_depths-depth_corr_window))+depth_corr_window];
depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

spike_binning_t = 0.05; % seconds
spike_binning_t_edges = nanmin(spike_times_timelite):spike_binning_t:nanmax(spike_times_timelite);

binned_spikes_depth = zeros(size(depth_corr_bins,2),length(spike_binning_t_edges)-1);
for curr_depth = 1:size(depth_corr_bins,2)
    curr_depth_templates_idx = ...
        find(template_depths >= depth_corr_bins(1,curr_depth) & ...
        template_depths < depth_corr_bins(2,curr_depth));
    
    binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timelite( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

mua_corr = corrcoef(binned_spikes_depth');

figure;
imagesc(depth_corr_bin_centers,depth_corr_bin_centers,mua_corr);
axis image;
clim([-1,1].*0.5);
colormap(ap.colormap('BWR'))


%% Cell raster

if contains(bonsai_workflow,{'passive','Image'})
    % (L/C/R passive)

    align_times_all = stimOn_times; 
    if isfield(trial_events.values,'TrialStimX')
        align_category_all = vertcat(trial_events.values.TrialStimX);
    elseif isfield(trial_events.values,'StimFrequence')
        align_category_all = vertcat(trial_events.values.StimFrequence);
    elseif isfield(trial_events.values,'PictureID')
        align_category_all = vertcat(trial_events.values.PictureID);
    end
    % (get only quiescent trials)
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        1:length(stimOn_times));
    align_times = align_times_all(quiescent_trials);
    align_category = align_category_all(quiescent_trials);

    ap.cellraster(align_times,align_category);

elseif contains(bonsai_workflow,'stim_wheel')
    % (regular task)
    wheel_starts = timelite.timestamps(diff([0;wheel_move]) == 1);
    wheel_stops = timelite.timestamps(diff([0;wheel_move]) == -1);
    
    % (get wheel starts when no stim on screen: not sure this works yet)
    iti_move_idx = interp1(photodiode_times, ...
        photodiode_values,wheel_starts,'previous') == 0;
    [~,iti_move_sortidx] = sort(wheel_stops(iti_move_idx) - ...
        wheel_starts(iti_move_idx));

    [~,rxn_sort_idx] = sort(stim_to_move);

    % (if multiple task types - split by task type)
    if ~isfield(trial_events.values,'TaskType')
        trial_sort_idx = rxn_sort_idx;
        reward_sort_idx = 1:length(reward_times);
    else
        [~,type_rxn_sort_idx] = sortrows([vertcat(trial_events.values(1:n_trials).TaskType),stim_to_move]);
        trial_sort_idx = [vertcat(trial_events.values(1:n_trials).TaskType),type_rxn_sort_idx,rxn_sort_idx];

        reward_sort_idx = vertcat(trial_events.values([trial_events.values.Outcome]==1).TaskType);
    end

    ap.cellraster({stimOn_times(1:n_trials),stim_move_time,wheel_starts(iti_move_idx),reward_times}, ...
        {trial_sort_idx,trial_sort_idx,iti_move_sortidx,reward_sort_idx});
end

set(gcf,'Name',sprintf('%s, %s',animal,rec_day));


%% PSTH - units

% PSTH for all conditions
% (get quiescent trials)
stim_window = [0,0.5];
quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
    timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
    timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
    (1:length(stimOn_times))');

if contains(bonsai_workflow,'lcr')
    % (vis passive)
    stim_x = vertcat(trial_events.values.TrialStimX);
    align_times = cellfun(@(x) stimOn_times(stim_x(1:length(stimOn_times)) == x & quiescent_trials),num2cell(unique(stim_x)),'uni',false);
elseif contains(bonsai_workflow,'hml')
    % (aud passive)
    stim_x = vertcat(trial_events.values.StimFrequence);
    align_times = cellfun(@(x) stimOn_times(stim_x == x & quiescent_trials),num2cell(unique(stim_x)),'uni',false);
elseif contains(bonsai_workflow,'stim_wheel')
    % (task)
    align_times = {stimOn_times,stim_move_time,reward_times};
end

[unit_psth,~,unit_psth_t] = ...
    ap.psth(spike_times_timelite,align_times,spike_templates, ...
    'smoothing',100,'norm_window',[-0.5,0],'softnorm',1);

% Plot sorted
% (sort depth)
[~,sort_idx] = sort(template_depths);
% 
% % (sort max time)
% [~,max_idx] = max(max(abs(unit_psth),[],3),[],2);
% [~,sort_idx] = sort(max_idx);

ap.imscroll(unit_psth(sort_idx,:,:));
clim([-5,5]);
colormap(AP_colormap('BWR'));


%% |--> Get responsive units

% Set event to get response
% (get quiescent trials)
stim_window = [0,0.5];
quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
    timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
    timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
    (1:length(stimOn_times))');

if contains(bonsai_workflow,'lcr')
    % (vis passive)
    stim_type = vertcat(trial_events.values.TrialStimX);
    use_align = stimOn_times(stim_type(1:length(stimOn_times)) == 90 & quiescent_trials);
elseif contains(bonsai_workflow,'hml')
    % (aud passive)
    stim_type = vertcat(trial_events.values.StimFrequence);
    use_align = stimOn_times(stim_type == 8000 & quiescent_trials);
elseif contains(bonsai_workflow,'stim_wheel')
    % (task)
    use_align = stimOn_times;
end

response_t = [0,0.2];

window_center = mean(response_t).*[-1,1];
bin_size = diff(response_t);
[~,event_spikes] = ap.psth(spike_times_timelite,use_align,spike_templates, ...
    'window',window_center,'bin_size',bin_size);

event_response = squeeze(mean(diff(event_spikes,[],2),1));

n_shuff = 1000;
event_response_shuff = cell2mat(arrayfun(@(shuff) ...
    squeeze(mean(diff(ap.shake(event_spikes,2),[],2),1)), ...
    1:n_shuff,'uni',false));

event_response_rank = tiedrank(horzcat(event_response,event_response_shuff)')';
event_response_p = event_response_rank(:,1)./(n_shuff+1);

% Plot responsive units by depth
unit_dots = ap.plot_unit_depthrate(spike_templates,template_depths,probe_areas);
unit_dots.CData = +([1,0,0].*(event_response_p > 0.95)) + ([0,0,1].*(event_response_p < 0.05));

% Plot rasters of responsive units (from above - if done)
psth_use_t = unit_psth_t >= response_t(1) & unit_psth_t <= response_t(2);

responsive_units = find(event_response_p < 0.05 | event_response_p > 0.95);

% % (sort by max amplitude from avg across alignments)
% [~,sort_idx] = sort(nanmean(unit_psth(responsive_units,psth_use_t,:),[2,3]));
% 
% % (sort by max time in single alignment)
% sort_align = 1;
% [~,max_t] = max(unit_psth(responsive_units,:,sort_align),[],2);
% [~,sort_idx] = sort(max_t);
%
% (sort by depth)
[~,sort_idx] = sort(template_depths(responsive_units));

ap.imscroll(unit_psth(responsive_units(sort_idx),:,:));
colormap(AP_colormap('BWR'));
clim([-5,5]);


%% |--> Plot TAN PSTH

% Classify cell types (if not done already)
if ~exist('tans','var')
    AP_longstriatum_classify_striatal_units
end

[~,tan_sort_idx] = sort(template_depths(striatum_celltypes.tan));
tan_idx = find(tans);
ap.imscroll(unit_psth(tan_idx(tan_sort_idx),:,:))
colormap(AP_colormap('BWR'));
clim([-2,2]);


%% PSTH - MUA by depth

mua_method = 'even'; % depth, click

switch mua_method

    case 'even'
        % (to group multiunit by evenly spaced depths)
        n_depths = 20;
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
            unique(spike_templates),20,'k','filled');
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

% Make depth groups
n_depths = length(depth_group_edges) - 1;
depth_group = discretize(spike_depths,depth_group_edges);

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% PSTH for all conditions
if contains(bonsai_workflow,'passive')
    if contains(bonsai_workflow,'lcr')
        % (L/C/R visual passive)
        stim_type = vertcat(trial_events.values.TrialStimX);
    elseif contains(bonsai_workflow,'hml')
        % (H/M/L auditory passive)
        stim_type = vertcat(trial_events.values.StimFrequence);
    end
    stim_window = [0,0.5];
    quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
        timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
        timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
        1:length(stimOn_times))';

    align_times = cellfun(@(x) stimOn_times(quiescent_trials & stim_type == x), ...
        num2cell(unique(stim_type)),'uni',false);

elseif contains(bonsai_workflow,'stim_wheel')
    % (task)
    align_times = {stimOn_times,stim_move_time,reward_times};
end

depth_psth = nan(n_depths,length(t_bins)-1,2);
for curr_align = 1:length(align_times)
    use_align = align_times{curr_align};
    t_peri_event = bsxfun(@plus,use_align,t_bins);
    for curr_depth = 1:n_depths
        curr_spikes = spike_times_timelite(depth_group == curr_depth);

        curr_spikes_binned = cell2mat(arrayfun(@(x) ...
            histcounts(curr_spikes,t_peri_event(x,:)), ...
            [1:size(t_peri_event,1)]','uni',false))./psth_bin_size;
        curr_mean_psth = mean(curr_spikes_binned,1);

        depth_psth(curr_depth,:,curr_align) = curr_mean_psth;
    end
end

smooth_size = 100;
depth_psth_smooth = smoothdata(depth_psth,2,'gaussian',smooth_size);

figure('units','normalized','position',[0.02,0.15,0.2,0.7]);
h = tiledlayout(1,2);

% Plot units and depths
ax_h = nexttile; set(gca,'YDir','reverse');hold on;
ap.plot_unit_depthrate(spike_templates,template_depths,probe_areas,ax_h)

% Plot MUA
depth_psth_smooth_baseline = nanmean(depth_psth_smooth(:,t_centers<0,:),2);
depth_psth_smooth_norm = (depth_psth_smooth-depth_psth_smooth_baseline)./depth_psth_smooth_baseline;

nexttile; set(gca,'YDir','reverse'); hold on;
align_col = lines(length(align_times));
yscale = 0.9*max(diff(depth_group_edges))/max(depth_psth_smooth_norm,[],'all');

for curr_align = 1:length(align_times)

    curr_stackplot = -yscale*depth_psth_smooth_norm(:,:,curr_align) + ...
        reshape(depth_group_centers,[],1);
    plot(t_centers,curr_stackplot','color',align_col(curr_align,:),'linewidth',2);

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
for px_x = 1:nX
    for px_y = 1:nY
        align_times = ...
            stim_times(find( ...
            (noise_locations(px_y,px_x,1:end-1) == 128 & ...
            noise_locations(px_y,px_x,2:end) == 255) | ...
            (noise_locations(px_y,px_x,1:end-1) == 128 & ...
            noise_locations(px_y,px_x,2:end) == 0))+1);

        stim_times_grid{px_y,px_x} = align_times;
    end
end       


% Vectorize stim times by y,x position
[stim_x,stim_y] = meshgrid(1:nX,1:nY);
stim_positions = cellfun(@(x,y,times) repmat([y,x],length(times),1), ...
    num2cell(stim_x),num2cell(stim_y),stim_times_grid,'uni',false);

% Pick unit
% use_spikes = spike_times_timelite(spike_templates == 276);
% use_spikes = spike_times_timelite(spike_templates == 180);
use_spikes = spike_times_timelite(ismember(spike_templates, ...
    find(template_depths > 1500 & template_depths < 2500)));

% Get stim times vector (x,y)
stim_aligned_avg = cell(nY,nX);
raster_window = [-0.05,0.4];
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

mua_method = 'even'; % striatum, even, click, define

switch mua_method

    case 'striatum'
        AP_longstriatum_find_striatum_depth
        depth_length = 200;
        depth_group_edges = striatum_depth(1):200:striatum_depth(end);
        [~,~,depth_group] = histcounts(spike_depths,depth_group_edges);
        depth_groups_used = unique(depth_group);
        depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

    case 'even'
        % (to group multiunit by evenly spaced depths)
        n_depths = 10;
        depth_group_edges = round(linspace(0,4000,n_depths+1));
        [~,~,depth_group] = histcounts(spike_depths,depth_group_edges);
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
            unique(spike_templates),20,'k','filled');
        multiunit_lines = arrayfun(@(x) line(xlim,[0,0],'linewidth',2,'visible','off'),1:2);
        xlim(unit_axes,[-0.1,1]);
        ylim([-50, max(channel_positions(:,2))+50]);
        ylabel('Depth (\mum)')
        xlabel('Normalized log rate')

        title('Click MUA borders');
        user_click_coords = ginput;
        close(h);
        depth_group_edges = user_click_coords(:,2);

    case 'define'
        depth_group_edges = cellfun(@str2num,inputdlg({'MUA start','MUA end'}));

    otherwise
        warning('No valid depth selection chosen');
        return
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
    unique(spike_templates),20,'k','filled');
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
    curr_spike_times = spike_times_timelite(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

% use_svs = 1:200;
% kernel_t = [-1,1];
% kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
% lambda = 100;
% zs = [false,false];
% cvfold = 5;
% return_constant = false;
% use_constant = true;

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
    ap.regresskernel(fVdf_deconv_resample, ...
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
r_px_max = squeeze(r_px(:,:,kernel_frames==0,:));% squeeze(max(r_px,[],3));
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

% Plot max maps vertically stacked
figure;
imagesc(reshape(permute(r_px_max,[1,3,2]),[],size(r_px_max,2)));
clim(max(abs(clim)).*[-1,1]*0.5);
colormap(AP_colormap('BWR'));
axis image off


%% |--> Predicted PSTH by MUA

% % (passive)
% stim_window = [0,0.5];
% quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
%     timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
%     timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
%     1:length(stimOn_times))';
% 
% stim_x = vertcat(trial_events.values.TrialStimX);
% align_times = stimOn_times(quiescent_trials & stim_x == 90);

% (task)
align_times = stimOn_times;

% Set times for PSTH
surround_window = [-0.2,1];
surround_samplerate = sample_rate;

t = surround_window(1):1/surround_samplerate:surround_window(2);
peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

aligned_trace_measured = reshape(interp1(time_bin_centers,binned_spikes_std',peri_event_t,'linear'), ...
    length(align_times),length(t),[]);
aligned_trace_predicted = reshape(interp1(time_bin_centers,predicted_spikes',peri_event_t,'linear'), ...
    length(align_times),length(t),[]);

aligned_trace_measured_mean = permute(nanmean(aligned_trace_measured,1),[3,2,1]);
aligned_trace_predicted_mean = permute(nanmean(aligned_trace_predicted,1),[3,2,1]);

% Normalize to baseline
softnorm = 0;
aligned_trace_measured_mean_baseline = nanmean(aligned_trace_measured_mean(:,t < 0),2);
aligned_trace_measured_mean_norm = (aligned_trace_measured_mean-aligned_trace_measured_mean_baseline)./(aligned_trace_measured_mean_baseline+softnorm);

aligned_trace_predicted_mean_baseline = nanmean(aligned_trace_predicted_mean(:,t < 0),2);
aligned_trace_predicted_mean_norm = (aligned_trace_predicted_mean-aligned_trace_predicted_mean_baseline)./(aligned_trace_predicted_mean_baseline+softnorm);

% Plot
yscale = max(aligned_trace_measured_mean_norm,[],'all');

figure; hold on;
AP_stackplot(aligned_trace_measured_mean_norm',t,yscale,[],'k');
AP_stackplot(aligned_trace_predicted_mean_norm',t,yscale,[],'r');
xline(0);

%% Widefield/ephys regression (JUST ITI) DOESN'T WORK
% THIS DOESNT WORK 

% Set upsample value for regression
upsample_factor = 1;
sample_rate = (1/mean(diff(wf_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = wf_t(find(wf_t > skip_seconds,1)):1/sample_rate:wf_t(find(wf_t-wf_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% (group multiunit by evenly spaced depths)
n_depths = 8;
depth_group_edges = round(linspace(0,4000,n_depths+1));
[depth_group_n,~,depth_group] = histcounts(spike_depths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

% Bin spikes
n_depths = length(depth_group_edges) - 1;
depth_group = discretize(spike_depths,depth_group_edges);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timelite(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

% Set times to use (no movement, no stim)
stim_on_epochs = interp1(photodiode_times,photodiode_values,time_bin_centers','previous');
wheel_move_epochs = interp1(timelite.timestamps,+wheel_move,time_bin_centers','previous');
use_t = ~stim_on_epochs & ~wheel_move_epochs;

use_svs = 1:100;
kernel_t = [-0.2,0.2];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 50;
zs = [false,false];
cvfold = 5;
return_constant = false;
use_constant = true;
discontinuities = [0;diff(use_t) == 1];

fVdf_deconv_resample = interp1(wf_t,wf_V(use_svs,:)',time_bin_centers)';

% Regress cortex to spikes
[k,~,explained_var] = ...
    ap.regresskernel(fVdf_deconv_resample(:,use_t), ...
    binned_spikes_std(:,use_t),kernel_frames, ...
    lambda,zs,cvfold,return_constant,use_constant,discontinuities(use_t));

% Apply kernel (from partial data) to whole dataset
% (DOESN'T WORK - maybe easier to modify function to have things included
% in regression and things not)
predicted_spikes = cell2mat(arrayfun(@(str) sum(cell2mat(arrayfun(@(v) ...
    conv(fVdf_deconv_resample(v,:),k(v,:,str),'same'), ...
    1:size(k,1),'uni',false)'),1),1:size(k,3),'uni',false)');

% Convert kernel to pixel space
r_px = plab.wf.svd2px(wf_U(:,:,use_svs),k);

AP_imscroll(r_px,kernel_frames/wf_framerate);
clim([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
colormap(AP_colormap('BWR'));
axis image;






%% Widefield/ephys regression maps (all single units)

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
    curr_spike_times = spike_times_timelite(spike_templates == use_templates(curr_unit_idx));
    binned_spikes(curr_unit_idx,:) = histcounts(curr_spike_times,time_bins);
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

use_svs = 1:100;
kernel_t = [-0.1,0.1];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 10;
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
% clim([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)]) % (takes too long)
clim([-0.02,0.02]);
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

%% |--> Predicted PSTH by unit

% LCR passive: align to quiescent stim onset
stim_window = [0,0.5];
quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
    timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
    timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
    1:length(stimOn_times))';

stim_x = vertcat(trial_events.values.TrialStimX);
align_times = stimOn_times(quiescent_trials & stim_x == 90);

% Set times for PSTH
surround_window = [-0.2,1];
surround_samplerate = sample_rate;

t = surround_window(1):1/surround_samplerate:surround_window(2);
peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

aligned_trace_measured = reshape(interp1(time_bin_centers,binned_spikes_std',peri_event_t,'linear'), ...
    length(align_times),length(t),[]);
aligned_trace_predicted = reshape(interp1(time_bin_centers,predicted_spikes',peri_event_t,'linear'), ...
    length(align_times),length(t),[]);

aligned_trace_measured_mean = permute(nanmean(aligned_trace_measured,1),[3,2,1]);
aligned_trace_predicted_mean = permute(nanmean(aligned_trace_predicted,1),[3,2,1]);

% Normalize to baseline
softnorm = 0.3;
aligned_trace_measured_mean_baseline = nanmean(aligned_trace_measured_mean(:,t < 0),2);
aligned_trace_measured_mean_norm = (aligned_trace_measured_mean-aligned_trace_measured_mean_baseline)./(aligned_trace_measured_mean_baseline+softnorm);

aligned_trace_predicted_mean_baseline = nanmean(aligned_trace_predicted_mean(:,t < 0),2);
aligned_trace_predicted_mean_norm = (aligned_trace_predicted_mean-aligned_trace_predicted_mean_baseline)./(aligned_trace_predicted_mean_baseline+softnorm);

% Plot depth-sorted
[~,sort_idx] = sort(template_depths);
figure;
subplot(1,2,1);
imagesc(t,[],aligned_trace_measured_mean_norm(sort_idx,:,:));
xline(0);
clim([-2,2]);
title('Measured');

subplot(1,2,2);
imagesc(t,[],aligned_trace_predicted_mean_norm(sort_idx,:,:));
xline(0);
clim([-2,2]);
title('Predicted');

colormap(AP_colormap('BWR'));




%% ~~~~~~~~ BATCH

%% Grab and plot histology pictures (pre-SMZ)

animals = {'AP022'};

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    histology_path = plab.locations.filename('server',animal,[],[],'histology','processed');

    if exist(histology_path,'dir')

        histology_dir = dir(fullfile(histology_path,'**','slice_*.tif'));
        [~,sort_idx] = natsortfiles({histology_dir.name});
        
        histology_im = cell(length(histology_dir),1);
        for curr_slice = 1:length(histology_dir)
            histology_im{curr_slice} = imread(fullfile( ...
                histology_dir(sort_idx(curr_slice)).folder, ...
                histology_dir(sort_idx(curr_slice)).name));
        end

        figure('Name',animal);
        montage(histology_im);
        title(animal);
        drawnow;
    else
        fprintf('No histology found: %s\n',animal);
    end
end

%% Grab and plot histology pictures (SMZ)

animal = 'AP025';

% Just load all images
histology_path = plab.locations.filename('server',animal,[],[],'histology');
histology_dir = dir(fullfile(histology_path,'*.tif'));

histology_filenames = cellfun(@(path,name) fullfile(path,name), ...
    {histology_dir.folder},{histology_dir.name},'uni',false);
[~,sort_idx] = natsortfiles(histology_filenames);

histology_im = cell(size(histology_dir));
for curr_im = 1:length(sort_idx)
    histology_im{curr_im} = tiffreadVolume( ...
        histology_filenames{sort_idx(curr_im)});
end

n_chan = size(histology_im{1},3);

% % Plot channel montage separately
% figure('Name',animal); tiledlayout(1,n_chan);
% for curr_chan = 1:n_chan
%     nexttile;
%     m = montage(cellfun(@(im) im(:,:,curr_chan),histology_im,'uni',false));
%     clim(mean(m.CData,[1,2])*[0.25,4]);
% end

% Grab image montage and display as RGB
im_montage = uint16([]);
chan_cols = [0,1,0;1,0,0];
for curr_chan = 1:n_chan
    h = figure;
    m = montage(cellfun(@(im) im(:,:,curr_chan),histology_im,'uni',false));
    im_montage = cat(3,im_montage,m.CData);
    close(h);
end

m_clim = arrayfun(@(chan) double(prctile(im_montage(:,:,chan),[20,90],'all')).*[1;1.5],1:n_chan,'uni',false);

im_montage_rgb = min(sum(cell2mat(arrayfun(@(chan) ...
    mat2gray(im_montage(:,:,chan),double(m_clim{chan})).* ...
    permute(chan_cols(chan,:),[1,3,2]), ...
    permute(1:n_chan,[1,3,4,2]),'uni',false)),4),1);

im_montage_rgb = min(sum(cell2mat(arrayfun(@(chan) ...
    mat2gray(im_montage(:,:,chan),double(m_clim{chan})).* ...
    permute(chan_cols(chan,:),[1,3,2]), ...
    permute(1:n_chan,[1,3,4,2]),'uni',false)),4),1);

figure;image(im_montage_rgb);axis image off;
title(animal);

%% Histology scroller (NEW)

animal = 'AM010';


% (pre-SMZ)
% histology_path = plab.locations.filename('server',animal,[],[],'histology','raw_combined');
% (SMZ)
histology_path = plab.locations.filename('server',animal,[],[],'histology','raw');


% Histology scroller
ap_histology.histology_scroll(histology_path)



%% Plot units by area

animal = 'AM027';
recordings = plab.find_recordings(animal);
recordings = recordings([recordings.ephys]);

figure('Name',animal);
h = tiledlayout(1,length(recordings));
for curr_recording = 1:length(recordings)

    % Grab pre-load vars
    preload_vars = who;

    % Load data
    rec_day = recordings(curr_recording).day;
    rec_time = recordings(curr_recording).recording{end};
   
    load_parts = struct;
    load_parts.ephys = true;
    ap.load_recording;

    unit_axes = nexttile; hold on;
    unit_axes.YDir = 'reverse';

    % Plot units (depth vs normalized rate) with areas
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

    unit_dots = scatter( ...
        norm_spike_n,template_depths(unique(spike_templates)),20,'k','filled');
    multiunit_lines = arrayfun(@(x) line(xlim,[0,0],'linewidth',2,'visible','off'),1:2);
    xlim(unit_axes,[-0.1,1]);
    ylim([-50, max(channel_positions(:,2))+50]);
    ylabel('Depth (\mum)')
    xlabel('Normalized log rate')
    title(rec_day);

    drawnow;

end

title(h,animal);

%% Batch MUA by depth

animal = 'AP024';

use_workflow = 'lcr_passive';
% use_workflow = 'hml_passive_audio';
% use_workflow = 'stim_wheel*';
recordings = plab.find_recordings(animal,[],use_workflow);
recordings = recordings([recordings.ephys]);

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t_centers = conv2(t_bins,[1,1]/2,'valid');

% Set depth groups
n_depths = 8;
depth_group_edges = round(linspace(0,4000,n_depths+1));

day_mua = cell(length(recordings),1);

for curr_recording = 1:length(recordings)

    % Grab pre-load vars
    preload_vars = who;

    % Load data
    rec_day = recordings(curr_recording).day;
    rec_time = recordings(curr_recording).recording{end};
   
    load_parts = struct;
    load_parts.ephys = true;
    ap.load_recording;

    % Make depth groups
    [depth_group_n,~,depth_group] = histcounts(spike_depths,depth_group_edges);
    depth_groups_used = unique(depth_group);
    depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

    n_depths = length(depth_group_edges) - 1;
    depth_group = discretize(spike_depths,depth_group_edges);

    % PSTH for all conditions
    if contains(bonsai_workflow,'passive')
        if contains(bonsai_workflow,'lcr')
            % (L/C/R visual passive)
            stim_type = vertcat(trial_events.values.TrialStimX);
        elseif contains(bonsai_workflow,'hml')
            % (H/M/L auditory passive)
            stim_type = vertcat(trial_events.values.StimFrequence);
        end
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            1:length(stimOn_times))';

        align_times = cellfun(@(x) stimOn_times(quiescent_trials & ...
            stim_type(1:length(stimOn_times)) == x), ...
            num2cell(unique(stim_type)),'uni',false);

    elseif contains(bonsai_workflow,'stim_wheel')
        % (task)
        align_times = {stimOn_times,stim_move_time,reward_times};
    end

    depth_psth = nan(n_depths,length(t_bins)-1,length(align_times));
    for curr_align = 1:length(align_times)

        use_align = align_times{curr_align};
        if isempty(use_align)
            continue
        end

        t_peri_event = bsxfun(@plus,use_align,t_bins);

        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timelite(depth_group == curr_depth);

            curr_spikes_binned = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./psth_bin_size;
            curr_mean_psth = mean(curr_spikes_binned,1);

            depth_psth(curr_depth,:,curr_align) = curr_mean_psth;
        end
    end

    smooth_size = 100;
    depth_psth_smooth = smoothdata(depth_psth,2,'gaussian',smooth_size);

    depth_psth_smooth_baseline = nanmean(nanmean(depth_psth_smooth(:,t_centers<0,:),2),3);

    softnorm = 2;
    depth_psth_smooth_norm = (depth_psth_smooth - depth_psth_smooth_baseline)./ ...
        (depth_psth_smooth_baseline+softnorm);

    % Store MUA
    day_mua{curr_recording} = depth_psth_smooth_norm;

    % Prep for next loop
    AP_print_progress_fraction(curr_recording,length(recordings));
    clearvars('-except',preload_vars{:});

end


% Plot MUA
figure; 
h = tiledlayout(1,length(day_mua));
for curr_day = 1:length(day_mua)

    nexttile(h); hold on;
    align_col = lines(size(day_mua{curr_day},3));
    for curr_align = 1:size(day_mua{curr_day},3)
        AP_stackplot(day_mua{curr_day}(:,:,curr_align)',t_centers, ...
            2,1,align_col(curr_align,:));
    end
    xline(0,'k');
    title(recordings(curr_day).day);

    drawnow;
end
linkaxes(h.Children);
title(h,sprintf('%s: %s',animal,use_workflow),'interpreter','none');



