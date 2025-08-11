function expscroll
% expscroll
%
% Scroll data (mousecam, widefield, ephys)

fprintf('s: save\np: autoplay\n');

% Pull variables from base workspace
gui_data = struct;
gui_data.plot_mousecam = false;
gui_data.plot_widefield = false;
gui_data.plot_ephys = false;
gui_data.plot_behavior = false;
gui_data.plot_screen = false;

% (timelite - required for timestamps)
try
    timelite = evalin('base','timelite');
catch me
    error('Timelite data not in workspace: %s',me);
end

% (mousecam)
try
    gui_data.mousecam_fn = evalin('base','mousecam_fn');
    gui_data.mousecam_t = evalin('base','mousecam_times');
    gui_data.plot_mousecam = true;
end

% (widefield)
try 
    gui_data.wf_U = evalin('base','wf_U');
    gui_data.wf_V = evalin('base','wf_V');
    gui_data.wf_t = evalin('base','wf_t');
    gui_data.plot_widefield = true;
end

% (ephys)
try 
    spike_times_timeline = evalin('base','spike_times_timelite');
    spike_templates = evalin('base','spike_templates');
    template_depths = evalin('base','template_depths');

    % Bin ephys into MUA by depth
    depth_corr_window = 50; % MUA window in microns
    depth_corr_window_spacing = 50; % MUA window spacing in microns

    max_depths = 3840;

    depth_corr_bins = [0:depth_corr_window_spacing:(max_depths-depth_corr_window); ...
        (0:depth_corr_window_spacing:(max_depths-depth_corr_window))+depth_corr_window];
    depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

    spike_binning_t = 0.005; % seconds
    spike_binning_t_edges = min(timelite.timestamps):spike_binning_t:max(timelite.timestamps);

    binned_spikes_depth = zeros(size(depth_corr_bins,2),length(spike_binning_t_edges)-1);
    for curr_depth = 1:size(depth_corr_bins,2)
        curr_depth_templates_idx = ...
            find(template_depths >= depth_corr_bins(1,curr_depth) & ...
            template_depths < depth_corr_bins(2,curr_depth));

        binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timeline( ...
            ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
    end

    % Smooth spikes
    smooth_size = 10;
    binned_spikes_depth_smooth = smoothdata(binned_spikes_depth,2,'gaussian',smooth_size);

    gui_data.ephys_t = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;
    gui_data.ephys = binned_spikes_depth_smooth;    

    gui_data.plot_ephys = true;
end

%(behavior)
try
    timelite = evalin('base','timelite');
    wheel_velocity = evalin('base','wheel_velocity');
    photodiode_trace = evalin('base','photodiode_trace');
    reward_thresh = evalin('base','reward_thresh');
    lick_thresh = evalin('base','lick_thresh');
    wheel_move= evalin('base','wheel_move');

    gui_data.behavior_t=timelite.timestamps;

    wheel_velocity_norm = (wheel_velocity - min(wheel_velocity)) / (max(wheel_velocity) - min(wheel_velocity));...

    gui_data.behavior_matrix=[(photodiode_trace>3)+4.4,...
        reward_thresh+3.3,...
        wheel_move+2.2,...
        wheel_velocity_norm+1.1,...
        lick_thresh];
    gui_data.plot_behavior = true;

end

%（screen）
try
    timelite = evalin('base','timelite');
    % stim_pos = evalin('base','stim_pos');
    gray_bg = evalin('base','gray_bg');
    % grating_pattern= evalin('base','grating_pattern');
    % screen_w= evalin('base','screen_w');
    % screen_h= evalin('base','screen_h');
    gui_data.screen_t=timelite.timestamps;
    gui_data.plot_screen = true;

end




% Set up time
gui_data.t = prctile(timelite.timestamps,[0,100]);
gui_data.t_curr = 0;

% Set up plots
gui_fig = figure('color','w'); colormap(gray)
set(gui_fig,'WindowScrollWheelFcn',{@imgSlider_MouseWheel, gui_fig});
set(gui_fig, 'KeyPressFcn', {@im_keypress, gui_fig});

tiledlayout(gui_fig,'flow');






% (mousecam)
if gui_data.plot_mousecam
    gui_data.mousecam_axis = nexttile;

    gui_data.mousecam_vr = VideoReader(gui_data.mousecam_fn);
    mousecam_im = read(gui_data.mousecam_vr,1);
    gui_data.mousecam_im = imagesc(gui_data.mousecam_axis,mousecam_im); 
    axis(gui_data.mousecam_axis,'off','image')
    colormap(gui_data.mousecam_axis,'gray');
    clim(gui_data.mousecam_axis,[0,255]);
end

% (behavior)
if gui_data.plot_behavior
    
    gui_data.behavior_axis = nexttile;
    plot(gui_data.behavior_axis,gui_data.behavior_t,gui_data.behavior_matrix)
    behavior_xlim = gui_data.t_curr + [-1,1];
     xlim(gui_data.behavior_axis,behavior_xlim);
     yticks([0 1.1 2.2 3.3 4.4])
     yticklabels({'lick','wheel velocity','wheel move','reward','stim'})
       gui_data.behavior_tmark = xline(gui_data.t_curr,'r');

ylim([-0.1 5.5])
end

% (widefield)
if gui_data.plot_widefield
    gui_data.wf_axis = nexttile;

    wf_im = plab.wf.svd2px(gui_data.wf_U,gui_data.wf_V(:,1));
    gui_data.wf_im = imagesc(gui_data.wf_axis,wf_im);
    axis(gui_data.wf_axis,'off','image');
    clim(gui_data.wf_axis,[-1,1].*median(std(gui_data.wf_V,[],2))*2);
    colormap(gui_data.wf_axis,ap.colormap('PWG'));
end

% (screen)
if gui_data.plot_screen
    gui_data.screen_axis = nexttile;
    frame=gray_bg;
    gui_data.screen_im = imagesc(gui_data.screen_axis,frame,[0 1]);
    axis(gui_data.screen_axis,'off','image');
    % colormap(gui_data.screen_axis,'gray');
end

% (ephys)
if gui_data.plot_ephys
    gui_data.ephys_axis = nexttile;

    ephys_time_surround = 1; % in seconds
    gui_data.ephys_plot = imagesc(gui_data.ephys_axis, ...
        gui_data.ephys_t,depth_corr_bin_centers,gui_data.ephys);
    ephys_xlim = gui_data.t_curr + [-1,1].*ephys_time_surround;
    gui_data.ephys_tmark = xline(gui_data.t_curr,'r');
    xlim(gui_data.ephys_axis,ephys_xlim);
    clim(gui_data.ephys_axis,[0,max(gui_data.ephys,[],'all')*0.5]);
    colormap(gui_data.ephys_axis,ap.colormap('WK'));

end

% Set up scrollbar (use timer function to prevent lag)
ypos = [0 0 1 0.05];
gui_data.imgSlider = uicontrol('style','slider','units','normalized','position',ypos,'min',0,'max',1,'value',0);
set(gui_data.imgSlider,'Callback',{@imgSlider_Listener_timer,gui_fig});
gui_data.scroll_timer = timer('ExecutionMode','singleShot','TimerFcn',{@imgSlider_Listener,gui_fig});

set(gui_data.imgSlider,'Min',0);
set(gui_data.imgSlider,'Max',1);
set(gui_data.imgSlider,'Value',0.01);
set(gui_data.imgSlider,'SliderStep',[0.01, 0.1]);

% Set up time title
gui_data.time_text = uicontrol('Style','text','String', ...
    ['Time: ' num2str(gui_data.t_curr) ,'s'],'FontSize',14,'Units', ...
    'Normalized','Position',[0.3,0.93,0.4,0.07],'BackgroundColor','w');

% Set up autoplay
gui_data.autoplay = timer('TimerFcn',{@autoplay,gui_fig}, ...
            'Period', 0.01, 'ExecutionMode','fixedSpacing', ...
            'TasksToExecute', inf);

% Update guidata
guidata(gui_fig, gui_data);

drawnow;


function imgSlider_Listener_timer(currentObject, eventdata, gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Only run if not already running
if strcmp(gui_data.scroll_timer.Running,'off')
    start(gui_data.scroll_timer);
end


function imgSlider_Listener(currentObject, eventdata, gui_fig)
% Executes whenever the slider is pressed

% Get guidata
gui_data = guidata(gui_fig);

% Get frame number from slider, round appropriately
t_frac = get(gui_data.imgSlider,'Value');
t_new = diff(gui_data.t)*t_frac + gui_data.t(1);

% Update the images
update_im(gui_data,gui_fig,t_new);


function imgSlider_MouseWheel(currentObject, eventdata, gui_fig)
% Executes when mouse wheel is scrolled in figure

% Get guidata
gui_data = guidata(gui_fig);

% Update current frame based on mouse wheel
t_step = 0.05;
mouse_wheel_count = eventdata.VerticalScrollCount;
t_new = gui_data.t_curr + mouse_wheel_count*t_step;

% Set limit within t bounds
if t_new < gui_data.t(1)
    t_new = gui_data.t(1);
end
if t_new > gui_data.t(2)
    t_new = gui_data.t(2);
end

% Set the slider
t_prct = t_new/(gui_data.t(1) + diff(gui_data.t));
set(gui_data.imgSlider,'Value',t_prct);

% Update the images
update_im(gui_data,gui_fig,t_new);


function im_keypress(currentObject, eventdata, gui_fig)
% Executes when a key is pressed

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    case 'p'
        % Toggle auto-play
        if strcmp(gui_data.autoplay.Running,'off')
            start(gui_data.autoplay);
        elseif strcmp(gui_data.autoplay.Running,'on')
            stop(gui_data.autoplay);
        end

    case 's'
        % Save section of images as movie

        % Get options
        disp('Preparing to make movie:');
        movie_t_limit = input('Start/stop time (e.g. [0 5]): ');
        movie_framerate = input('Framerate: ');
        [save_file,save_path] = uiputfile('.avi','Choose save location');
        save_filename = [save_path save_file];
        
        movie_t = movie_t_limit(1):1/movie_framerate:movie_t_limit(2);
        n_movie_frames = length(movie_t);
        
        % Run through selected frames and save
        disp('Recording...')
        movie_frames(n_movie_frames) = struct('cdata',[],'colormap',[]);
        for curr_t_idx = 1:n_movie_frames
            curr_t = movie_t(curr_t_idx);
            % Update images
            update_im(gui_data,gui_fig,curr_t);
            movie_frames(curr_t_idx) = getframe(gui_fig);
        end
        
        % Write movie
        disp('Saving...')
        writerObj = VideoWriter(save_filename);
        writerObj.FrameRate = movie_framerate;
        open(writerObj);
        writeVideo(writerObj,movie_frames);
        close(writerObj);
        
        disp('Done.')
        
end

% Update guidata
guidata(gui_fig, gui_data);


function update_im(gui_data, gui_fig, t_new)

% Get guidata
gui_data = guidata(gui_fig);

% Update the images (with closest to time)
if gui_data.plot_mousecam
    [~,mousecam_frame] = min(abs(t_new - gui_data.mousecam_t));
    mousecam_im = read(gui_data.mousecam_vr,mousecam_frame);
    set(gui_data.mousecam_im,'Cdata',mousecam_im);
end

if gui_data.plot_widefield
    [~,wf_frame] = min(abs(t_new - gui_data.wf_t));
    wf_im = plab.wf.svd2px(gui_data.wf_U,gui_data.wf_V(:,wf_frame));
    set(gui_data.wf_im,'CData',wf_im);
end

if gui_data.plot_ephys
    xrange = 2.*[-0.5,0.5] + t_new;
    xlim(gui_data.ephys_axis,xrange);
    set(gui_data.ephys_tmark,'Value',t_new);
end

if gui_data.plot_behavior
    xrange = [-1,1] + t_new;
    xlim(gui_data.behavior_axis,xrange);
   set(gui_data.behavior_tmark,'Value',t_new);
end

if gui_data.plot_screen
    gray_bg = evalin('base','gray_bg');
    stim_pos = evalin('base','stim_pos');
    screen_w = evalin('base','screen_w');
    screen_h = evalin('base','screen_h');
    radius = evalin('base','radius');
    grating_pattern= evalin('base','grating_pattern');

        [~,screen_frame] = min(abs(t_new - gui_data.screen_t));
         
        frame = gray_bg;  % 每次重新生成背景
    if ~isnan(stim_pos(screen_frame,1))
        cx = round(stim_pos(screen_frame,1));
        cy = round(stim_pos(screen_frame,2));
        x_range = (cx-radius):(cx+radius);
        y_range = (cy-radius):(cy+radius);
        
        valid_x = x_range > 0 & x_range <= screen_w;
        valid_y = y_range > 0 & y_range <= screen_h;
        
        frame(y_range(valid_y), x_range(valid_x)) = ...
            grating_pattern(valid_y, valid_x);
    end

   set(gui_data.screen_im,'CData',frame);
end

% Update the time text
set(gui_data.time_text,'String',sprintf('%.2f s',t_new));

% Update current time
gui_data.t_curr = t_new;

% Update guidata
guidata(gui_fig, gui_data);

drawnow;

function autoplay(obj,event,gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

t_step = 0.05;
t_new = gui_data.t_curr + t_step;

% Set the slider
t_prct = t_new/(gui_data.t(1) + diff(gui_data.t));
set(gui_data.imgSlider,'Value',t_prct);

% Update the images
update_im(gui_data,gui_fig,t_new);






