function expscroll(U,V,wf_t,mousecam_fn,mousecam_t,trace,trace_t)
% expscroll(U,V,wf_t,mousecam_fn,mousecam_t,trace,trace_t)
%
% U = wf U
% V = wf V
% wf_t = wf frames in timeline time
%
% mousecam_fn = filename of mousecam movie
% mousecam_t = mousecam frames in timelite time
%
% trace = trace to plot (> 3 traces plots as image)
% trace_t = trace times

disp('Press s to save section of movies');

if ~exist('wf_t','var') || isempty(wf_t)
    wf_t = 1:size(V,2);
end

if exist('mousecam_fn','var') && exist('mousecam_t','var') && ~isempty(mousecam_t)
   handles.plot_mousecam = true; 
else
    handles.plot_mousecam = false;
end

if exist('trace','var') && exist('trace_t','var')
   handles.plot_trace = true; 
else
    handles.plot_trace = false;
end

% Set up plots
gui_fig = figure; colormap(gray)
set(gui_fig,'WindowScrollWheelFcn',{@imgSlider_MouseWheel, gui_fig});
set(gui_fig, 'KeyPressFcn', {@im_keypress, gui_fig});

tiledlayout(gui_fig,'flow');

handles.wf_axis = nexttile;
if handles.plot_mousecam
    handles.mousecam_axis = nexttile;
end
if handles.plot_trace
    handles.trace_axis = nexttile;
end

% Set up widefield data and image
handles.V = V; 
handles.U = U;
n_frames = size(V,2);
framerate = median(diff(wf_t));

handles.wf_t = wf_t;

wf_frame = 1;
handles.t = wf_t(wf_frame);

wf_im = plab.wf.svd2px(U,V(:,wf_frame));
handles.wf_im = imagesc(handles.wf_axis,wf_im); axis(handles.wf_axis,'off','image');

% Set up color axis (std non-linear so not accurate, but ballpark)
avg_im = plab.wf.svd2px(U,nanstd(V,[],2));
clim(handles.wf_axis,(prctile(abs(avg_im(:)),95)*5)*[-1,1]);

% Set up videoreaders, relative times, and images of cameras
discretize_times = [wf_t;wf_t(end)+framerate];
wf_frame_idx = wf_t;

if handles.plot_mousecam
    handles.mousecam_vr = VideoReader(mousecam_fn);
    handles.mousecam_frame_idx = discretize(mousecam_t,discretize_times);
        
    [~,mousecam_frame] = min(abs(wf_frame - handles.mousecam_frame_idx));
    mousecam_im = read(handles.mousecam_vr,mousecam_frame);
    handles.mousecam_im = imagesc(handles.mousecam_axis,mousecam_im); axis(handles.mousecam_axis,'off','image');
    caxis(handles.mousecam_axis,[0,255]);
end

% Set up trace
trace_time_surround = 1; % in seconds
if handles.plot_trace
    if size(trace,2) < 3
        handles.trace_plot = plot(handles.trace_axis,trace_t,trace,'linewidth',3);
        handles.trace_xlim = [handles.t-trace_time_surround,handles.t+trace_time_surround];
        handles.trace_tmark = line([handles.t,handles.t],ylim,'color','k');
    else
        handles.trace_plot = imagesc(handles.trace_axis,trace_t,1:size(trace,2),trace');
        handles.trace_xlim = [handles.t-trace_time_surround,handles.t+trace_time_surround];
        handles.trace_tmark = line([handles.t,handles.t],ylim,'color','r');
    end    
    xlim(handles.trace_axis,handles.trace_xlim);
end

% Set up scrollbar (use timer function to prevent lag)
ypos = [0 0 1 0.05];
handles.imgSlider = uicontrol('style','slider','units','normalized','position',ypos,'min',0,'max',1,'value',0);
set(handles.imgSlider,'Callback',{@imgSlider_Listener_timer,gui_fig});
handles.scroll_timer = timer('ExecutionMode','singleShot','TimerFcn',{@imgSlider_Listener,gui_fig});

set(handles.imgSlider,'Min',1);
set(handles.imgSlider,'Max',n_frames);
set(handles.imgSlider,'Value',1);
set(handles.imgSlider,'SliderStep',[10/n_frames, 100/n_frames]);

% Set up time title
handles.time_text = uicontrol('Style','text','String', ...
    ['Time: ' num2str(handles.t) ,'s'],'FontSize',14,'Units', ...
    'Normalized','Position',[0.3,0.93,0.4,0.07]);

% Update guidata
handles.wf_frame = wf_frame;
guidata(gui_fig, handles);

drawnow;


function imgSlider_Listener_timer(currentObject, eventdata, gui_fig)
% Get guidata
handles = guidata(gui_fig);

% Only run if not already running
if strcmp(handles.scroll_timer.Running,'off')
    start(handles.scroll_timer);
end


function imgSlider_Listener(currentObject, eventdata, gui_fig)
% Executes whenever the slider is pressed

% Get guidata
handles = guidata(gui_fig);

% Get frame number from slider, round appropriately
wf_frame = get(handles.imgSlider,'Value');
wf_frame = round(wf_frame);
set(handles.imgSlider,'Value',wf_frame);

% Update the images
update_im(handles,gui_fig,wf_frame);


function imgSlider_MouseWheel(currentObject, eventdata, gui_fig)
% Executes when mouse wheel is scrolled in figure

% Get guidata
handles = guidata(gui_fig);

% Get current frame
wf_frame = handles.wf_frame;

% Update current frame based on mouse wheel
mouse_wheel_count = eventdata.VerticalScrollCount;
wf_frame = wf_frame + mouse_wheel_count;
if wf_frame > length(handles.wf_t)
    wf_frame = length(handles.wf_t);
elseif wf_frame < 1
    wf_frame = 1;
end

% Set the slider
set(handles.imgSlider,'Value',wf_frame);

% Update the images
update_im(handles,gui_fig,wf_frame);


function im_keypress(currentObject, eventdata, gui_fig)
% Executes when a key is pressed

% Get guidata
handles = guidata(gui_fig);

switch eventdata.Key
    
    % Save section of images as movie
    case 's'
        
        % Get options
        disp('Preparing to make movie:');
        movie_t = input('Start/stop time (e.g. [0 5]): ');
        movie_framerate = input('Framerate: ');
        [save_file,save_path] = uiputfile('.avi','Choose save location');
        save_filename = [save_path save_file];
        
        movie_wf_frames = find(handles.wf_t > movie_t(1) & ...
            handles.wf_t < movie_t(2));
        n_movie_frames = length(movie_wf_frames);
        
        % Run through selected frames and save
        disp('Recording...')
        movie_frames(n_movie_frames) = struct('cdata',[],'colormap',[]);
        for curr_movie_frame_idx = 1:n_movie_frames
            curr_movie_frame = movie_wf_frames(curr_movie_frame_idx);
            % Update images
            update_im(handles,gui_fig,curr_movie_frame);
            movie_frames(curr_movie_frame_idx) = getframe(gui_fig);
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
guidata(gui_fig, handles);


function update_im(handles, gui_fig, wf_frame)

handles.t = handles.wf_t(wf_frame);

% Update the images (with the closest frame in case skip/offset)
if handles.plot_mousecam
    [~,mousecam_frame] = min(abs(wf_frame - handles.mousecam_frame_idx));
    mousecam_im = read(handles.mousecam_vr,mousecam_frame);
    set(handles.mousecam_im,'Cdata',mousecam_im);
end

wf_im = plab.wf.svd2px(handles.U,handles.V(:,wf_frame));

set(handles.wf_im,'CData',wf_im);

% Update trace
if handles.plot_trace
    xrange = diff(handles.trace_xlim);
    handles.trace_xlim = [handles.t-xrange/2,handles.t+xrange/2];
    xlim(handles.trace_axis,handles.trace_xlim);
    set(handles.trace_tmark,'XData',[handles.t,handles.t]);
end

% Update the time text
set(handles.time_text,'String',['Time: ' num2str(handles.t) ,'s']);

% Update guidata
handles.wf_frame = wf_frame;
guidata(gui_fig, handles);

drawnow;


