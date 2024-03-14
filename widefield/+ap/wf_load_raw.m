function im = wf_load_raw(animal,rec_day,rec_time,frames,scroll_color)
% im = wf_load_raw(animal,rec_day,rec_time,frames,scroll_color)
% Grab raw data captured from (captured with plab.rig.widefield)
%
% frames = number: Load raw frames from binary 
% frames = 'scroll': scroll through raw data
%   in this mode: scroll_color = [n colors, plot color]

data_path = plab.locations.filename('server',animal,rec_day,[],'widefield');
if ~exist(data_path,'dir')
    error('Data path not found: %s',data_path);
end

%% Get metadata

% Get filename
metadata_fn = fullfile(data_path,sprintf('widefield_%s_metadata.bin',rec_time));

% Load metadata
% Metadata format: 9 x n frames
% [image height, image width, frame number, timestamp (y,m,d,H,M,S)]
curr_metadata_fid = fopen(metadata_fn,'r');
curr_metadata = reshape(fread(curr_metadata_fid,'double'),9,[]);
fclose(curr_metadata_fid);

% Get image info
n_frames = size(curr_metadata,2);
im_size = reshape(curr_metadata(1:2,1),1,[]);
frame_timestamps = datetime( ...
    curr_metadata(4,:),curr_metadata(5,:),curr_metadata(6,:), ...
    curr_metadata(7,:),curr_metadata(8,:),curr_metadata(9,:));

%% Load images (if frames specified)

if isnumeric(frames)
    % Open file for reading
    data_fn = fullfile(data_path,sprintf('widefield_%s_data.bin',rec_time));
    data_fid = fopen(data_fn,'r');

    % Load all selected frames
    im = zeros([im_size,length(frames)],'uint16');
    for curr_frame_idx = 1:length(frames)
        curr_frame_location = prod(im_size)*(frames(curr_frame_idx)-1)*2; % uint16: *2 bytes
        fseek(data_fid,curr_frame_location,-1);
        im(:,:,curr_frame_idx) = reshape(fread(data_fid,prod(im_size),'uint16=>uint16'),im_size);
    end
end


%% Set up GUI (if frames='scroll');

if strcmp(frames,'scroll')

    % Create figure for scrolling and ROIs
    gui_fig = figure;
    set(gui_fig, ...
        'WindowScrollWheelFcn',{@im_change, gui_fig}, ...
        'closeRequestFcn',{@close_gui,gui_fig});
    gui_data = struct;

    % Open file for reading
    data_fn = fullfile(data_path,sprintf('widefield_%s_data.bin',rec_time));
    gui_data.data_fid = fopen(data_fn,'r');

    % Store image info
    if ~exist('scroll_color','var')
        scroll_color = [1,1];
    end
    gui_data.scroll_color = scroll_color;
    gui_data.n_frames = n_frames;
    gui_data.im_size = im_size;

    % Set up scrollbar
    ypos = [0 0 1 0.05];
    gui_data.imgSlider = uicontrol('style','slider','units','normalized','position',ypos,'min',0,'max',1,'value',0);
    set(gui_data.imgSlider,'Callback',{@im_change, gui_fig});

    set(gui_data.imgSlider,'Enable','on');
    set(gui_data.imgSlider,'Min',1);
    set(gui_data.imgSlider,'Max',n_frames);
    set(gui_data.imgSlider,'Value',1);
    set(gui_data.imgSlider,'SliderStep',[10/n_frames, 10/n_frames]);

    % Create axes, plot first image
    gui_data.curr_frame = scroll_color(2);
    gui_data.im = imagesc(0);
    clim([0,30000]);
    axis image off;
    colormap(gray);

    % Update guidata
    guidata(gui_fig,gui_data);

    % Draw first image
    update_im(gui_fig);

end

end

function update_im(gui_fig)

% Get guidata
gui_data = guidata(gui_fig);

% Load current frame
curr_frame_location = prod(gui_data.im_size)*(gui_data.curr_frame-1)*2; % uint16: *2 bytes
fseek(gui_data.data_fid,curr_frame_location,-1);
curr_im = reshape(fread(gui_data.data_fid,prod(gui_data.im_size), ...
    'uint16=>uint16'),gui_data.im_size);

% Update image
gui_data.im.CData = curr_im;

end


function im_change(obj, eventdata, gui_fig)
% Executes when mouse wheel is scrolled in figure

% Get guidata
gui_data = guidata(gui_fig);

% Update image based on control modality
if isprop(eventdata,'VerticalScrollCount')
    % Mouse wheel was moved
    mouse_wheel_count = eventdata.VerticalScrollCount;
    new_frame = gui_data.curr_frame + mouse_wheel_count*gui_data.scroll_color(1);
elseif obj == gui_data.imgSlider
    % If slider was moved
    curr_slider = round(get(gui_data.imgSlider,'Value'));
    new_frame = (curr_slider-mod(curr_slider,gui_data.scroll_color(1)))+gui_data.scroll_color(2);
end

if new_frame > gui_data.n_frames
    new_frame = gui_data.n_frames;
elseif new_frame < 1
    new_frame = 1;
end
gui_data.curr_frame = new_frame;

% Update slider
set(gui_data.imgSlider,'Value',new_frame);

% Update guidata
guidata(gui_fig, gui_data);

% Update image
update_im(gui_fig);

end


function close_gui(obj, eventdata, gui_fig)

% Get guidata
gui_data = guidata(gui_fig);

% Close file
fclose(gui_data.data_fid);

% Close figure
delete(gui_fig);

end



