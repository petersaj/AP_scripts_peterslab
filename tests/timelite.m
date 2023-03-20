function timelite
% Light-weight NI-DAQ acquisition GUI (lite version of Timeline)
%
% Requires: 
% data acquisition toolbox
% instrument control toolbox

%%%% TO DO: 
% - display daq/channel info?
% - set filename dynamically (append to current experiments)
% - make listen mode to set filename and run
% - upload to server

%% Make GUI

gui_fig = figure('Name','Timelite','Units','Normalized', ...
    'Position',[0,0.1,0.1,0.8],'color','w','menu','none');

% Control buttons
clear controls_h
% (text)
controls_h(1) = uicontrol('Parent',gui_fig,'Style','text', ...
    'String','Setting up DAQ...','FontSize',12, ...
    'HorizontalAlignment','left', ...
    'Units','normalized','BackgroundColor','w','Position',[0,0.5,1,0.5]);

% (buttons)
controls_h(end+1) = uicontrol('Parent',gui_fig,'Style','togglebutton', ...
    'String','Manual','Callback',{@daq_manual,gui_fig});
controls_h(end+1) = uicontrol('Parent',gui_fig,'Style','togglebutton', ...
    'String','Listen','Callback',{@daq_listen,gui_fig});
set(controls_h(2:end),'Units','normalized','FontSize',16, ...
    'BackgroundColor','w','Position',[0,0,1,0.5/(length(controls_h)-1)]);

set(gcf,'Children',flipud(get(gcf,'Children')));
align(fliplr(controls_h),'center','fixed',0);

% Disable controls until DAQ is set up
set(controls_h,'Enable','off');

% Draw initial controls (so user knows startup is in progress)
drawnow;

%% Set up DAQ

try
    % Set up DAQ according to local config file
    daq_device = plab.local_rig.timelite_config;
catch me
    % Error if DAQ could not be configured
    error('Timelite - DAQ device could not be configured: \n %s',me.message)
end

% Display DAQ channel info
channel_text = [ ...
    {sprintf('Sample rate:  %d',daq_device.Rate)}
    {sprintf('Sample rate:  %d',daq_device.Rate)}; ...
    {sprintf('Samples/upload:  %d \n',daq_device.ScansAvailableFcnCount)}
    {'Channels: '}; ...
    join([{daq_device.Channels.ID}',{daq_device.Channels.Name}',],': ')];

set(controls_h(1),'String',channel_text);

% Set DAQ update function
daq_device.ScansAvailableFcn = @(src,evt,x) daq_upload(src,evt,gui_fig);

% Initialize and upload gui data
gui_data = struct;
gui_data.controls_h = controls_h;
gui_data.daq_device = daq_device;

guidata(gui_fig,gui_data);

% When setup is complete, enable controls
set(controls_h,'Enable','on');

end

function daq_manual(h,eventdata,gui_fig)

switch h.Value
    case 1
        % Manual recording is turned on

        % Get gui data
        gui_data = guidata(gui_fig);

        % Change button display and disable other buttons
        h.String = 'Stop';
        h.BackgroundColor = [0.8,0,0];
        h.ForegroundColor = 'w';
        set(gui_data.controls_h(gui_data.controls_h ~= h),'Enable','off');

        % User choose mouse name
        mouse_name = cell2mat(inputdlg('Mouse name'));
        if isempty(mouse_name)
            % (if no mouse entered, do nothing)
            return
        end

        % Set save filename
        save_dir = "C:\Users\petersa\Desktop";
        gui_data.save_filename = fullfile(save_dir,[mouse_name,'_timelite.mat']);

        % Update gui data
        guidata(gui_fig,gui_data);

        % Start DAQ acquisition
        daq_start(gui_fig);

    case 0
        % Manual recording is turned off

        % Stop recording
        daq_stop(gui_fig)

        % Get gui data
        gui_data = guidata(gui_fig);

        % Change button display and disable other buttons
        h.String = 'Manual';
        h.BackgroundColor = 'w';
        h.ForegroundColor = 'k';
        set(gui_data.controls_h,'Enable','on');

end

end

function daq_listen(h,eventdata,gui_fig)
% Listen for start from experiment controller

tcpclient("163.1.249.17",plab.locations.timelite_port);
% TO DO: pull from bonsai_server to do the same thing

end

function daq_start(gui_fig)

% Get gui data
gui_data = guidata(gui_fig);

% Create save file
% (daq information)
daq_info = struct( ...
    'rate',gui_data.daq_device.Rate, ...
    'device',gui_data.daq_device.Channels(1).Device.Model, ...
    'type',{gui_data.daq_device.Channels.Type}, ...
    'channel',{gui_data.daq_device.Channels.ID}, ...
    'measurement_type',{gui_data.daq_device.Channels.MeasurementType}, ...
    'channel_name',{gui_data.daq_device.Channels.Name});
% (daq data - to be filled during streaming)
[timestamps,data] = deal([]);
% (save initial variables and keep open for streaming)
save(gui_data.save_filename,'daq_info','timestamps','data','-v7.3');
gui_data.save_file_mat = matfile(gui_data.save_filename,'Writable',true);

% Start DAQ acquisition
start(gui_data.daq_device,'continuous');

% Create live plot (unclosable)
live_plot_fig = figure('CloseRequestFcn',[],'color','w', ...
    'Units','normalized','Position',[0.5,0,0.5,1], ...
    'Name','Timelite live plot','menu','none');
gui_data.live_plot_fig = live_plot_fig;

% Update gui data
guidata(gui_fig,gui_data);

end

function daq_stop(gui_fig)

% Get gui data
gui_data = guidata(gui_fig);

% Stop DAQ acquisition
stop(gui_data.daq_device)

% Delete live plot
delete(gui_data.live_plot_fig);

% Update gui data
guidata(gui_fig,gui_data);

end

function daq_upload(obj,event,gui_fig)

% Get gui data
gui_data = guidata(gui_fig);

% Read buffered data
[daq_data,daq_timestamps] = ...
    read(obj,obj.ScansAvailableFcnCount,'OutputFormat','Matrix');

% Counter data: convert from unsigned to signed integer type
% (allow for negative values, rather than overflow)
% (assumes 32-bit counter)
position_channels_idx = strcmp({obj.Channels.MeasurementType},'Position');
daq_data(:,position_channels_idx) = ...
    double(typecast(uint32(daq_data(:,position_channels_idx)),'int32'));

% Plot data
daq_plot(obj,gui_data,daq_data,gui_fig);

% Save data (stream: append to mat file)
% (get index for new data)
curr_data_size = size(gui_data.save_file_mat.timestamps,1);
new_data_size = length(daq_timestamps);
new_data_row_idx = curr_data_size+1:curr_data_size+new_data_size;
% (write data appended to old data)
gui_data.save_file_mat.timestamps(new_data_row_idx,1) = daq_timestamps;
gui_data.save_file_mat.data(new_data_row_idx,1:size(daq_data,2)) = daq_data;

end

function daq_plot(obj,gui_data,daq_data,gui_fig)

plot_data_t = 10; % seconds of data to plot

if ~isfield(gui_data,'live_plot_traces') || ~any(isgraphics(gui_data.live_plot_traces))
    % If nothing plotted, create traces and plot data at end
    blank_data = zeros(plot_data_t*obj.Rate,size(daq_data,2));
    gui_data.live_plot_traces = stackedplot(gui_data.live_plot_fig,(0:size(blank_data,1)-1)/obj.Rate,blank_data);
    gui_data.live_plot_traces.DisplayLabels = {gui_data.daq_device.Channels.Name};
    % Update gui data
    guidata(gui_fig,gui_data);
end

% Shift off old data, swap in new data
old_plot_data = get(gui_data.live_plot_traces,'YData');
new_plot_data = circshift(old_plot_data,-size(daq_data,1),1);
new_plot_data(end-size(daq_data,1)+1:end,:) = daq_data;

% Draw new data to plot
gui_data.live_plot_traces.YData = new_plot_data;
gui_data.live_plot_traces.DisplayLabels = {gui_data.daq_device.Channels.Name};

% Update gui data
guidata(gui_fig,gui_data);

end








