function daq_test

% Requires: 
% data acquisition toolbox
% instrument control toolbox

%% Make GUI

gui_fig = figure('Units','Normalized','Position',[0.01,0.2,0.1,0.2],'color','w');

% Control buttons
clear controls_h
controls_h(1) = uicontrol('Parent',gui_fig,'Style','pushbutton', ...
    'String','Run','Callback',{@daq_run,gui_fig});
controls_h(end+1) = uicontrol('Parent',gui_fig,'Style','togglebutton', ...
    'String','Listen','Callback',{@daq_listen,gui_fig});
controls_h(end+1) = uicontrol('Parent',gui_fig,'Style','pushbutton', ...
    'String','Stop','Callback',{@daq_stop,gui_fig});
set(controls_h,'Units','normalized','FontSize',12, ...
    'BackgroundColor','w','Position',[0,0,1,1/length(controls_h)]);

set(gcf,'Children',flipud(get(gcf,'Children')));
align(fliplr(controls_h),'center','fixed',0);

%% Set up DAQ

daq_sample_rate = 4000;
daq_samples_per_notify = 1/daq_sample_rate; % collect 1s of data per cycle

% Find DAQ
daqs_available = daqlist;
% (for now: assume first available DAQ is the one to use)
use_daq = 1;
daq_device = daq(daqs_available.VendorID(use_daq));

% Set DAQ properties 
daq_device.Rate = daq_sample_rate;
daq_device.ScansAvailableFcn = @(src,evt,x) daq_live_plot_test(src,evt,gui_fig);
daq_device.ScansAvailableFcnCount = daq_sample_rate;

ch = addinput(daq_device,daqs_available.DeviceID(use_daq),'ai0','Voltage');
ch.TerminalConfig = 'SingleEnded';
ch.Name = 'ai_test';

ch = addinput(daq_device,daqs_available.DeviceID(use_daq),'ctr0','Position');
ch.EncoderType = 'X4';
ch.Name = 'wheel';

% Initialize and upload gui data
gui_data = struct;
gui_data.daq_device = daq_device;

guidata(gui_fig,gui_data);

end

function daq_run(h,eventdata,gui_fig)

% Get gui data
gui_data = guidata(gui_fig);

% Start DAQ acquisition
% pause(0.5);
start(gui_data.daq_device,"continuous");

gui_data.live_plot = axes(figure);

% Update gui data
guidata(gui_fig,gui_data);

end

function daq_stop(h,eventdata,gui_fig)

% Get gui data
gui_data = guidata(gui_fig);

% Stop DAQ acquisition
stop(gui_data.daq_device)

close(gui_data.live_plot.Parent);

% Update gui data
guidata(gui_fig,gui_data);

end

function daq_live_plot_test(obj,event,gui_fig)

% Get gui data
gui_data = guidata(gui_fig);

% Read buffered data
data = read(obj,obj.ScansAvailableFcnCount);

plot_data_t = 10; % seconds of data to plot

if ~isfield(gui_data,'live_plot_traces') || ~any(isgraphics(gui_data.live_plot_traces))
    % If nothing plotted, create traces and plot data at end
    blank_data = nan(plot_data_t*obj.Rate,size(data,2));
    gui_data.live_plot_traces = plot(gui_data.live_plot,[0:size(blank_data,1)-1]/obj.Rate,blank_data);
end

% Shift off old data, swap in new data
old_plot_data = cell2mat(get(gui_data.live_plot_traces,'YData'))';
new_plot_data = circshift(old_plot_data,-size(data,1),1);
new_plot_data(end-size(data,1)+1:end,:) = data.Variables;

% Anything with a position measurement type (e.g. rotary encoder): 
% switch from unsigned to signed (remove overflow if negative)
position_channels_idx = strcmp({obj.Channels.MeasurementType},'Position');
counter_bits = 32;
signed_threshold = 2^(counter_bits-1);
new_plot_data(new_plot_data(:,position_channels_idx) > signed_threshold,position_channels_idx) = ...
    new_plot_data(new_plot_data(:,position_channels_idx) > signed_threshold,position_channels_idx) - 2^counter_bits;

% Draw new data to plot
new_plot_data_cell = mat2cell(new_plot_data,size(new_plot_data,1),ones(2,1));
[gui_data.live_plot_traces.YData] = deal(new_plot_data_cell{:});

% for wheel: storing here
%  plot(gui_data.live_plot,data.Time,typecast(uint32(data.wheel),'int32'))

% Update gui data
guidata(gui_fig,gui_data);

end



