% Exp Control overhaul

function expcontrol_test

% Create figure
gui_fig = uifigure('Name','Experiment controller', ...
    'units','normalized','position',[0.01,0.1,0.2,0.8]);
gui_grid = uigridlayout(gui_fig,[4,1], ...
    'RowHeight',{'1x','1x','1x','1x','5x'},'BackgroundColor','w');
handles = struct;

% Connections lamps
connection_panel = uipanel(gui_grid,'Title','Connection status');
connection_modalities = {'Timelite','Bonsai','Mousecam','Widefield','Ephys'};
connection_ports = [ ...
    plab.locations.timelite_port, ...
    plab.locations.bonsai_port, ...
    plab.locations.mousecam_port, ...
    plab.locations.widefield_port, ...
    37497, ... % Open Ephys: constant in software
    ];

connection_tcpservers = arrayfun(@(x) tcpserver("0.0.0.0",x),connection_ports);

% (Bonsai server connection: set up callback for remote stopping)
configureCallback( ...
    connection_tcpservers([connection_tcpservers.ServerPort] == plab.locations.bonsai_port), ...
    "terminator", @(src,event,x) read_incoming(src,event,gui_fig));

connection_panel_grid = uigridlayout(connection_panel,[2,length(connection_modalities)]);
arrayfun(@(x) uilabel(connection_panel_grid,'Text',x,'HorizontalAlignment','Center'),connection_modalities);
handles.connection_lamps = arrayfun(@(x) uilamp(connection_panel_grid,'Color','y'),1:length(connection_modalities));

% Recording settings
rec_settings_panel = uipanel(gui_grid,'Title','Recording settings');
rec_settings_panel_grid = uigridlayout(rec_settings_panel,[2,2], ...
    'ColumnWidth',{'1x','4x'});

uilabel(rec_settings_panel_grid,'Text','Mouse','HorizontalAlignment','Center');
uibutton(rec_settings_panel_grid,'Text','Protocol','ButtonPushedFcn',{@choose_protocol,gui_fig});

handles.mouse = uieditfield(rec_settings_panel_grid,'ValueChangingFcn',{@update_controls,gui_fig});
handles.protocol = ...
    uilabel(rec_settings_panel_grid,'Text','<No protocol>','HorizontalAlignment','Center','BackgroundColor','w');


% Control buttons
exp_control_panel = uipanel(gui_grid,'Title','Control recording');
exp_control_panel_grid = uigridlayout(exp_control_panel,[1,2]);

handles.start = uibutton(exp_control_panel_grid,'Text','Start', ...
    'ButtonPushedFcn',{@recording_start,gui_fig}, ...
    'BackgroundColor',[0.6,1,0.6],'Enable',false);

handles.stop = uibutton(exp_control_panel_grid,'Text','Stop', ...
    'ButtonPushedFcn',{@recording_stop,gui_fig,true}, ...
    'BackgroundColor',[1,0.6,0.6],'Enable',false);


% Ephys section?
ephys_panel = uipanel(gui_grid,'Title','Ephys');


% Notes section
notes_panel = uipanel(gui_grid,'Title','Notes');
notes_panel_grid = uigridlayout(notes_panel,[4,1], ...
    'RowHeight',{'1x','10x','1x','10x'});

uilabel(notes_panel_grid,'Text','Recording notes');
handles.notes_recording = uitextarea(notes_panel_grid);

uilabel(notes_panel_grid,'Text','Day notes');
handles.notes_day = uitextarea(notes_panel_grid);

% Timer function to check connections
connections_timer = timer('Period',1,'ExecutionMode','fixedSpacing', ...
    'TimerFcn',{@connections_check,gui_fig});

% Store guidata
gui_data = struct;

gui_data.connection_tcpservers = connection_tcpservers;
gui_data.handles = handles;

guidata(gui_fig,gui_data);
start(connections_timer);

end


function connections_check(source,eventdata,gui_fig)
gui_data = guidata(gui_fig);
[gui_data.handles.connection_lamps([gui_data.connection_tcpservers.Connected]).Color] = deal('g');
[gui_data.handles.connection_lamps(~[gui_data.connection_tcpservers.Connected]).Color] = deal('r');
end


function choose_protocol(source,event,gui_fig)
bonsai_workflow_path = fullfile(extractBetween(string(which('plab.rig.bonsai_server')), ...
    '','PetersLab_rigging','Boundaries','inclusive'),'bonsai_workflows','*.bonsai');

[protocol_name, protocol_path] = uigetfile(bonsai_workflow_path);
gui_data = guidata(gui_fig);
gui_data.handles.protocol.Text = erase(protocol_name,'.bonsai');
gui_data.handles.protocol.UserData = fullfile(protocol_path,protocol_name);
update_controls([],[],gui_fig);
end


function update_controls(source,event,gui_fig)
gui_data = guidata(gui_fig);
mouse = gui_data.handles.mouse.Value;
bonsai_filename = gui_data.handles.protocol.UserData;
gui_data.handles.start.Enable = ~isempty(mouse) && exist(bonsai_filename,'file');
end


function recording_start(source,event,gui_fig)

gui_data = guidata(gui_fig);

% Grab recording information
mouse = string(gui_data.handles.mouse.Value);
bonsai_filename = string(gui_data.handles.protocol.UserData);
rec_day = string(datetime('now','Format','yyyy-MM-dd'));
rec_time = string(datetime('now','Format','HHmm'));

% Create recording path
recording_path = plab.locations.filename('server',mouse,rec_day,rec_time);
if exist(plab.locations.server_data_path,'dir')
    mkdir(recording_path);
else
    error('Server inaccessible');
end

% Disable start, enable stop
gui_data.handles.start.Enable = false;
gui_data.handles.stop.Enable = true;

% Broadcast start messages
recording_info = struct( ...
    'mouse', mouse, ...
    'date', rec_day, ...
    'time', rec_time, ...
    'bonsai_filename', bonsai_filename);
recording_info_json = jsonencode(recording_info);

connected_tcp = [gui_data.connection_tcpservers.Connected];
arrayfun(@(tcp) writeline(tcp,recording_info_json), ...
    gui_data.connection_tcpservers(connected_tcp));

end

function recording_stop(source,event,gui_fig,user_confirm)

gui_data = guidata(gui_fig);

% User confirm (if selected)
if user_confirm
    user_confirm_choice = uiconfirm(gui_fig,'Stop recording?','Confirm stop');
    if ~strcmpi(user_confirm_choice,'ok')
        return
    end
end

% Broadcast stop message
connected_tcp = [gui_data.connection_tcpservers.Connected];
arrayfun(@(tcp) writeline(tcp,'stop'), ...
    gui_data.connection_tcpservers(connected_tcp));

% Enable start, disable stop
gui_data.handles.start.Enable = true;
gui_data.handles.stop.Enable = false;

end


function read_incoming(source,event,gui_fig)
gui_data = guidata(gui_fig);

incoming_message = readline(source);

if strcmp(incoming_message,'Bonsai finished')
    recording_stop(source,event,gui_fig,false)
end
end

       