classdef bonsai_server_test

    properties (SetAccess = protected)
        communication_handles
    end

    methods

        %% User functions
        function obj = bonsai_server_test
            % Initialize bonsai server

            % Set up structure for listeners and receivers
            obj.communication_handles = struct;

            % Set up experiment control communication
            obj.communication_handles.client_expcontrol = ...
                tcpclient(plab.local_rig.config.local.client,plab.locations.bonsai_port);

            % Set up Bonsai emitter
            obj.communication_handles.bonsai_udp = udpport("IPV4");

            % Set up Bonsai reciever
            % (Bonsai only sends in Open Sound Control (OSC) protocol)
            % (load the jar file)
            javaaddpath(fullfile(fileparts(which( ...
                'plab.rig.bonsai_server_helpers.bonsai_oscsend')), ...
                'javaosctomatlab.jar'));
            % (import java packages)
            import com.illposed.osc.*;
            import java.lang.String;

            % (set up OSC listener and receiver)
            oscport = 20000; % (defined in Bonsai workflow)
            obj.communication_handles.bonsai_receiver = OSCPortIn(oscport);
            obj.communication_handles.bonsai_listener = MatlabOSCListener();
            obj.communication_handles.bonsai_receiver.addListener(String('/stop'), ...
                obj.communication_handles.bonsai_listener);
            obj.communication_handles.bonsai_receiver.startListening();

            % Set up callback for experiment controller
            configureCallback(obj.communication_handles.client_expcontrol, ...
                "terminator",@obj.read_expcontrol_message);    

            disp('~~ BONSAI SERVER CREATED ~~')

        end

        function delete(obj)
            % Delete bonsai server, close connections

            % Get OSC handler and stop listeners
            obj.communication_handles.bonsai_receiver.stopListening();
            obj.communication_handles.bonsai_receiver.close();

            % Clear TCP client
            delete(obj.communication_handles.client_expcontrol)

        end

    end

    methods (Access = protected)

        function obj = read_expcontrol_message(obj,source,eventdata)

            incoming_message = readline(source);

            if strcmp(incoming_message, 'stop')
                % If stop, send STOP to bonsai
                plab.rig.bonsai_server_helpers.bonsai_oscsend( ...
                    obj.communication_handles.bonsai_udp,'/stop', ...
                    "localhost",30000,'s','stop');
            else
                % Otherwise, assume recording info and start bonsai
                obj = obj.run_bonsai(incoming_message);
            end
        end

        function obj = run_bonsai(obj,recording_info_json)

            % Get recording information
            recording_info = jsondecode(recording_info_json);

            % Copy Bonsai workflow folder into local recording path
            local_save_path = ...
                plab.locations.filename('local', ...
                recording_info.mouse,recording_info.date,recording_info.time,'bonsai');

            mkdir(local_save_path);

            [~,bonsai_name,bonsai_ext] = fileparts(recording_info.bonsai_filename);
            [~,bonsai_folder] = fileparts(fileparts(recording_info.bonsai_filename));

            local_bonsai_folder = fullfile(local_save_path,bonsai_folder);
            local_bonsai_filename = fullfile(local_bonsai_folder,[bonsai_name,bonsai_ext]);

            copyfile(fileparts(recording_info.bonsai_filename),local_bonsai_folder);

            % Create and start timer function for bonsai listener
            bonsai_timerfcn = timer('TimerFcn',{@obj.read_bonsai_message,local_save_path}, ...
                'Period',1,'ExecutionMode','fixedSpacing');
            start(bonsai_timerfcn)

            % Run Bonsai workflow (locally)
            command = ...
                sprintf("cd %s & ",fileparts(local_bonsai_filename)) +  ... % cd to Bonsai folder (to see extensions)
                sprintf("%s %s ","%LOCALAPPDATA%\Bonsai\Bonsai.exe", ... % Bonsai executable
                local_bonsai_filename) +  ... % Open local workflow
                " --start" + ... % Run workflow
                " --no-editor" + ... % No editor/command window
                sprintf(" -p %s=%s","SavePath",local_save_path); % Set save path variable

            system(command,'-echo');
           
        end

        function read_bonsai_message(obj,source,eventdata,local_save_path)

            incoming_message = obj.communication_handles.bonsai_listener.getMessageArgumentsAsString();
            
            % Recieved stop from Bonsai
            if strcmp(incoming_message,'stop')
                % Stop and delete timer function
                stop(source); 
                delete(source);

                % Send stop message to experiment controller
                writeline(obj.communication_handles.client_expcontrol, 'Bonsai finished');

                % Move data to server
                obj.move_data_to_server(local_save_path);
            end

        end

        function move_data_to_server(obj,local_save_path)
            % Move data from local to server

            % Check if the server is available
            if ~exist(plab.locations.server_data_path,'dir')
                warning('Server not accessible at %s',plab.locations.server_data_path)
                return
            end

            % Move local data directories to server
            curr_data_path_server = strrep(local_save_path, ...
                plab.locations.local_data_path,plab.locations.server_data_path);
            [status,message] = movefile(local_save_path,curr_data_path_server);
            if ~status
                warning('Failed copying to server: %s',message);
            else
            end

            % Clean local data folder (don't print results)
            plab.rig.clean_local_data_folder(false)

        end
    end
end