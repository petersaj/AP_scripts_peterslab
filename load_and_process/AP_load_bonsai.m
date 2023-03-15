function trial_events = AP_load_bonsai(fn)

if ~exist('fn','var') || isempty(fn)
    fn = uigetfile('*.csv','Choose Bonsai file');
end

% Load Bonsai CSV file
data = readtable(fn,'Delimiter',',');

% Convert Bonsai timestamp to matlab timestamp
bonsai_timestamp_format = 'yyyy-MM-dd''T''HH:mm:ss.SSS';
bonsai_timestamp_length = length(bonsai_timestamp_format)-2; % -2 for the 'T'

matlab_timetamps = cellfun(@(x) datetime(x(1:bonsai_timestamp_length), ...
    'InputFormat',bonsai_timestamp_format),data.Timestamp);

% Create nested structure for trial events
trial_events = struct('parameters',cell(1),'values',cell(1),'timestamps',cell(1)); 

% Save anything in "Trial 0" as a parameter
parameter_idx = data.Trial == 0;
unique_parameters = unique(data.Event(parameter_idx));
for curr_parameter = unique_parameters'
    curr_parameter_idx = parameter_idx & strcmp(data.Event,curr_parameter);
    trial_events.parameters.(cell2mat(curr_parameter)) = data.Value(curr_parameter_idx);
end

% Loop through trials, save values and timestamps for all events
n_trials = max(data.Trial);
unique_events = unique(data.Event);
for curr_trial = 1:n_trials
    for curr_event = unique_events'
        curr_event_idx = data.Trial == curr_trial & strcmp(data.Event,curr_event);
        trial_events.values(curr_trial).(cell2mat(curr_event)) = data.Value(curr_event_idx);
        trial_events.timestamps(curr_trial).(cell2mat(curr_event)) = matlab_timetamps(curr_event_idx);
    end
end















