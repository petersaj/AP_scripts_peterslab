function trial_events = AP_load_bonsai_events(fn)
% Load Bonsai events file saved with CSVWRITER, output structure of events
%
% INPUTS: 
% fn = Bonsai csv file (prompts if nothing entered)
% 
% OUTPUTS: 
% trial_events = structure of events for each trial
%   .parameters = things set once at start (Trial 0)
%   .values = event values (as double)
%   .timestamps = event timestamps (as datetime)
% 
% Conventions: 
% Headers = Trial, Event, Value, Timestamp
% Trial 0 = Parameters set once at the workflow start
% Value = Only numbers (can't load mixed, e.g. [5;TRUE]) 

if ~exist('fn','var') || isempty(fn)
    fn = uigetfile('*.csv','Choose Bonsai file');
end

% Set Bonsai timestamp format
bonsai_table_opts = detectImportOptions(fn);
bonsai_table_opts = setvaropts(bonsai_table_opts,'Timestamp','Type','datetime', ...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSZ','TimeZone','local', ...
    'DatetimeFormat','yyyy-MM-dd HH:mm:ss.SSS');

% Load Bonsai CSV file
data = readtable(fn,bonsai_table_opts);

% Check for NaT timestamps, throw warning if any
if any(isnat(data.Timestamp))
    warning('Bonsai file ends improperly: %s',fn);
end

% Create nested structure for trial events
trial_events = struct('parameters',cell(1),'values',cell(1),'timestamps',cell(1)); 

% Save anything in "Trial 0" as a parameter
parameter_idx = data.Trial == 0;
unique_parameters = unique(data.Event(parameter_idx));
for curr_parameter = unique_parameters'
    curr_parameter_idx = parameter_idx & strcmp(data.Event,curr_parameter);
    trial_events.parameters.(cell2mat(curr_parameter)) = data.Value(curr_parameter_idx);
end

% Loop through trials (excluding 0), save values and timestamps for all events
% (exclude entries with empty events - happens sometimes on last entry)
empty_events = cellfun(@isempty,data.Event);
n_trials = max(data.Trial);
unique_events = unique(data.Event(~parameter_idx & ~empty_events));
for curr_trial = 1:n_trials
    for curr_event = unique_events'
        curr_event_idx = data.Trial == curr_trial & strcmp(data.Event,curr_event);
        trial_events.values(curr_trial).(cell2mat(curr_event)) = data.Value(curr_event_idx);
        trial_events.timestamps(curr_trial).(cell2mat(curr_event)) = data.Timestamp(curr_event_idx);
    end
end





