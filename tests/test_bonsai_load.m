fn = "C:\Users\petersa\Desktop\test.csv";

data = readtable(fn,'Delimiter',',');

bonsai_timestamp_format = 'yyyy-MM-dd''T''HH:mm:ss.SSS';
bonsai_timestamp_length = length(bonsai_timestamp_format)-2; % -2 for the 'T'

matlab_timetamps = cellfun(@(x) datetime(x(1:bonsai_timestamp_length), ...
    'InputFormat',bonsai_timestamp_format),data.Timestamp);


% in progress: make nested structure grouped by trial
n_trials = max(data.Trial);
unique_events = unique(data.Event);

trial_events = cell2struct(cell(length(unique_events),n_trials),unique_events);

for curr_trial = 1:n_trials
    for curr_event = reshape(unique_events,1,[])
        curr_event_idx = data.Trial == curr_trial & strcmp(data.Event,curr_event);
        trial_events(curr_trial).(cell2mat(curr_event)).Value = data.Value(curr_event_idx);
        trial_events(curr_trial).(cell2mat(curr_event)).Timestamp = matlab_timetamps(curr_event_idx);
    end
end


% trial_events = cell2struct(cell(2,n_trials),{'value','time'});
% 
% for curr_trial = 1:n_trials
%     for curr_event = reshape(unique_events,1,[])
%         curr_event_idx = data.Trial == curr_trial & strcmp(data.Event,curr_event);
%         trial_events(curr_trial).value.(cell2mat(curr_event)) = data.Value(curr_event_idx);
%         trial_events(curr_trial).time.(cell2mat(curr_event)) = matlab_timetamps(curr_event_idx);
%     end
% end


% data_struct = cell2struct(cell(length(unique_events)*2,n_trials), ...
%     [append(unique_events,'_value'),append(unique_events,'_time')]);
% 
% for curr_trial = 1:n_trials
%     for curr_event = reshape(unique_events,1,[])
%         curr_event_idx = data.Trial == curr_trial & strcmp(data.Event,curr_event);
%         data_struct(curr_trial).(append(cell2mat(curr_event),'_value')) = data.Value(curr_event_idx);
%         data_struct(curr_trial).(append(cell2mat(curr_event),'_time')) = matlab_timetamps(curr_event_idx);
%     end
% end














