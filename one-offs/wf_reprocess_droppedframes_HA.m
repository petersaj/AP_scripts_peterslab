 
%% Re-SVD HA mice after specific date
% (this was week during frame drop detection with bug)

clear

animals = {'HA009','HA010','HA011','HA012'};
bad_day = '2025-05-20';

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    recordings = plab.find_recordings(animal);

    check_recordings = datenum({recordings.day}) >= datenum(bad_day) & ...
        datenum({recordings.day}) < today;

    for curr_recording_idx = find(check_recordings)'

        rec_day = recordings(curr_recording_idx).day;

        dropped_frame_flag = false;
        for curr_rec_time_idx = 1:length(recordings(curr_recording_idx).recording)

            curr_rec_time = recordings(curr_recording_idx).recording{curr_rec_time_idx};

            widefield_metadata_fn = ...
                plab.locations.filename('server',animal,rec_day,[], ...
                'widefield',sprintf('widefield_%s_metadata.bin',curr_rec_time));

            % (find dropped frames with the bad and fixed method, flag
            % either for re-SVD)
            dropped_frames_bad = find_dropped_frames_bad(widefield_metadata_fn,false);
            dropped_frames_fixed = plab.wf.find_dropped_frames(widefield_metadata_fn,false);

            if ~isempty(dropped_frames_bad) || ~isempty(dropped_frames_fixed)
                dropped_frame_flag = true;
            end
        end

        % If any dropped frames in the day, re-svd
        if dropped_frame_flag
            fprintf('Re-SVDing: %s %s\n',animal,rec_day);

            % Copy widefield locally
            server_wf_path = plab.locations.filename('server',animal,rec_day,[],'widefield');
            local_wf_path = strrep(server_wf_path,plab.locations.server_data_path,plab.locations.local_data_path);
            copyfile(server_wf_path,local_wf_path);

            % Re-SVD data
            plab.wf.batch_widefield_preprocess;

        end
    end
end


%% Find and plot dropped frames by day from set of mice 

clear 

animals = {'HA009','HA010','HA011','HA012'};

dropped_frame_n = cell(size(animals));
dropped_frame_days = cell(size(animals));

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    recordings = plab.find_recordings(animal);
    check_recordings = find(datetime({recordings.day}) < datetime('today'));

    dropped_frame_n{curr_animal} = zeros(size(check_recordings));
    dropped_frame_days{curr_animal} = {recordings(check_recordings).day};

    for curr_recording_idx = check_recordings

        rec_day = recordings(curr_recording_idx).day;

        for curr_rec_time_idx = 1:length(recordings(curr_recording_idx).recording)

            curr_rec_time = recordings(curr_recording_idx).recording{curr_rec_time_idx};

            widefield_metadata_fn = ...
                plab.locations.filename('server',animal,rec_day,[], ...
                'widefield',sprintf('widefield_%s_metadata.bin',curr_rec_time));

            % dropped_frames = find_dropped_frames_bad(widefield_metadata_fn,false);
            dropped_frames = plab.wf.find_dropped_frames(widefield_metadata_fn,false);
            
            if ~isempty(dropped_frames)
                dropped_frame_n{curr_animal}(curr_recording_idx) = ...
                    dropped_frame_n{curr_animal}(curr_recording_idx)+length(dropped_frames);
            end
        end
    end
end

figure;
nexttile; hold on
arrayfun(@(x) plot(datetime(dropped_frame_days{x}),dropped_frame_n{x}), ...
    1:length(animals));
ylabel('Dropped frames');
xline(datetime('2025-05-21'),'r','Sophos off');
legend(animals);

nexttile;
[x,x_grp] = ap.groupfun(@sum,horzcat(dropped_frame_n{:}), ...
    datenum(horzcat(dropped_frame_days{:})));
plot(datetime(x_grp,'ConvertFrom','datenum'),x,'k','linewidth',2);
ylabel('Dropped frames');
xline(datetime('2025-05-21'),'r','Sophos off');
xline(datetime('2025-05-30'),'r','Sophos internet-only');

ap.prettyfig;