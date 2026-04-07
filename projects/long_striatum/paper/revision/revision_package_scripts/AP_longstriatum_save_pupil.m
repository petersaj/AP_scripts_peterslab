%% Load and package passive pupil data
% Modified from PG_pupillometryiwthplots

save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Marica_2026\data';

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

% Confidence threshold for SLEAP point scores
scoreThresh = 0.0;

% Define time to extract
mousecam_framerate = 30;
pupil_t = -0.5:1/mousecam_framerate:1.5;

% Initialize data cell
data_all = cell(length(animals),1);

% Loop across all mice
for animal_idx = 1:length(animals)

    animal = animals{animal_idx};
    data_animal = table;

    % Find passive recording days that also have task
    workflow_passive = {'lcr_passive'};
    recordings_passive = plab.find_recordings(animal, [], workflow_passive);
    workflow_task = {'stim_wheel_right*'};
    recordings_task = plab.find_recordings(animal, [], workflow_task);
    training_days = ismember({recordings_passive.day}, {recordings_task.day});
    train_rec_passive = recordings_passive(training_days);

    %% Load and process SLEAP data

    % Preallocate outputs. Cells stop lengths from being an issue.
    diameterPx_all   = cell(length(train_rec_passive),1);
    diameterZ_all    = cell(length(train_rec_passive),1);
    radius_all       = cell(length(train_rec_passive),1);
    center_all       = cell(length(train_rec_passive),1);
    fitRmse_all      = cell(length(train_rec_passive),1);

    % Big loop to get the above for each video
    for curr_rec = 1:length(train_rec_passive)

        % Get pupil SLEAP filename
        rec_day = train_rec_passive(curr_rec).day;
        rec_time = train_rec_passive(curr_rec).recording{end};

        pupil_sleap_dir = dir(plab.locations.filename('server',animal,rec_day,rec_time,'mousecam','sleap','pupil_v1','*.h5'));
        pupil_sleap_fn = fullfile(pupil_sleap_dir.folder,pupil_sleap_dir.name);

        % Read datasets
        tracks = h5read(pupil_sleap_fn,'/tracks');       % frames x nodes x 2
        pointScores = h5read(pupil_sleap_fn,'/point_scores')'; % frames x nodes
        instanceScores = h5read(pupil_sleap_fn, '/instance_scores')'; % transpose to 1 x frames

        sz = size(tracks);
        numFrames = sz(1);
        numNodes= sz(2);

        X = squeeze(tracks(:,:,1))'; % becomes [numFrames x numNodes]
        Y = squeeze(tracks(:,:,2))'; % [numFrames x numNodes]

        framesWithRawLabels = sum(any(~isnan(X),1));
        pct = (framesWithRawLabels/numFrames)*100;
        if pct < 50
            warning('%d/%d frames contain raw labels (%.2f%% of video)', ...
                framesWithRawLabels, numFrames, pct);
        else
            fprintf('%d/%d frames contain raw labels (%.2f%% of video)\n', ...
                framesWithRawLabels, numFrames, pct);
        end

        % Mask low-confidence instances
        X(:,  instanceScores< scoreThresh) = NaN;
        Y(:,  instanceScores< scoreThresh) = NaN;

        fprintf('%d bad instances dropped, %g3%% of video.', sum(sum(instanceScores<scoreThresh)), ((sum(sum(instanceScores<scoreThresh)))./numFrames)*100);

        % Mask low-confidence points
        X(pointScores < scoreThresh) = NaN;
        Y(pointScores < scoreThresh) = NaN;

        fprintf(' %d additional bad nodes dropped \n \n', sum(sum(pointScores<scoreThresh)));

        % Fit circles to frames
        pupil_radius = nan(numFrames,1);
        pupil_center = nan(numFrames,2);
        pupil_diameterPx = nan(numFrames,1);
        pupil_fitRmse = nan(numFrames,1);

        for f = 1:numFrames
            xpts = X(:,f);
            ypts = Y(:,f);
            valid = ~isnan(xpts) & ~isnan(ypts);
            if nnz(valid) < 3
                % Not enough points to fit a circle
                pupil_radius(f) = NaN;
                pupil_center(f,:) = [NaN NaN];
                pupil_diameterPx(f) = NaN;
                pupil_fitRmse(f) = NaN;
                continue;
            end

            xg = xpts(valid);
            yg = ypts(valid);

            % Could replace here with median difference b/t opposed points

            % Algebraic least-squares fit:
            % Solve for a,b,c in x^2 + y^2 + a*x + b*y + c = 0
            A = [xg, yg, ones(length(xg),1)];
            bvec = -(xg.^2 + yg.^2);

            % Solve linear system (least squares)
            p = A \ bvec;        % p = [a; b; c]
            a = p(1); bpar = p(2); c = p(3);

            xc = -a/2;
            yc = -bpar/2;
            radTerm = (a^2 + bpar^2)/4 - c;
            if radTerm <= 0
                % numerical degeneracy -> treat as invalid
                pupil_radius(f) = NaN;
                pupil_center(f,:) = [NaN NaN];
                pupil_diameterPx(f) = NaN;
                pupil_fitRmse(f) = NaN;
                continue;
            end
            R = sqrt(radTerm);

            % Store results
            pupil_radius(f) = R;
            pupil_center(f,:) = [xc yc];
            pupil_diameterPx(f) = 2*R;

            % Compute RMSE of radial residuals as a fit quality metric. Might have
            % to fit ovals instead if the fit is bad
            dists = hypot(xg - xc, yg - yc);
            residuals = dists - R;
            pupil_fitRmse(f) = sqrt(mean(residuals.^2));
        end

        % Now Z-score
        mu = nanmean(pupil_diameterPx);
        sig = nanstd(pupil_diameterPx);
        if sig == 0 || isnan(sig)
            pupil_diameterZ = nan(size(pupil_diameterPx));
        else
            pupil_diameterZ = (pupil_diameterPx - mu) ./ sig;
        end

        % Store results
        diameterPx_all{curr_rec} = pupil_diameterPx;
        diameterZ_all{curr_rec}  = pupil_diameterZ;
        radius_all{curr_rec}     = pupil_radius;
        center_all{curr_rec}     = pupil_center;
        fitRmse_all{curr_rec}    = pupil_fitRmse;

        % Lowpass and smooth pupil data
        pupil_diameterPx_filt = sgolayfilt(lowpass(fillmissing(pupil_diameterZ',"linear"),4,30),3,15);
        % (remove originally missing data)
        pupil_diameterPx_filt(isnan(pupil_diameterZ)) = NaN;


        %% Align pupillometry to stimuli
        load_parts.mousecam = true;
        ap.load_recording;

        % Get stim positions
        trial_stim_x_all = vertcat(trial_events.values.TrialStimX);
        trial_stim_x = trial_stim_x_all(1:length(stimOn_times));

        % Get quiescent trials
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');

        % Extract pupil traces around each stim time
        % (SLEAP frames can be < mousecam frames if empty data at end)
        extract_times = stimOn_times + pupil_t;
        pupil_diameter_trial = interp1(mousecam_times(1:length(pupil_diameterPx_filt)), ...
            pupil_diameterPx_filt,extract_times);

        % Save data in table
        data_animal.animal(curr_rec) = {animal};
        data_animal.rec_day(curr_rec) = {rec_day};

        data_animal.trial_stim_values(curr_rec) = {trial_stim_x(quiescent_trials)};
        data_animal.pupil_diameter(curr_rec) = {pupil_diameter_trial(quiescent_trials,:)};

        fprintf('Finished %s rec %d/%d\n',animal,curr_rec,length(train_rec_passive))
    end

    % Add current animal to full dataset
    data_all{animal_idx} = data_animal;

end

% Concatenate data into one table and save
pupil = vertcat(data_all{:});
save_filename = fullfile(save_path, 'pupil_passive');
save(save_filename, "pupil", "-v7.3");
fprintf('Saved: %s\n',save_filename);

