% Modified from PG_pupillometryiwthplots

save_path = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\Papers\Marica_2025\data';

animals = { ...
    'AM011','AM012','AM014','AM015','AM016','AM017', ...
    'AM018','AM019','AM021','AM022','AM026','AM029', ...
    'AP023','AP025'};

% Confidence threshold for SLEAP point scores
scoreThresh = 0.0;

% Define time to extract
mousecam_framerate = 30;
pupil_t = -0.5:1/mousecam_framerate:1.5;

% pooled outputs for quiescent trials only
psthData_allQui       = [];
orientationIdx_allQui = [];
learnDayIdx_allQui    = [];
mouseIdx_allQui       = [];

% Loop across all mice
for animal_idx = 1:length(animals)

    animal = animals{animal_idx};

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
        radius = nan(numFrames,1);
        center = nan(numFrames,2);
        diameterPx = nan(numFrames,1);
        fitRmse = nan(numFrames,1);

        for f = 1:numFrames
            xpts = X(:,f);
            ypts = Y(:,f);
            valid = ~isnan(xpts) & ~isnan(ypts);
            if nnz(valid) < 3
                % Not enough points to fit a circle
                radius(f) = NaN;
                center(f,:) = [NaN NaN];
                diameterPx(f) = NaN;
                fitRmse(f) = NaN;
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
                radius(f) = NaN;
                center(f,:) = [NaN NaN];
                diameterPx(f) = NaN;
                fitRmse(f) = NaN;
                continue;
            end
            R = sqrt(radTerm);

            % Store results
            radius(f) = R;
            center(f,:) = [xc yc];
            diameterPx(f) = 2*R;

            % Compute RMSE of radial residuals as a fit quality metric. Might have
            % to fit ovals instead if the fit is bad
            dists = hypot(xg - xc, yg - yc);
            residuals = dists - R;
            fitRmse(f) = sqrt(mean(residuals.^2));
        end

        % Now Z-score
        mu = nanmean(diameterPx);
        sig = nanstd(diameterPx);
        if sig == 0 || isnan(sig)
            diameterZ = nan(size(diameterPx));
        else
            diameterZ = (diameterPx - mu) ./ sig;
        end

        % Store results
        diameterPx_all{curr_rec} = diameterPx;
        diameterZ_all{curr_rec}  = diameterZ;
        radius_all{curr_rec}     = radius;
        center_all{curr_rec}     = center;
        fitRmse_all{curr_rec}    = fitRmse;

        % Lowpass and smooth pupil data
        diameterPx_filt = sgolayfilt(lowpass(fillmissing(diameterZ',"linear"),4,30),3,15);
        % (remove originally missing data)
        diameterPx_filt(isnan(diameterZ)) = NaN;


        %% Align pupillometry to stimuli
        load_parts.mousecam = true;
        ap.load_recording;




    
    diameterPxAllFlip = cellfun(@transpose, diameterPx_all, 'UniformOutput', false);
    diameterZAllFlip = cellfun(@transpose, diameterZ_all, 'UniformOutput', false);

    %Lowpass filter to get rid of jitteriness
    diameterZAllFlipFilt = cellfun(@(x) lowpass(fillmissing(x, "linear"), 4, 30), diameterZAllFlip, 'UniformOutput', false);
    diameterPxAllFlipFilt = cellfun(@(x) lowpass(fillmissing(x, "linear"), 4, 30), diameterPxAllFlip, 'UniformOutput', false);

    % try savgol after
    diameterZAllFlipFiltSav = cellfun(@(x) sgolayfilt(x, 3, 15), diameterZAllFlipFilt, 'UniformOutput', false);
    diameterPxAllFlipFiltSav = cellfun(@(x) sgolayfilt(x, 3, 15), diameterPxAllFlipFilt, 'UniformOutput', false);

    % for just savgol
    % diameterZAllFlipSav = cellfun(@(x) sgolayfilt(x, 3, 15), diameterZAllFlip, 'UniformOutput', false);
    % diameterPxAllFlipSav = cellfun(@(x) sgolayfilt(x, 3, 15), diameterPxAllFlip, 'UniformOutput', false);

    % Output to use
    pupilPerFile = diameterZAllFlipFiltSav;

    % Mask fillmissing'ed values back out
    basemasks = cellfun(@isnan, diameterPx_all, 'UniformOutput',false);
    for nanmask = 1:numel(pupilPerFile)
        pupilPerFile{nanmask}(basemasks{nanmask}) = NaN;
    end





    end

    %% Align pupillometry to stimuli

    % Grab the recording days and time for the mouse

    allPassives = cell(numel(fileList),2);

    passiveM = plab.find_recordings(animalID, [], 'lcr_passive');

    % limit to recording days, not habituation
    if nFiles ~= size(passiveM, 2)
        fprintf(animalID, ' mismatch between recording files and hdf5 files. Excluding habitation days may be broken')
    end

    taskIdx = find(strcmp(animalID,bhv{:,"animal"}));
    taskDays = bhv{:,"rec_day"}(taskIdx);
    taskMask = ismember({passiveM.day}, {taskDays{:}});
    passiveM = passiveM(taskMask);

    diameterPx_all = diameterPx_all(taskMask);
    diameterZ_all  = diameterZ_all(taskMask);

    [s,recs] = size(passiveM);

    % loop to get relevant variables from each recording
    frameStims = cell(recs,5);
    for curr_rec = 1:recs
        animal = passiveM(curr_rec).animal;
        rec_day = passiveM(curr_rec).day;
        rec_time = passiveM(curr_rec).recording{end};
        verbose = false; % this turns on/off progress display in command line

        % Loading components separately. Could maybe be more efficient but is
        % better than calling ap.load_recording
        ap.load_timelite;
        ap.load_mousecam;
        ap.load_bonsai;

        % Get quiescence mask
        stim_window = [0,0.5];
        quiescent_trials = arrayfun(@(x) ~any(wheel_move(...
            timelite.timestamps >= stimOn_times(x)+stim_window(1) & ...
            timelite.timestamps <= stimOn_times(x)+stim_window(2))), ...
            (1:length(stimOn_times))');

        % Grab the frame times and put them in one cell
        frameStims(curr_rec,1) = {mousecam_times};
        frameStims(curr_rec,2) = {stimOn_times'};
        frameStims(curr_rec,3) = {(vertcat(trial_events.values.TrialStimX))'};
        frameStims(curr_rec,4) = {quiescent_trials'};
    end




    diameterPxAllFlip = cellfun(@transpose, diameterPx_all, 'UniformOutput', false);
    diameterZAllFlip = cellfun(@transpose, diameterZ_all, 'UniformOutput', false);

    %Lowpass filter to get rid of jitteriness
    diameterZAllFlipFilt = cellfun(@(x) lowpass(fillmissing(x, "linear"), 4, 30), diameterZAllFlip, 'UniformOutput', false);
    diameterPxAllFlipFilt = cellfun(@(x) lowpass(fillmissing(x, "linear"), 4, 30), diameterPxAllFlip, 'UniformOutput', false);

    % try savgol after
    diameterZAllFlipFiltSav = cellfun(@(x) sgolayfilt(x, 3, 15), diameterZAllFlipFilt, 'UniformOutput', false);
    diameterPxAllFlipFiltSav = cellfun(@(x) sgolayfilt(x, 3, 15), diameterPxAllFlipFilt, 'UniformOutput', false);

    % for just savgol
    % diameterZAllFlipSav = cellfun(@(x) sgolayfilt(x, 3, 15), diameterZAllFlip, 'UniformOutput', false);
    % diameterPxAllFlipSav = cellfun(@(x) sgolayfilt(x, 3, 15), diameterPxAllFlip, 'UniformOutput', false);

    % Output to use
    pupilPerFile = diameterZAllFlipFiltSav;

    % Mask fillmissing'ed values back out
    basemasks = cellfun(@isnan, diameterPx_all, 'UniformOutput',false);
    for nanmask = 1:numel(pupilPerFile)
        pupilPerFile{nanmask}(basemasks{nanmask}) = NaN;
    end

    % Prepare accumulators
    psthData       = [];   % will become [nTrials x (pre+post+1)]
    orientationIdx = [];   % will become [nTrials x 1] (1=left,2=center,3=right)
    dayIdx         = [];   % will become [nTrials x 1] (recording index n)
    orientationLab = {};
    recDayLabels   = {};

    nRecs = min(nFiles, size(frameStims,1));

    % Getting peri-stimulus pupil diameter values
    for this_rec = 1:nRecs
        % grab stimulus/frame info for this recording
        frameTimes = frameStims{this_rec,1}';
        stimTimes  = frameStims{this_rec,2}';
        stimPosVec = frameStims{this_rec,3}';
        quiescent_trials = frameStims{this_rec,4}';

        % pupil trace for this recording
        trace = pupilPerFile{this_rec}(:); % ensure column vector
        nFrames = numel(trace);

        % only take positions with a stimtime and frames with a pupil
        stimPosVec = stimPosVec(1:numel(stimTimes));
        frameTimes = frameTimes(1:numel(trace));

        % filter to only valid trials
        validTrials = ~isnan(stimTimes) & ~isnan(stimPosVec);
        stimTimes = stimTimes(validTrials);
        stimPosVec = stimPosVec(validTrials);
        quiescent_trials = quiescent_trials(validTrials);

        % Extract pupil traces around each stim time
        extract_times = stimTimes + pupil_t;
        localMat = interp1(frameTimes, trace, extract_times);

        % Orientation labels
        localOrient = nan(numel(stimPosVec), 1);
        localOrient(stimPosVec == -90) = 1;   % left
        localOrient(stimPosVec == 0)   = 2;   % center
        localOrient(stimPosVec == 90)  = 3;   % right

        % readable labels for trialInfo
        for tt = 1:numel(localOrient)
            if localOrient(tt) == 1, ol = 'left';
            elseif localOrient(tt) == 2, ol = 'center';
            elseif localOrient(tt) == 3, ol = 'right';
            else ol = 'other'; end
            orientationLab{end+1,1} = ol; %#ok<SAGROW>
            % store the rec day string if passiveM has it
            if isfield(passiveM, 'day') && numel(passiveM) >= this_rec
                recDayLabels{end+1,1} = passiveM(this_rec).day; %#ok<SAGROW>
            else
                recDayLabels{end+1,1} = sprintf('rec%d', this_rec); %#ok<SAGROW>
            end
        end

        % Grab the learning day from aniDates to produce LearnDayIdx
        AL = find(strcmp(animalID, aniDates(:,1)));
        if ~isnan(aniDates{AL,2})
            learnDay = datetime(aniDates{AL,2});
            idxLearn = find({passiveM.day} == learnDay);
            learnDayIdx = dayIdx - idxLearn;
        else
            learnDayIdx = nan([size(dayIdx)]);
            idxLearn = NaN;
        end

        % Now to make the megatrial matrix
        nTrialsLocal = size(localMat,1);
        localMouse = repmat(animal_idx, nTrialsLocal,1);

        % append for for quiescent
        localMatQui    = localMat(quiescent_trials,:);
        localOrientQui = localOrient(quiescent_trials,:);
        localMouseQui  = localMouse(quiescent_trials,:);

        psthData_allQui       = [psthData_allQui; localMatQui];                %#ok<AGROW>
        orientationIdx_allQui = [orientationIdx_allQui; localOrientQui];       %#ok<AGROW>
        mouseIdx_allQui       = [mouseIdx_allQui; localMouseQui]; %#ok<AGROW>

        % learn-day relative index for this recording (n is recording index)
        if exist('idxLearn','var') && ~isempty(idxLearn) && ~isnan(idxLearn)
            rel = this_rec - idxLearn;
        else
            rel = NaN;
        end
        localDay = repmat(rel, nTrialsLocal,1);
        localDayQui = localDay(frameStims{this_rec,4}',:);
        learnDayIdx_allQui = [learnDayIdx_allQui; localDayQui];
    end
end


% Package into a beautiful quiescent table

animalIDQui = {};
for i = 1:numel(mouseIdx_allQui)
    animalIDQui{i,1} = mouseIDs{mouseIdx_allQui(i)};
end

orientationDirQui = [];
for y = 1:numel(orientationIdx_allQui)
    if orientationIdx_allQui(y,1) == 1
        orientationDirQui(y,1) = -90;    % left
    elseif orientationIdx_allQui(y,1) == 2
        orientationDirQui(y,1) = 0;    % center
    elseif orientationIdx_allQui(y,1) == 3
        orientationDirQui(y,1) = 90;     % right
    end
end

pupilDiamByFrameQui = psthData_allQui;

pupilsQuiescent = table(animalIDQui, mouseIdx_allQui, learnDayIdx_allQui, orientationDirQui, pupilDiamByFrameQui);


