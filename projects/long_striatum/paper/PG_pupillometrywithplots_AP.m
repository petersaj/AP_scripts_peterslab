
%% Script to get pupil and plot diameter across trials from labelled SLEAP HDF5 files

% Directory for sorted SLEAP output HDF5 files
sleapDir = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Peter_Gorman\Pupils_temporary\Outputs_Sorted\lcr_passive';

% Confidence threshold for SLEAP point scores
scoreThresh = 0.0;    

% PSTH window
pre = 15;
post = 45;
x = -pre:post; 


%% Load behavior, add 'learned_days' and 'days_from_learning' fields
data_path = fullfile(plab.locations.server_path,'Lab','Papers','Marica_2026','data');
% Set stat and p-value to define learning day
use_stat = 'firstmove_mean';
learn_stat_p = 0.05;
% Load behavior
load(fullfile(data_path,'bhv'));
% Set "learned_days" and "days_from_learning"
bhv.learned_days = cellfun(@(x) x < learn_stat_p,bhv.(['stimwheel_pval_',use_stat]));
for curr_animal = unique(bhv.animal)'
    curr_rows = strcmp(bhv.animal,curr_animal);
    bhv.days_from_learning(curr_rows) = ...
        (1:sum(curr_rows))' - ...
        max([NaN,find(bhv.learned_days(curr_rows),1)]);
end

%% Grab the learning day for each animal
% See unique animals in bhv
aniUniq = unique(bhv.animal);
% initialize cell to store learning days
aniDates = cell(numel(aniUniq),2);
for curr_animal = 1:numel(aniUniq); 
    id = aniUniq{curr_animal};
    aniDates(curr_animal,1) = {id};
    % Get only recordings for this animal
    daysID = bhv(strcmp(bhv.animal, id),:);
    if any(daysID.days_from_learning == 0) % check if there was a learning day
        ld = daysID(daysID.days_from_learning == 0, :);
        aniDates(curr_animal,2) = ld.rec_day; % set learning date
    else
        aniDates(curr_animal,2) = {NaN}; % if no learning date, NaN
    end
end

%% Pupil responses to lcr_passive stimuli across days, for multiple mice

mouseIDs = unique(bhv.animal); % grabbing mice in question from bhv table

% pooled outputs 
psthData_all       = [];
orientationIdx_all = [];
learnDayIdx_all    = [];
mouseIdx_all       = [];    

% pooled outputs for quiescent trials only
psthData_allQui       = [];
orientationIdx_allQui = [];
learnDayIdx_allQui    = [];
mouseIdx_allQui       = []; 

% Loop across all mice
for m = 1:numel(mouseIDs)
    animalID = mouseIDs{m};
    fprintf(animalID); fprintf(' \n \n');
    dataDir = fullfile(sleapDir, animalID);  
    pattern = fullfile(dataDir, '*.analysis.h5');  
    files = dir(pattern);                          
    fileList = fullfile({files.folder}, {files.name})';
    nFiles = numel(fileList);
    
    % Preallocate outputs. Cells stop lengths from being an issue.
    diameterPx_all   = cell(nFiles,1);
    diameterZ_all    = cell(nFiles,1);
    radius_all       = cell(nFiles,1);
    center_all       = cell(nFiles,1);
    fitRmse_all      = cell(nFiles,1);
    
    % Big loop to get the above for each video
    
    for idx = 1:nFiles
    
        h5file = fileList{idx};
    
        % Read datasets
    
        tracks = h5read(h5file,'/tracks');       % frames x nodes x 2
        pointScores = h5read(h5file,'/point_scores')'; % frames x nodes
        instanceScores = h5read(h5file, '/instance_scores')'; % transpose to 1 x frames
    
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
        diameterPx_all{idx} = diameterPx;
        diameterZ_all{idx}  = diameterZ;
        radius_all{idx}     = radius;
        center_all{idx}     = center;
        fitRmse_all{idx}    = fitRmse;
    end

    %% Aligning pupillometry to stimuli
    % Follows on from processed pupil videos using hdf5 files (in same order as
    % below)
    
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

    %% Output to use
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
    
    nWindow = numel(-pre:post);
    nRecs = min(nFiles, size(frameStims,1)); 

    % Getting peri-stimulus pupil diameter values
    for this_rec = 1:nRecs
        % grab stimulus/frame info for this recording 
        frameTimes = frameStims{this_rec,1};
        stimTimes  = frameStims{this_rec,2};
        stimPosVec = frameStims{this_rec,3};    
        
        % pupil trace for this recording
        trace = pupilPerFile{this_rec}(:); % ensure column vector
        nFrames = numel(trace);
        
        % convert stim times to nearest (previous) frame indices
        stimFrameIdx = interp1(frameTimes, (1:numel(frameTimes))', stimTimes, 'previous', NaN);
    
        % if lengths differ for some reason, trim to the minimum and warn
        if numel(stimFrameIdx) ~= numel(stimPosVec)
            warning('Recording %d: stimTimes (%d) and stimPosVec (%d) differ in length — trimming to min length.', ...
                    this_rec, numel(stimFrameIdx), numel(stimPosVec));
            L = min(numel(stimFrameIdx), numel(stimPosVec));
            stimFrameIdx = stimFrameIdx(1:L);
            stimPosVec   = stimPosVec(1:L);
        end
    
        % build logical mask of valid trials 
        validTrials = ~isnan(stimFrameIdx) & ~isnan(stimPosVec);
        
        % filter to only valid trials
        stimFrameIdx = stimFrameIdx(validTrials);
        stimPosVec   = stimPosVec(validTrials);
        
        % also filter quiescent mask to match the same trial subset
        quiescent_trials = frameStims{this_rec,4};
        quiescent_trials = quiescent_trials(validTrials);
        
        if isempty(stimFrameIdx)
            continue;
        end
        
        % preallocate local trial matrix (trials x window)
        nTrialsLocal = numel(stimFrameIdx);
        localMat = nan(nTrialsLocal, nWindow);
        localOrient = nan(nTrialsLocal,1);
        
        winOffsets = -pre:post;
        for tt = 1:nTrialsLocal
            centerIdx = stimFrameIdx(tt);
            windowIdx = centerIdx + winOffsets;
            
            % clip and fill with NaN where outside bounds
            valid = windowIdx >= 1 & windowIdx <= nFrames;
            tmp = nan(1, nWindow);
            tmp(valid) = trace(windowIdx(valid));
            localMat(tt,:) = tmp;
            
            % map orientation label to index
            pos = stimPosVec(tt);
            if pos == -90
                localOrient(tt) = 1;    % left
            elseif pos == 0
                localOrient(tt) = 2;    % center
            elseif pos == 90
                localOrient(tt) = 3;    % right
            else
                % unknown orientation: add as NaN and store as a separate code
                localOrient(tt) = NaN;
            end
        end
        
        % append to global arrays
        psthData       = [psthData; localMat]; %#ok<AGROW>
        orientationIdx = [orientationIdx; localOrient]; %#ok<AGROW>
        dayIdx         = [dayIdx; repmat(this_rec, size(localMat,1), 1)]; %#ok<AGROW>
        
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
        
        % append to pooled arrays
        localMouse = repmat(m, nTrialsLocal,1);

        psthData_all       = [psthData_all; localMat];                %#ok<AGROW>
        orientationIdx_all = [orientationIdx_all; localOrient];       %#ok<AGROW>
        mouseIdx_all       = [mouseIdx_all; localMouse]; %#ok<AGROW>
       

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
        learnDayIdx_all = [learnDayIdx_all; localDay];

        localDayQui = localDay(frameStims{this_rec,4}',:);
        learnDayIdx_allQui = [learnDayIdx_allQui; localDayQui];
    end
end

%% Package into a beautiful table 

animalID = {};
for i = 1:numel(mouseIdx_all)
    animalID{i,1} = mouseIDs{mouseIdx_all(i)};
end

orientationDir = [];
for y = 1:numel(orientationIdx_all)
    if orientationIdx_all(y,1) == 1
            orientationDir(y,1) = -90;    % left
    elseif orientationIdx_all(y,1) == 2
            orientationDir(y,1) = 0;    % center
    elseif orientationIdx_all(y,1) == 3
            orientationDir(y,1) = 90;     % right
    end
end

pupilDiamByFrame = psthData_all;

pupils = table(animalID, mouseIdx_all, learnDayIdx_all, orientationDir, pupilDiamByFrame);

%% Package into a beautiful quiescent table 

animalIDQui = {};
for i = 1:numel(mouseIdx_allQui)
    animalIDQui{i,1} = mouseIDs{mouseIdx_allQui(i)};
end

orientationDirQui = [];
for y = 1:numel(orientationIdx_allQui)
    if orientationIdx_allQui(y,1) == 1
            orientationDirQui(y,1) = -90;    % left
    elseif orientationIdx_all(y,1) == 2
            orientationDirQui(y,1) = 0;    % center
    elseif orientationIdx_all(y,1) == 3
            orientationDirQui(y,1) = 90;     % right
    end
end

pupilDiamByFrameQui = psthData_allQui;

pupilsQuiescent = table(animalIDQui, mouseIdx_allQui, learnDayIdx_allQui, orientationDirQui, pupilDiamByFrameQui);

%% Script for videographic plots 

col = [0.066666666666667	0.266666666666667	0.611764705882353;  0.580392156862745	0.427450980392157	0.203921568627451; 0.890196078431372	0.309803921568627	0.309803921568627];

dayBins = {
    struct('mask', learnDayIdx_allQui <= -3,                         'title', '<=-3')
    struct('mask', learnDayIdx_allQui >= -2 & learnDayIdx_allQui <= -1, 'title', '-2:-1')
    struct('mask', learnDayIdx_allQui >= 0 & learnDayIdx_allQui <= 1, 'title', '0:1')
    struct('mask', learnDayIdx_allQui >= 2,                          'title', '>=2')
};


%% deriv with avg and sem across mice, no baseline subtraction

figure;
tiledlayout(1,4,'TileSpacing','compact','Padding','compact')

mouseIDs = unique(mouseIdx_allQui);

for g = 1:numel(dayBins)
    nexttile; hold on
    axis square

    for orientation = 1:3
        mouseTraces = [];   

        for m = 1:numel(mouseIDs)
            trialMask = orientationIdx_allQui == orientation & ...
                        dayBins{g}.mask & ...
                        mouseIdx_allQui == mouseIDs(m);

            subset = psthData_allQui(trialMask, :);

            % Drop sparse trials
            trialCoverage = mean(~isnan(subset), 2);
            subset = subset(trialCoverage >= 0.99, :);

            % Derivative per trial
            subsetDeriv = diff(subset, 1, 2);

            % Baseline subtract each trial
            % trialBase = mean(subsetDeriv(:, 20:21), 2, 'omitnan');
            % subsetDeriv = subsetDeriv - trialBase;

            % Mouse-level mean derivative trace
            mouseTrace = mean(subsetDeriv, 1, 'omitnan');
            mouseTraces = [mouseTraces; mouseTrace]; %#ok<AGROW>
        end

        % Mean and SEM across mice
        avg = mean(mouseTraces, 1, 'omitnan');
        nValid = sum(~isnan(mouseTraces), 1);
        sem = std(mouseTraces, 0, 1, 'omitnan') ./ sqrt(nValid);
        
        c = col(orientation,:); 
        % plot individual traces
        % for m = 1:size(mouseTraces,1)
        %     plot(mouseTraces(m,:), ...
        %         'Color', [c 0.4], ...   
        %         'LineWidth', 1);
        % end

        x = 1:numel(avg);
        upper = avg + sem;
        lower = avg - sem;

        fill([x fliplr(x)], [upper fliplr(lower)], c, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(x, avg, 'Color', c, 'LineWidth', 4)
    end

    ylim([-0.03, 0.035]);
    xline(15.5, 'k--')
    title(dayBins{g}.title)
    hold off
end


%% AUC quantification from derivative traces, averaged across mice, no baseline subtraction


basewindow = 15:16;   % pre-stimulus derivative bins
peakwindow = 25:45;  % AUC window on derivative trace

mouseIDs = unique(mouseIdx_allQui);

aucMouse = nan(numel(mouseIDs), numel(dayBins), 3);
aucMeans = nan(numel(dayBins), 3);
aucSems  = nan(numel(dayBins), 3);

for g = 1:numel(dayBins)
    for orientation = 1:3
        for m = 1:numel(mouseIDs)
            trialMask = orientationIdx_allQui == orientation & dayBins{g}.mask & mouseIdx_allQui == mouseIDs(m);

            subset = psthData_allQui(trialMask, :);

            % Drop sparse trials
            trialCoverage = mean(~isnan(subset), 2);
            subset = subset(trialCoverage >= 0.99, :);

            % Derivative per trial
            subsetDeriv = diff(subset, 1, 2);

            % Baseline subtract each trial 
            % trialBase = mean(subsetDeriv(:, basewindow), 2, 'omitnan');
            % subsetDeriv = subsetDeriv - trialBase;

            % Trial-wise AUC in the window
            trialAUC = trapz(peakwindow, subsetDeriv(:, peakwindow), 2);

            % Mouse-level mean AUC
            aucMouse(m, g, orientation) = mean(trialAUC, 'omitnan');
            %figure;hold on;for n = 1:size(subsetDeriv,1) plot(subsetDeriv(n,:));end
        end

        vals = aucMouse(:, g, orientation);
        aucMeans(g, orientation) = mean(vals, 'omitnan');
        nMouse = sum(~isnan(vals));
        aucSems(g, orientation) = std(vals, 0, 'omitnan') ./ sqrt(nMouse);
    end
end

% Plot across day bins 
figure; hold on
x = 1:numel(dayBins);

for orientation = 1:3
    y = aucMeans(:, orientation);
    e = aucSems(:, orientation);
    c = col(orientation,:);

    % plot individual mice
    % for g = 1:numel(dayBins)
    %     vals = aucMouse(:, g, orientation);
    %     xj = g + (rand(size(vals))-0.5)*0.15; % jitter
    % 
    %     scatter(xj, vals, 30, c, ...
    %         'filled', 'MarkerFaceAlpha', 0.5);
    % end

    errorbar(x, y, e, ...
        'Color', c, ...
        'LineWidth', 3, ...
        'CapSize', 0);
end

xlim([1 numel(dayBins)])
ylim([-0.3 0.3])  
set(gca, 'XTick', x, ...
         'XTickLabel', cellfun(@(s) s.title, dayBins, 'UniformOutput', false))
axis square
axis padded
xlabel('Day bin')
ylabel('Mean derivative AUC')

hold off


%% pairwise comparisons

for g = 1:numel(dayBins)
    vals = squeeze(aucMouse(:, g, :));  % [mouse x orientation]

    % remove mice with any NaNs
    valid = all(~isnan(vals), 2);
    vals = vals(valid, :);

    % pairwise comparisons
    [~, p12] = ttest(vals(:,1), vals(:,2));
    [~, p13] = ttest(vals(:,1), vals(:,3));
    [~, p23] = ttest(vals(:,2), vals(:,3));

    fprintf('Day %d: p12=%.3g, p13=%.3g, p23=%.3g\n', g, p12, p13, p23);
end