
%% Script to get pupil and plot diameter across trials from labelled SLEAP HDF5 files

% Directory for sorted SLEAP output HDF5 files
sleapDir = fullfile(plab.locations.server_path,'Users\Peter_Gorman\Pupils_temporary\Outputs_Sorted\lcr_passive');

% Confidence threshold for SLEAP point scores
scoreThresh = 0.0;    

% Define time to extract
mousecam_framerate = 30;
pupil_t = -0.5:1/mousecam_framerate:1.5;

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
        localMouse = repmat(m, nTrialsLocal,1);

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


%% Package into a beautiful quiescent table 

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


%% Script for videographic plots 

col = [0.38823529411	0.38823529411	0.38823529411;  0.066666666666667	0.266666666666667	0.611764705882353; 0.890196078431372	0.309803921568627	0.309803921568627];

t = pupil_t;

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

tDeriv = t(1:end-1) + mean(diff(t))/2;

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
            subset = subset(trialCoverage == 1, :);

            % Derivative per trial
            subsetDeriv = diff(subset, 1, 2);

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
        %     plot(tDeriv, mouseTraces(m,:), 'Color', 0.7*c + 0.3, 'LineWidth', 1);
        % end
        
        ap.errorfill(tDeriv, avg, sem, c, 0.2, 1, 4);
    end

    ylim([-0.03, 0.035]);
    xline(0)
    title(dayBins{g}.title)
    hold off
end


%% AUC quantification from derivative traces, averaged across mice, no baseline subtraction

t = pupil_t;
tDeriv = t(1:end-1) + mean(diff(t))/2;

baseMask = tDeriv >= -0.02 & tDeriv <= 0.02;   
peakMask = tDeriv >= 0.32 & tDeriv <= 0.98;    

mouseIDs = unique(mouseIdx_allQui);

aucMouse = nan(numel(mouseIDs), numel(dayBins), 3);
aucMeans = nan(numel(dayBins), 3);
aucSems  = nan(numel(dayBins), 3);

for g = 1:numel(dayBins)
    for orientation = 1:3
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

            % Trial-wise AUC in the window
            trialAUC = trapz(tDeriv(peakMask), subsetDeriv(:, peakMask), 2);

            % Mouse-level mean AUC
            aucMouse(m, g, orientation) = mean(trialAUC, 'omitnan');
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


%% permutation tests: orientation comparisons within each day bin

nPerm = 10000;

oriPairs = [1 2;
            1 3;
            2 3];

oriResults = table();

for g = 1:numel(dayBins)
    vals = squeeze(aucMouse(:, g, :));   % [mouse x orientation]

    for p = 1:size(oriPairs,1)
        o1 = oriPairs(p,1);
        o2 = oriPairs(p,2);

        x = vals(:, o1);
        y = vals(:, o2);

        [pval, obsDiff, nPairs] = permtest_paired_signflip(x, y, nPerm);

        newRow = table( ...
            g, string(dayBins{g}.title), o1, o2, nPairs, obsDiff, pval, ...
            'VariableNames', {'dayIdx','dayLabel','ori1','ori2','nPairs','obsDiff','pValue'} );

        oriResults = [oriResults; newRow]; %#ok<AGROW>
    end
end

disp(oriResults)

%% permutation tests: consecutive day comparisons across all mice

aucMouseDay = squeeze(mean(aucMouse, 3, 'omitnan'));   % [mouse x day]

nPerm = 10000;

dayConsecResults = table();

for g = 1:numel(dayBins)-1
    x = aucMouseDay(:, g);
    y = aucMouseDay(:, g+1);

    [pval, obsDiff, nPairs] = permtest_paired_signflip(x, y, nPerm);

    newRow = table( ...
        g, string(dayBins{g}.title), ...
        g+1, string(dayBins{g+1}.title), ...
        nPairs, obsDiff, pval, ...
        'VariableNames', {'day1Idx','day1Label','day2Idx','day2Label','nPairs','obsDiff','pValue'} );

    dayConsecResults = [dayConsecResults; newRow]; %#ok<AGROW>
end

disp(dayConsecResults)

%% permutation tests: day comparisons within each orientation

nPerm = 10000;

dayResults = table();

for orientation = 1:3
    for g1 = 1:numel(dayBins)-1
        for g2 = g1+1:numel(dayBins)

            x = aucMouse(:, g1, orientation);
            y = aucMouse(:, g2, orientation);

            [pval, obsDiff, nPairs] = permtest_paired_signflip(x, y, nPerm);

            newRow = table( ...
                orientation, g1, string(dayBins{g1}.title), ...
                g2, string(dayBins{g2}.title), ...
                nPairs, obsDiff, pval, ...
                'VariableNames', {'orientation','day1Idx','day1Label','day2Idx','day2Label','nPairs','obsDiff','pValue'} );

            dayResults = [dayResults; newRow]; %#ok<AGROW>
        end
    end
end

disp(dayResults)

oriResults.pBonf = min(oriResults.pValue * 3, 1);
dayConsecResults.pBonf = min(dayConsecResults.pValue * (numel(dayBins)-1), 1);

function [p, obsDiff, nPairs] = permtest_paired_signflip(x, y, nPerm)
%PERMTEST_PAIRED_SIGNFLIP Paired permutation test using sign flips.
%   x, y      paired observations (can contain NaNs)
%   nPerm     number of permutations
%
%   p         two-sided permutation p-value
%   obsDiff   observed mean paired difference
%   nPairs    number of paired observations used

    if nargin < 3 || isempty(nPerm)
        nPerm = 10000;
    end

    % keep only complete pairs
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid);
    y = y(valid);

    d = x(:) - y(:);
    nPairs = numel(d);

    if nPairs == 0
        p = NaN;
        obsDiff = NaN;
        return
    end

    obsDiff = mean(d);

    permStats = nan(nPerm,1);
    for i = 1:nPerm
        signs = sign(rand(nPairs,1) - 0.5);
        signs(signs == 0) = 1;
        permStats(i) = mean(d .* signs);
    end

    p = mean(abs(permStats) >= abs(obsDiff));
end