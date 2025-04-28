%% EEG Experiment

% Housekeeping
clear; close all

%% Collect subject ID
cfg.subjectID = input('Enter subject number: ', 's');
cfg.date = datetime('today');

%% General Experiment Configuration
testrun = false;

cfg.debugTimingFactor = 0.2; % Must be 1 for accurat timing (< 1 will give faster timing)
cfg.imageDuration = 0.25 * cfg.debugTimingFactor;  % Image presentation time in seconds
cfg.iti = 0.75 * cfg.debugTimingFactor;  % Inter-trial interval in seconds
cfg.startPad = 2 * cfg.debugTimingFactor;  % Time before the first trial in seconds
cfg.endPad = 2 * cfg.debugTimingFactor;  % Time after the last trial in seconds
cfg.ImageFileFormat = 'tif';

numBlocks = 2;
categories = {'kitchen', 'bathroom'};

%% Paths
stimPath = fullfile(pwd,'..', 'stimuli');
outputPath = fullfile(pwd,'..', 'sourcedata', ['sub-', cfg.subjectID], 'beh');
functionPath = fullfile(pwd,'utilities');

% add functions folder to path
addpath(functionPath)

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

%% Image Loading
kitchenImages = dir(fullfile(stimPath, 'kitchen', ['*.', cfg.ImageFileFormat]));
bathroomImages = dir(fullfile(stimPath, 'bathroom', ['*.', cfg.ImageFileFormat]));

if testrun
    kitchenImages = kitchenImages(1:10);
    bathroomImages = bathroomImages(1:10);
end
numKitchenTrials = length(kitchenImages)*2; % ?
numBathroomTrials = length(bathroomImages)*2;

% %% Help File Setup
% helpPath = fullfile(pwd,'..','helpdata');
% runHelpFile = fullfile(helpPath, 'permutations.tsv');
% % Make sure helpPath exists
% if ~exist(helpPath, 'dir')
%     mkdir(helpPath);
% end
% % Check if file name exists already to avoid overwriting 
% if exist(runHelpFile, 'file')
%     error('File name exists already')
% end 
% helpfileID = fopen(runHelpFile, 'w');
% 
% % Write header for the log file
% fprintf(helpfileID, 'block\timgnumber\tcategory\timagename\n');

%% Output File Setup
% BIDS-compliant log file
runOutputFile = fullfile(outputPath, sprintf('sub-%s_task-main_events.tsv', cfg.subjectID));

% Check if file name exists already to avoid overwriting 
if exist(runOutputFile, 'file')
    error('File name exists already')
end 
fileID = fopen(runOutputFile, 'w');

% Write header for the log file
fprintf(fileID, 'subject\tblock\ttrial\ttexture\tcategory\timage\ttrialOnset\titiOnset\ttrialEnd\ttriggerDate\ttriggerTimeStamp\tOnsets\tresponseKey\tresponseTime\taccuracy\n');

%% Initialize Psychtoolbox
Screen('Preference', 'SkipSyncTests', 1);  % Skip sync tests for demo purposes (remove in actual experiment)
PsychDefaultSetup(2)

% Define colors
Color.white = [255 255 255]; Color.black = [0 0 0]; Color.gray=(Color.black+Color.white)/2;
Color.red = [255 0 0]; Color.yellow= [255 255 0];

screenNumber = max(Screen('Screens'));
screencount=size(Screen('screens'),2);
if screencount>1
    windowrect=Screen(1,'rect');
    screenNumber=1; %%%%%%%%%%%%%%% 2
else
    windowrect=Screen(0,'rect');
    screenNumber=0;
end

%[window, windowRect] = Screen('OpenWindow', screenNumber, [Color.gray], [], 32, 2,[], [],  kPsychNeed32BPCFloat);
[window, windowRect] = Screen('OpenWindow', screenNumber, [Color.gray]); % für Mac
[xCenter, yCenter] = RectCenter(windowRect);
ifi = Screen('GetFlipInterval', window);  % Get refresh interval
Screen('TextSize', window, 40) % define font size
DrawFormattedText(window, 'Loading...', 'center', 'center', [0 0 0]);
Screen('Flip', window);  % Show fixation cross
HideCursor; % Hide mouse cursor

%% Load all images as textures
%bathroomTextures = cell(1, length(bathroomImages));
for i = 1:length(bathroomImages)
    bathroomTextures(i) = Screen('MakeTexture', window, imread(fullfile(bathroomImages(i).folder, bathroomImages(i).name)));
end

%kitchenTextures = cell(1, length(kitchenImages));
for i = 1:length(kitchenImages)
    kitchenTextures(i) = Screen('MakeTexture', window, imread(fullfile(kitchenImages(i).folder, kitchenImages(i).name)));
end

%% Calculate size for desired degree of visual angle

% define visual angle
x_degree = 8;
y_degree = 6;

% Get the screen resolution and viewing distance
viewing_dist = 50;
screen_hor=31; %Convert to cm (horizontal)
screen_vert=22; %Convert to cm (vertial)

% Calculate pixels per centimeter
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
pixPerCmX = screenXpixels / screen_hor;
pixPerCmY = screenYpixels / screen_vert;

% Calculate the size in cm for the given visual angles
sizeCmX = 2 * viewing_dist * tan(deg2rad(x_degree) / 2);
sizeCmY = 2 * viewing_dist * tan(deg2rad(y_degree) / 2);

% Convert the size from cm to pixels
sizePixX = round(sizeCmX * pixPerCmX);
sizePixY = round(sizeCmY * pixPerCmY);

% get rectangle for image of correct size
image_rect = CenterRectOnPointd([0 0 sizePixX sizePixY], xCenter, yCenter);

%% Initialize keyboard
KbName('UnifyKeyNames');
abortKey = KbName('ESCAPE');
% define target keys
presentKey = 79;%37;
absentKey = 80;%39;

%% Initialize overall accuracy
overallAccuracy = [];

%% start for-loop
try
% Images
for i = 1:numBlocks
    if mod(i,2) == 0
        currentCategory = 'kitchen';
        blockImages = kitchenImages;
        runTextures = kitchenTextures;
    else
        currentCategory = 'bathroom';
        blockImages = bathroomImages;
        runTextures = bathroomTextures;
    end

    numTrials = length(blockImages)*2;

    % randomize trial order
    % Ensure reproducible order across participants
    rng(i);  % Run number as seed

    % Initialize table
    blkImgs = table;
    blkImgs.texture = runTextures';
    blkImgs.category = repmat({currentCategory}, 1, length(runTextures))';
    blkImgs.imgName = {blockImages.name}';
    %blkImgs.iti = 

    % Shuffle rows
    blkImgs1 = blkImgs(randperm(height(blkImgs)),:);
    blkImgs2 = blkImgs(randperm(height(blkImgs)),:);
    
    blkImgs = [blkImgs1; blkImgs2];

%     % help file loggen
%     for iImg = 1:height(blkImgs)
%     fprintf(helpfileID, '%d\t%d\t%s\t%s\n', ...
%         i, iImg, char(blkImgs.category(iImg)), char(blkImgs.imgName(iImg)));
%     end

    % get targets
    load(fullfile(functionPath, 'targets.mat'), 'targetStruct')
    if testrun
        targetNum = numel(targetStruct);
    else
        targetNum = i; %str2double(str2double(i)); % = i
    end
    
    % init accuracy vector
    accuracy = nan(1, numel(targetStruct(targetNum).imgName));
    
    
    % Experiment Start
%         KbWait; % Wait for any key press in dummy mode

        trialOnsets = nan(1, numTrials);  % Store trial onset times
        itiOnsets = nan(1, numTrials);  % Store ITI onset times
        trialEnd = nan(1, numTrials);  % Store trial end
        responseTimes = nan(1, numTrials);  % Store response times
        responseKeys = cell(1, numTrials);  % Store response keys
        trialAccuracy = nan(1, numTrials);  % Store response accuracy
    
        % Display fixation cross before the trial
        DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Black fixation cross
        Screen('Flip', window);  % Show fixation cross
    
        % Wait for initial start pad (4s)
        WaitSecs(cfg.startPad);
    
        % Loop through all trials in the block
        for iImg = 1:height(blkImgs)
    
            % Initialize trial
            %cfg.iti = 0.75 * cfg.... + rand
            trialDuration = cfg.imageDuration + cfg.iti;
            responseFlag = false;
            itiFlag = false;
            responseTimes(iImg) = NaN;
            responseKeys{iImg} = 'none';
            trialAccuracy(iImg) = NaN;
            elapsedTime = 0;
            triggerTimeStamp = i;
            triggerDate = blkImgs.texture(iImg);
    
            % Present image
            Screen('DrawTexture', window, blkImgs.texture(iImg), [], image_rect);
            DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Black fixation cross
            trialOnsets(iImg) = Screen('Flip', window);  % Get trial onset time
    
            % trial timing
            while elapsedTime < trialDuration - ifi * 0.75
    
                [keyIsDown, ~, keyCode] = KbCheck;
    
                % Abort experiment when ESC is pressed
                if keyCode(abortKey)
                    error('Experiment aborted by user')
                end
    
                % Show fixation cross after 250 ms
                if ~itiFlag && elapsedTime > cfg.imageDuration - ifi * 0.5
                    % Inter-trial interval (ITI)
                    DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Show fixation cross during ITI
                    itiOnsets(iImg) = Screen('Flip', window);  % Show fixation cross
                    itiFlag = true;
                end
    
                % Update timer
                elapsedTime = GetSecs - trialOnsets(iImg);
            end
    
            % Log trial end time
            trialEnd(iImg) = GetSecs;
    
            % check if target trial
            if ismember(iImg, targetStruct(targetNum).trialNum)
    
                % get target trial
                targetIdx = find(iImg == targetStruct(targetNum).trialNum);
    
                % check if target image is correct
                if strcmp(targetStruct(targetNum).imgName{targetIdx}, ...
                        blkImgs.imgName{iImg})
                    disp('Correct target was selected')
                else
                    disp(['Selected target: ', targetStruct(targetNum).imgName{targetIdx}])
                    disp(['Current trial: ',  blkImgs.imgName{iImg}])
                    error('Target names do not match')
                end
                         
                % show target message
                targetMsg = 'Was the following image the same as before?';
                DrawFormattedText(window, targetMsg, 'center', 'center', [0 0 0]);
                qTime = Screen('Flip', window);
                WaitSecs(2);

                % show target
                targetImg = targetStruct(targetNum).showedImgName{targetIdx};
                rowIdx = find(strcmp(blkImgs.imgName, targetImg));
                targetImgTexture = blkImgs.texture(rowIdx(1));

                Screen('DrawTexture', window, targetImgTexture, [], image_rect);
                DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Black fixation cross
                qTime = Screen('Flip', window);
        
                % Wait for response
                responseFlag = false;
                while ~responseFlag
                    % Check for button press (target detection task)
                    [keyIsDown, responseTime, keyCode] = KbCheck;
    
                    if keyIsDown
                    % Abort experiment when ESC is pressed
                        if keyCode(abortKey)
                            error('Experiment aborted by user')
                        end
                    end
    
                    % Store first button press
                    if keyIsDown && ~responseFlag
                        responseTimes(iImg) = responseTime;
                        responseKeys{iImg} = num2str(find(keyCode));
                        
                        % show fixation cross after response is provided
                        DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Show fixation cross during ITI
                        Screen('Flip', window); 
                        responseFlag = true;
                    end
    
                end
    
                % get accuracy if response was recorded
                if strcmp(responseKeys{iImg}, 'none')
                    accuracy(targetIdx) = NaN;
                else
                    if targetStruct(targetNum).targetPresent(targetIdx) == 1
                        accuracy(targetIdx) = str2double(responseKeys{iImg}) == presentKey;
                    else
                        accuracy(targetIdx) = str2double(responseKeys{iImg}) == absentKey;
                    end
                    trialAccuracy(iImg) = accuracy(targetIdx);
                end
            end
    
                % Display fixation cross before the next trial
                DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Black fixation cross
                Screen('Flip', window);  % Show fixation cross
                WaitSecs(0.75 * cfg.debugTimingFactor); % hier stand 2
        end
    
            fprintf(fileID, '%s\t%s\t%d\t%d\t%s\t%s\t%.6f\t%.6f\t%.6f\t%s\t%.6f\t%.6f\t%s\t%.6f\t%.6f\n', ...
                cfg.subjectID, i, iImg,...
                blkImgs.texture(iImg), char(blkImgs.category(iImg)), char(blkImgs.imgName(iImg)), ...
                trialOnsets(iImg), itiOnsets(iImg), trialEnd(iImg), ...
                char(triggerDate), triggerTimeStamp, trialOnsets(iImg) - triggerTimeStamp,...
                responseKeys{iImg}, responseTimes(iImg), trialAccuracy(iImg)); 
%         end
    
    % add accuracy
    overallAccuracy = [overallAccuracy accuracy];

    % Pause after block
        pauseText = sprintf('End of block %d (%s).\n You can take a short break.\nPress any key to continue.', i, currentCategory);
        DrawFormattedText(window, pauseText, 'center', 'center', [0 0 0]);
        Screen('Flip', window);
        KbWait;

end % end for-loop

% show feedback message
feedbackMsg = ['The end', newline,...
    'You answered ', num2str(round(mean(overallAccuracy, 'omitnan')*100)), ...
    '% of the questions correctly', newline, ...
    'Thank you!'];
DrawFormattedText(window, feedbackMsg, 'center', 'center', [0 0 0]);
Screen('Flip', window);

% Wait for end pad
WaitSecs(5 * cfg.debugTimingFactor);

% Close log file
fclose(fileID);
% fclose(helpfileID);

% Close Psychtoolbox Screen
Screen('CloseAll');
ShowCursor;

% Show experiment duration 
totalLength = GetSecs - triggerTimeStamp;
disp(['The run duration was: ', char(minutes(totalLength/60))]);

catch ME
        % Handle errors
        fprintf('An error occurred: %s\n', ME.message);
        Screen('CloseAll');  % Ensure the screen is closed in case of an error
        fclose(fileID);  % Close log file if open
%         fclose(helpfileID);  % Close log file if open
        ShowCursor;
end
