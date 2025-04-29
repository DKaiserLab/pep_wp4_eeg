%% EEG Experiment

% Housekeeping
clear; close all

%% Collect subject ID
dat.subjectID = input('Enter subject number: ', 's');

% Additional metadata
dat.date = datestr(now, 'yyyy-mm-dd');
dat.time = datestr(now, 'HH:MM:SS');
dat.recordingModality = 'eye-tracking and eeg';
dat.viewing_dist_cm = 68;
dat.recordingLocation = 'math. dept. JLU Giessen';
dat.project = 'PEP_WP4_eeg';

%% General Experiment Configuration
testrun = false;
eegMode = false;
eyeTrackingMode = false;
debugTimingFactor = 1; % Must be 1 for accurat timing (< 1 will give faster timing)
imageDuration = 0.25 * debugTimingFactor;  % Image presentation time in seconds
startPad = 2 * debugTimingFactor;  % Time before the first trial in seconds
endPad = 2 * debugTimingFactor;  % Time after the last trial in seconds
ImageFileFormat = 'tif';
framAnticipation = 0.25;

numBlocks = 20;
categories = {'kitchen', 'bathroom'};

%% Paths
stimPath = fullfile(pwd,'..', 'stimuli');
outputPath = fullfile(pwd,'..', 'sourcedata', ['sub-', dat.subjectID], 'beh');
functionPath = fullfile(pwd,'utilities');

% add functions folder to path
addpath(functionPath)

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

%% Image Loading
kitchenImages = dir(fullfile(stimPath, 'kitchen', ['*.', ImageFileFormat]));
bathroomImages = dir(fullfile(stimPath, 'bathroom', ['*.', ImageFileFormat]));

if testrun
    kitchenImages = kitchenImages(1:20);
    bathroomImages = bathroomImages(1:20);
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
%
% Check if file name exists already to avoid overwriting
% if exist(runHelpFile, 'file')
%     error('File name exists already')
% end
% helpfileID = fopen(runHelpFile, 'w');
%
% % Write header for the log file
% fprintf(helpfileID, 'block\timgnumber\tcategory\timagename\n');

%% Output File Setup
% BIDS-compliant log file
runOutputFile = fullfile(outputPath, sprintf('sub-%s_task-eye-tracking-eeg_events.tsv', dat.subjectID));

% Check if file name exists already to avoid overwriting
if exist(runOutputFile, 'file')
    error('File name exists already')
end
fileID = fopen(runOutputFile, 'w');

% Write header for the log file
fprintf(fileID, 'subject\tblock\ttrial\ttexture\tcategory\timage\ttrialOnset\titiOnset\ttrialEnd\titi\tEEGtrigger\tresponseKey\tresponseTime\taccuracy\n');

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
[window, windowRect] = Screen('OpenWindow', screenNumber, [Color.gray]); % for Mac
[xCenter, yCenter] = RectCenter(windowRect);
ifi = Screen('GetFlipInterval', window);  % Get refresh interval
Screen('TextSize', window, 40) % define font size
DrawFormattedText(window, 'Loading...', 'center', 'center', [0 0 0]);
Screen('Flip', window);  % Show fixation cross
HideCursor; % Hide mouse cursor

%% Calculate size for desired degree of visual angle

% define visual angle
x_degree = 8;
y_degree = 6;

% Get the screen resolution in pixels per cm
[screen_hor, screen_vert] = Screen('DisplaySize', window); % width and height in mm
screen_hor = screen_hor / 10; % convert to cm
screen_vert = screen_vert / 10; % convert to cm

% Calculate pixels per centimeter
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
pixPerCmX = screenXpixels / screen_hor;
pixPerCmY = screenYpixels / screen_vert;

% Calculate the size in cm for the given visual angles
sizeCmX = 2 * dat.viewing_dist_cm * tan(deg2rad(x_degree) / 2);
sizeCmY = 2 * dat.viewing_dist_cm * tan(deg2rad(y_degree) / 2);

% Convert the size from cm to pixels
sizePixX = round(sizeCmX * pixPerCmX);
sizePixY = round(sizeCmY * pixPerCmY);

% get rectangle for image of correct size
image_rect = CenterRectOnPointd([0 0 sizePixX sizePixY], xCenter, yCenter);

% get smaller rectangle for targets (80% of original)
small_image_rect = CenterRectOnPointd([0 0 sizePixX*0.8 sizePixY*0.8], xCenter, yCenter);

%% Load all images as textures
for i = 1:length(bathroomImages)
    bathroomTextures(i) = Screen('MakeTexture', window, imread(fullfile(bathroomImages(i).folder, bathroomImages(i).name)));

    % add image information
    [~,file_name,ext] = fileparts(bathroomImages(i).name);
    stim_info(1,i).fInfo = bathroomImages(i);
    stim_info(1,i).fInfo.fname = file_name;
    stim_info(1,i).fInfo.ext = ext;
    stim_info(i).iInfo = imfinfo(fullfile(bathroomImages(i).folder, bathroomImages(i).name));
    stim_info(i).scrRect = image_rect;
end

for i = 1:length(kitchenImages)
    kitchenTextures(i) = Screen('MakeTexture', window, imread(fullfile(kitchenImages(i).folder, kitchenImages(i).name)));

    % add image information
    [~,file_name,ext] = fileparts(kitchenImages(i).name);
    stim_info(1,i + length(bathroomImages)).fInfo = kitchenImages(i);
    stim_info(1,i + length(bathroomImages)).fInfo.fname = file_name;
    stim_info(1,i + length(bathroomImages)).fInfo.ext = ext;
    stim_info(i + length(bathroomImages)).iInfo = imfinfo(fullfile(bathroomImages(i).folder, bathroomImages(i).name));
    stim_info(i + length(bathroomImages)).scrRect = image_rect;
end

%% Get setup struct and configure settings for eye-tracking
eyeTrackingMode
settings = Titta.getDefaults('Tobii Pro Fusion');
settings.debugMode = true; % Enable debug output
calViz = AnimatedCalibrationDisplay();
settings.cal.drawFunction = @calViz.doDraw;

% scale down the span of the calibration point
% (1.5 times as big as the presented simtuli)
scaling_factor = (sizePixX/screenXpixels) * 1.5;
centerPoint = 0.5;
top = 0.5 - 0.5 * scaling_factor;
bottom = 1 * scaling_factor + 0.5 - 0.5 * scaling_factor;

% get 9 calibration points (X-shaped)
calibrationPoints = [
    top, top;  % Top-left
    mean([top, centerPoint]), mean([top, centerPoint]);  % Intermediate top-left
    bottom, top;  % Top-right
    mean([bottom, centerPoint]), mean([top, centerPoint]);  % Intermediate top-right
    centerPoint, centerPoint;  % Center
    bottom, centerPoint;  % Middle-right
    mean([top, centerPoint]), mean([bottom, centerPoint]);  % Intermediate left-right
    top, bottom;  % Bottom-left
    mean([bottom, centerPoint]), mean([bottom, centerPoint]);  % Intermediate bottom-right
    bottom, bottom];  % Bottom-right

% get 5 validation points (diamant + center)
validationPoints = [
    centerPoint, centerPoint;  % Center
    mean([top, centerPoint]), centerPoint;  % Middle-left
    centerPoint, mean([top, centerPoint]);  % Middle-top
    mean([bottom, centerPoint]), centerPoint;  % Middle-right
    centerPoint, mean([bottom, centerPoint])];  % Middle-bottom

settings.cal.pointPos = calibrationPoints;
settings.val.pointPos = validationPoints;


%% Initialize Titta
EThndl = Titta(settings);
if ~eyeTrackingMode
    EThndl = EThndl.setDummyMode();
end
EThndl.init();

% add eye-tracker information to metadata
dat.recordingDevice = EThndl.deviceName;
dat.serialNumber = EThndl.serialNumber;
dat.samplingFrequency = EThndl.frequency;

%% Save participant meta data in a JSON file
datfilename = fullfile(outputPath, ['sub-', dat.subjectID, '_task-eye-tracking-eeg_participant.json']);
jsonText = jsonencode(dat);
fid = fopen(datfilename, 'w');
if fid == -1
    error('Cannot create JSON file');
end
fwrite(fid, jsonText, 'char');
fclose(fid);

%% Initialize keyboard
KbName('UnifyKeyNames');
abortKey = KbName('ESCAPE');

% define target keys
presentKey = 37;
absentKey = 39;

%% Initialize overall accuracy
overallAccuracy = [];

%% start experiment %%
try

    %% open serial port for EEG recording
    if eegMode
        SerialPortObj=serialport('COM4', 'TimeOut',1);%COM4 is the virtual serial output in the device manager
        SerialPortObj.BytesAvailableFcnMode='byte';
        SerialPortObj.BytesAvailableFcnCount=1;
        fopen(SerialPortObj);
        fwrite(SerialPortObj, 0, 'sync');
    end

    %% Initialize eye tracker calibration
    if eyeTrackingMode
        ListenChar(-1);
        tobii.calVal{1} = EThndl.calibrate(window);
        ListenChar(0);
    end

    %% Start recording eye-tracking data
    EThndl.buffer.start('gaze');
    WaitSecs(0.8);
    EThndl.sendMessage('start recording');

    %% start block-loop
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
        rng(i);  % block number as seed

        % Initialize table
        blkImgs = table;
        blkImgs.texture = runTextures';
        blkImgs.category = repmat({currentCategory}, 1, length(runTextures))';
        blkImgs.imgName = {blockImages.name}';

        % Shuffle rows
        blkImgs1 = blkImgs(randperm(height(blkImgs)),:);
        blkImgs2 = blkImgs(randperm(height(blkImgs)),:);
        blkImgs = [blkImgs1; blkImgs2];

        % get iti distribution
        itiDist = linspace(0.65, 0.85, height(blkImgs));
        itiDist = itiDist(randperm(length(itiDist)));
        itiDist = round(itiDist/ifi);
        itiDist = itiDist * ifi;
        blkImgs.iti = itiDist';

        % add triggercode (image ID = texture pointer - 10)
        blkImgs.EEGtrigger = blkImgs.texture - 10;

        %         % help file loggen
        %         for iImg = 1:height(blkImgs)
        %             fprintf(helpfileID, '%d\t%d\t%s\t%s\n', ...
        %                 i, iImg, char(blkImgs.category(iImg)), char(blkImgs.imgName(iImg)));
        %         end

        % get targets
        load(fullfile(functionPath, 'targets.mat'), 'targetStruct')
        if testrun
            targetNum = numel(targetStruct);
        else
            targetNum = i; %str2double(str2double(i));
        end

        % init accuracy vector
        accuracy = nan(1, numel(targetStruct(targetNum).imgName));

        % Experiment Start
        trialOnsets = nan(1, numTrials);  % Store trial onset times
        itiOnsets = nan(1, numTrials);  % Store ITI onset times
        trialEnd = nan(1, numTrials);  % Store trial end
        responseTimes = nan(1, numTrials);  % Store response times
        responseKeys = cell(1, numTrials);  % Store response keys
        trialAccuracy = nan(1, numTrials);  % Store response accuracy

        % Display fixation cross before the trial
        DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Black fixation cross
        Screen('Flip', window);  % Show fixation cross
        EThndl.sendMessage('FIX ON', GetSecs);

        % Wait for initial start pad
        WaitSecs(startPad);

        % Loop through all trials in the block
        for iImg = 1:height(blkImgs)

            % Initialize trial
            trialDuration = imageDuration + blkImgs.iti(i);
            responseFlag = false;
            itiFlag = false;
            responseTimes(iImg) = NaN;
            responseKeys{iImg} = 'none';
            trialAccuracy(iImg) = NaN;
            elapsedTime = 0;

            % Present image
            Screen('DrawTexture', window, blkImgs.texture(iImg), [], image_rect);
            DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Black fixation cross
            trialOnsets(iImg) = Screen('Flip', window);  % Get trial onset time

            % send EEG trigger
            if eegMode
                fwrite(SerialPortObj, blkImgs.EEGtrigger(iImg), 'sync');
                pause(0.01);
                fwrite(SerialPortObj, 0, 'sync');
            end

            % send message in eye-tracking file
            EThndl.sendMessage(sprintf('STIM ON: %d', blkImgs.EEGtrigger(iImg)), trialOnsets(iImg));

            % trial timing
            while elapsedTime < trialDuration - ifi * framAnticipation

                [keyIsDown, ~, keyCode] = KbCheck;

                % Abort experiment when ESC is pressed
                if keyCode(abortKey)
                    error('Experiment aborted by user')
                end

                % Show fixation cross after 250 ms
                if ~itiFlag && elapsedTime > imageDuration - ifi * framAnticipation
                    % Inter-trial interval (ITI)
                    DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);  % Show fixation cross during ITI
                    itiOnsets(iImg) = Screen('Flip', window);  % Show fixation cross
                    EThndl.sendMessage(sprintf('STIM OFF: %d', blkImgs.EEGtrigger(iImg)), itiOnsets(iImg));
                    EThndl.sendMessage('FIX ON', itiOnsets(iImg));
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
                targetMsg = 'Was this the same as the previous image?';
                DrawFormattedText(window, targetMsg, 'center', small_image_rect(2) - 50 , [0 0 0]);

                % show target
                targetImg = targetStruct(targetNum).targetImgName{targetIdx};
                rowIdx = find(strcmp(blkImgs.imgName, targetImg));
                targetImgTexture = blkImgs.texture(rowIdx(1));

                Screen('DrawTexture', window, targetImgTexture, [], small_image_rect);
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
            WaitSecs(0.75 * debugTimingFactor); % hier stand 2
            %end

            fprintf(fileID, ['%s\t%d\t%d\' ...
                't%d\t%s\t%s\' ...
                't%.6f\t%.6f\t%.6f\' ...
                't%d\t%d\' ...
                't%s\t%.6f\t%d\n'], ...
                dat.subjectID, i, iImg,...
                blkImgs.texture(iImg), char(blkImgs.category(iImg)), char(blkImgs.imgName(iImg)), ...
                trialOnsets(iImg), itiOnsets(iImg), trialEnd(iImg),...
                blkImgs.iti(iImg), blkImgs.EEGtrigger(iImg),...
                responseKeys{iImg}, responseTimes(iImg), trialAccuracy(iImg));
        end

        % send break message
        EThndl.sendMessage('START BREAK', GetSecs);

        % show feedback message
        feedbackMsg = sprintf('End of block %d.\n You answered %d %% of the questions correctly.', i, round(mean(accuracy, 'omitnan')*100));
        DrawFormattedText(window, feedbackMsg, 'center', 'center', [0 0 0]);
        Screen('Flip', window);
        WaitSecs(endPad)

        % Pause after block
        pauseText = 'You can take a short break.\n Press any key to continue.';
        DrawFormattedText(window, pauseText, 'center', 'center', [0 0 0]);
        Screen('Flip', window);
        KbWait;

        % send break end message
        EThndl.sendMessage('END BREAK', GetSecs);

        % add accuracy
        overallAccuracy = [overallAccuracy accuracy];

        %% Initialize eye tracker re-calibration after half the experiment
        if i == 10

            EThndl.sendMessage('RECALIBRATE', GetSecs);

            ListenChar(-1);
            tobii.calVal{1} = EThndl.calibrate(window);
            ListenChar(0);

            % Draw the fixation cross
            DrawFormattedText(window, '+', 'center', 'center', [0 0 0]);
            Screen('Flip', window);

            % start recording again
            EThndl.buffer.start('gaze');
            WaitSecs(0.8);
            EThndl.sendMessage('start recording');
        end

    end % end block-loop

    % show feedback message
    feedbackMsg = ['The end', newline,...
        'You answered ', num2str(round(mean(overallAccuracy, 'omitnan')*100)), ...
        '% of the questions correctly', newline, ...
        'Thank you!'];
    DrawFormattedText(window, feedbackMsg, 'center', 'center', [0 0 0]);
    Screen('Flip', window);

    % Wait for end pad
    WaitSecs(endPad);

    % Close log file
    fclose(fileID);
    %     fclose(helpfileID);

    % Close Psychtoolbox Screen
    Screen('CloseAll');
    ShowCursor;

    % send end message
    EThndl.sendMessage('END OF EXPERIMENT', GetSecs);

    % Show experiment duration
    totalLength = GetSecs - startTime;
    disp(['The run duration was: ', char(minutes(totalLength/60))]);

    % close serial port for EEG
    if eegMode
        fclose(SerialPortObj);
        delete(SerialPortObj);
        clear SerialPortObj;
    end

    %% Stop recording and save eye-tracking data
    EThndl.buffer.stop('gaze');
    WaitSecs(0.5)
    EThndl.sendMessage('STOP RECORDING', GetSecs);

    % save eye-tracking data
    ET_dat = EThndl.collectSessionData();
    ET_dat.expt.winRect = [0, 0, screenXpixels, screenYpixels]; % anaylsis scripts may need that information
    ET_dat.expt.resolution = [screenXpixels, screenYpixels];
    ET_dat.expt.stim = stim_info;
    EThndl.saveData(ET_dat, fullfile(outputPath, ['sub-', dat.subjectID, '_task-eye-tracking-eeg_physio']), true);

    % shut down
    EThndl.deInit();

catch ME

    % Handle errors
    fprintf('An error occurred: %s\n', ME.message);
    Screen('CloseAll');  % Ensure the screen is closed in case of an error
    fclose(fileID);  % Close log file if open
    %     fclose(helpfileID);  % Close log file if open
    ShowCursor;
    if eegMode
        fclose(SerialPortObj);
        delete(SerialPortObj);
        clear SerialPortObj;
    end

    % Stop and save recordingimgaussfilt
    EThndl.sendMessage('ERROR - PROGRAM ABORTED', GetSecs);
    EThndl.buffer.stop('gaze');
    WaitSecs(0.5)
    EThndl.sendMessage('STOP RECORDING', GetSecs);
    ET_dat = EThndl.collectSessionData();
    ET_dat.expt.winRect = [0, 0, screenXpixels, screenYpixels]; % anaylsis scripts may need that information
    ET_dat.expt.resolution = [screenXpixels, screenYpixels];
    ET_dat.expt.stim = stim_info;
    EThndl.saveData(ET_dat, fullfile(outputPath, ['sub-', dat.subjectID, '_task-eye-tracking-eeg_physio']), true);
    EThndl.deInit();

end
