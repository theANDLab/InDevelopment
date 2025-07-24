% Runs the visuals and interation with Psychtoolbox, playing the experiment
% with the characteristics of each trial in userStruct. expParams define
% some experiment parameters, while visualTypes the esthetics of the
% various visual types.
% Returns a numTrials-element array containing the user's response time for
% each trial, and -1 for no-response.
function [output] = RunTrial(window, SearchType, FeatureType, userStruct, structParams, visualTypes)

global EP;

% global iView
% global pSampleData


%%%%%%%%%%%%%%%%%%%%% Psychtoolbox Setup %%%%%%%%%%%%%%%%%%%%

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Let PTB know what keys to monitor.
KbName('UnifyKeyNames');
leftresponseKey = KbName('z'); % left key
rightresponseKey = KbName('m'); % right key
repeatKey = KbName('r');
escKey = KbName('ESCAPE'); %
spaceKey = KbName('Space');
activeKeys = [leftresponseKey, rightresponseKey, repeatKey, escKey, spaceKey];
RestrictKeysForKbCheck(activeKeys);

% Call GetSecs before trial-start to load Mex file for the first time. This
% avoids artificial time delays for later when timing really matters.
GetSecs();

% Set PTB verbosity to zero to improve performance. %TODO: Set to zero for
% production!
Screen('Preference', 'Verbosity', 3);



%%%%%%%%%%%%%%%%%%%%% Set up experiment start %%%%%%%%%%%%%%%%%%%%

% Get the centre coordinate of the window
screens = Screen('Screens');    % Get the screen numbers
screenNumber = max(screens);    % Draw to the external screen if avaliable
[xCenter, yCenter] = RectCenter(Screen('Rect', screenNumber));

% Calculate radius from center where dots should appear.
radius = EP.DotRadiusToCenter * Deg2Pix();

% Calculate the base position of all 4 possible circle locations.
baseLocationsX = radius * cos(2*pi* [1:4] ./4) + xCenter;
baseLocationsY = radius * sin(2*pi* [1:4] ./4) + yCenter;

% Rotate by three elements such that locations match "clock" orientation.
baseLocationsX = circshift(baseLocationsX, 2);
baseLocationsY = circshift(baseLocationsY, 2);

% Add nemo to circle center.
[nemoIm,~,alpha] = imread('media/nemo.png', 'png');
nemoIm(:,:,4) = alpha;
% Make 'texture' of nemo, to be later repeatedly drawn onto screen.
textPtr = Screen('MakeTexture', window, nemoIm);


% Load audio files
[nemoAudio, ~] = audioread('media/Nemo1.wav');
[doryAudio, Fs] = audioread('media/Dory1.wav');

% Upsample Nemo by factor of two to match Dory audio sampling period.
nemoAudio = repelem(nemoAudio,2,1);

% Add channel to Dory to move it from Mono to Stereo sound.
doryAudio = [doryAudio,doryAudio];

% Resample audio from 44.1kHz to 48kHz (for Windows audio driver).
nemoAudio = audioresample(nemoAudio, InputRate=Fs, OutputRate=48e3);
doryAudio = audioresample(doryAudio, InputRate=Fs, OutputRate=48e3);
Fs = 48e3;

InitializePsychSound(1); % Init sound driver and push for low-latency (1).
pahandle = PsychPortAudio('Open', [], 1, 1, Fs, 2);

numTrials = numel(userStruct);
userRT = -1 * ones(1,numTrials);
userResp = -1 * ones(1,numTrials);

% Maximum priority level
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);


%%%%%%%%%%%%%%%%%%% Run the experiment for each trial  %%%%%%%%%%%%%%%%%%%%%
keyIsDown = 1;
while keyIsDown
    [keyIsDown, ~, ~] = KbCheck;

    if ~keyIsDown
        keyIsDown = 0;

    end
end % while keyIsDown    

% set a switch to skip processing at end if esc is pressed 
processData = true; 

%begin running trials
for trialNum = 1:numTrials
    
    % Note: Per-trial preprocessing (this line until the while loop) takes
    % less than 10ms on average.
    
    % Nx4 matrix containing Rect pixel locations of each of N
    % visuals. The first visual is always the target, if one should be shown.                                   
    visualsRect = GetVisualsRect(userStruct(trialNum), baseLocationsX, ...
                                 baseLocationsY, visualTypes);
    
    % Nx2 matrix containing direction of each of N visuals. Eg [0,1]=vertical.
    motionRestrictor = MakeArrForVisuals(userStruct(trialNum), visualTypes, 4);
    motionRestrictor = [motionRestrictor, motionRestrictor];
    
    % Check if the groupDir field exists in the current userStructure.
    if isfield(userStruct(trialNum), 'groupDir')
        numVisuals = size(motionRestrictor,1);
        switch userStruct(trialNum).groupDir
            case 'horz'
                motionRestrictor = [ones(numVisuals,1), zeros(numVisuals,1)];
                motionRestrictor = [motionRestrictor, motionRestrictor];
            case 'vert'
                motionRestrictor = [zeros(numVisuals,1), ones(numVisuals,1)];
                motionRestrictor = [motionRestrictor, motionRestrictor];
        end
        
        motionType(trialNum) = {userStruct(trialNum).groupDir};
        
    elseif ~isfield(userStruct(trialNum), 'groupDir')
        
        motionType(trialNum) = {'none'};
        
    end
    
    visualsTexture = MakeArrForVisuals(userStruct(trialNum), visualTypes, 8);
    
    numVisuals = size(visualsRect,1);
    
    % Get maximum diameter that we will draw.
    %maxDotDiam = 2 * GetMaxDotRad(visualTypes);
    % no longer needed once Shape parameter is added
    
    % Dots will oscilate with a sine wave function in vert/horz motion. The
    % following parameters define the range of that motion and the frequency.
    motionAmpl = MakeArrForVisuals(userStruct(trialNum), visualTypes, 7);
    angFreq = 2 * pi * MakeArrForVisuals(userStruct(trialNum), visualTypes, 5);
    startPhase = MakeArrForVisuals(userStruct(trialNum), visualTypes, 6);
    time = 0;

    % check keyboard to make sure no key is down; don't move on until
    % released
    keyIsDown = 1;
    while keyIsDown
        [keyIsDown, ~, ~] = KbCheck;
        
        if ~keyIsDown
            keyIsDown = 0;

        end
    end % while keyIsDown  
    
    leftresponded = false;
    rightresponded = false;

    % Put tracker in idle/offline mode before recording. Eyelink('SetOfflineMode') is recommended 
    % however if Eyelink('Command', 'set_idle_mode') is used allow 50ms before recording as shown in the commented code:        
    % Eyelink('Command', 'set_idle_mode');% Put tracker in idle/offline mode before recording
    % WaitSecs(0.05); % Allow some time for transition        
    Eyelink('SetOfflineMode');% Put tracker in idle/offline mode before recording
    Eyelink('StartRecording'); % Start tracker recording
    WaitSecs(0.1); % Allow some time to record a few samples before presenting first stimulus
    
    % Get precise time at start of trial.
    preciseTime = GetSecs();     
    timeout = preciseTime + structParams.MaxTimePerTrial;
    loopIterCount = 0;
    
    % Loop the animation until timeout occurs or until key is released.
    while (GetSecs() < timeout) && ~leftresponded && ~rightresponded
        % Note: An iteration of this loop completes in approximately 17.5ms
        % (with a 60Hz refresh display: ~16.7ms per 'Flip').
        % Since keypress monitoring has been shown to take up almost none
        % of this time, the recorded response time can be delayed in the 
        % worst case by 17.5ms compared to the actual response.

        % Check that eye tracker is  still recording. Otherwise close and transfer copy of EDF file to Display PC
        err = Eyelink('CheckRecording');
        if(err ~= 0)
            fprintf('EyeLink Recording stopped!\n');
            % Transfer a copy of the EDF file to Display PC
            Eyelink('SetOfflineMode');% Put tracker in idle/offline mode
            Eyelink('CloseFile'); % Close EDF file on Host PC
            Eyelink('Command', 'clear_screen 0'); % Clear trial image on Host PC at the end of the experiment
            WaitSecs(0.1); % Allow some time for screen drawing
            % Transfer a copy of the EDF file to Display PC
            transferFile; % See transferFile function below)
            error('EyeLink is not in record mode when it should be. Unknown error. EDF transferred from Host PC to Display PC, please check its integrity.');
        end
        
        %%%%%%%%%%%%%%   Monitoring of keypresses   %%%%%%%%%%%%%%
        
        % Check if Escape or User keys are pressed
        [keyIsDown, ~, keyCode] = KbCheck;
        if(keyIsDown && keyCode(escKey))         
            % Write message to EDF file to mark the spacebar press time
            Eyelink('Message', 'EXP_TERMINATED') 
            break
        end % If escape pressed exit loop
        
        % userResp = -1, no response
        % userResp = 1, left response, z
        % userResp = 0, right reponse, m
        [keyIsDown, ~, keyCode] = KbCheck; 
        % If response key is down, store response time.
        if keyIsDown && keyCode(leftresponseKey)
            % Write message to EDF file to mark the spacebar press time
            Eyelink('Message', 'KEY_PRESSED');

            userRT(trialNum) = GetSecs() - preciseTime;
            userResp(trialNum) = 1;
            leftresponded = true;

        elseif keyIsDown && keyCode(rightresponseKey)
            % Write message to EDF file to mark the spacebar press time
            Eyelink('Message', 'KEY_PRESSED');
            
            userRT(trialNum) = GetSecs() - preciseTime;
            userResp(trialNum) = 0;
            rightresponded = true;
        end

        %%%%%%%%%%%%%%   Performing visual animation   %%%%%%%%%%%%%%
        
        % Calculate the OFFSET between each dot's base position and where
        % they must appear due to their oscillatory motion.
        motionDelta = motionAmpl .* repmat([sin(angFreq * time + startPhase), ...
                                            sin(angFreq * time + startPhase)], [1,2]) ...
                                 .* motionRestrictor;
   
        % Draw 'texture' of nemo on screen.
        Screen('DrawTexture', window, textPtr, [],[],[], 0);

        motionedRect = visualsRect + motionDelta;
        for i = 1:numVisuals
            Screen('DrawTexture', window, visualsTexture(i), [], motionedRect(i,:),[], 0);
        end

        % Send message here that trials started (include trial parameters) 
        % Flip to the newly-described screen.

        % Write message to EDF file to mark the start time of stimulus presentation.
        Eyelink('Message', 'STIM_ONSET');

        Screen('Flip', window);

        time = time + ifi;  % Increment the time for dot-movement.
        
        loopIterCount = loopIterCount+1;
        
    end

    % If escape was pressed BREAK to the end of the trials (ie stop
    % experiment).

    % Send message that ESC was pressed
    if(keyIsDown && keyCode(escKey)), processData = false; break; end

    % Play 'nemo' if correct, 'dory' if wrong.
    trialHasTarget = userStruct(trialNum).targetLocation ~= 0;
    if ~leftresponded == ~rightresponded
        PsychPortAudio('FillBuffer', pahandle, doryAudio');
    elseif trialHasTarget == leftresponded
        PsychPortAudio('FillBuffer', pahandle, nemoAudio');
    elseif trialHasTarget == rightresponded
        PsychPortAudio('FillBuffer', pahandle, doryAudio');
    elseif ~trialHasTarget == rightresponded
        PsychPortAudio('FillBuffer', pahandle, nemoAudio');
    elseif ~trialHasTarget == leftresponded
        PsychPortAudio('FillBuffer', pahandle, doryAudio');
    end

    PsychPortAudio('Start', pahandle);
    
    % Present nemo-only screen until the time since the trialStart exceed
    % (MaxTimePerTrial+MinTimePerWait) seconds.

    % Send message that key was pressed
    % Send message that fixation ITI is presented
    Screen('DrawTexture', window, textPtr, [],[],[], 0);

    % Write message to EDF file to mark the start time of stimulus presentation.
    Eyelink('Message', 'STIM_ONSET');

    Screen('Flip', window);

    waitTimeout = timeout + structParams.MinTimePerWait;
    WaitSecs('untilTime', waitTimeout);
    PsychPortAudio('Stop', pahandle);% Stop sound playback

    % Stop recording eye movements at the end of each trial
    WaitSecs(0.1); % Add 100 msec of data to catch final events before stopping
    Eyelink('StopRecording'); % Stop tracker recording

   % Write !V TRIAL_VAR messages to EDF file: creates trial variables in DataViewer
    % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Trial Message Commands
    Eyelink('Message', '!V TRIAL_VAR iteration %d', trialNum); % Trial iteration
    % Eyelink('Message', '!V TRIAL_VAR image %s', imgName); % Image name
    WaitSecs(0.001); % Allow some time between messages. Some messages can be lost if too many are written at the same time
    Eyelink('Message', '!V TRIAL_VAR rt %d', round(userRT(trialNum)*1000)); % Reaction time
    Eyelink('Message', '!V TRIAL_VAR acc %d', userResp(trialNum)); % User response (left or right button)
    WaitSecs(0.001); % Allow some time between messages. Some messages can be lost if too many are written at the same time
    % Write TRIAL_RESULT message to EDF file: marks the end of a trial for DataViewer
    % See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Defining the Start and End of a Trial
    Eyelink('Message', 'TRIAL_RESULT 0');
    WaitSecs(0.01); % Allow some time before ending the trial
end

if processData

targetPresentArr = [userStruct.targetLocation] ~= 0;
userAccuracyArr = targetPresentArr == userResp;

setSize = arrayfun(@(x)...
    length(userStruct(x).distractorLocation)+(userStruct(x).targetLocation > 0), ...
    1:numTrials);

userData = array2table(transpose([...
1:numTrials; ...
setSize; ...
targetPresentArr; ...
double(userAccuracyArr); ...
round(userRT*1000,0)]));

% make table of condition information 
search = cellstr(repmat(SearchType,numTrials,1)); 
feature = cellstr(repmat(FeatureType,numTrials,1));
groupMotion = motionType';
conditionTable = [search feature groupMotion];

output = [conditionTable userData];

output.Properties.VariableNames = {...
    'SearchType', 'Feature', 'GroupMotion', ...
    'Trial', 'SetSize', ...
    'TargetPresent', 'Accuracy', 'RTms'};

else
    
    Screen('CloseAll')
    error('trials ended early with ESC')

end


% Disable real-time mode.
Priority(0);

PsychPortAudio('Close', pahandle);% Close the audio device:
end


% Produce an Nx4 matrix of the bounds on position for each of N visuals.
%   Bounds vector for each visual is [xMin, yMin, xMax, yMax], denoting the
%   corner edge points of each visual.
function visualsRect = GetVisualsRect(trialStruct, baseLocsX, baseLocsY,...
                                      visualTypes)

% Nx1 matrix containing the radius of each of N visuals we want to display.
visualsRadii = MakeArrForVisuals(trialStruct, visualTypes, 3);

% Returns true if a target should be shown in this trial.
targetIsPresent = trialStruct.targetLocation ~= 0;

% Determine pixel-location of target, if any should appear.
visualsRect = [];
if targetIsPresent
    targetX = baseLocsX(trialStruct.targetLocation);
    targetY = baseLocsY(trialStruct.targetLocation);
    visualsRect = [targetX-visualsRadii(1), targetY-visualsRadii(1),...
                  targetX+visualsRadii(1), targetY+visualsRadii(1)];
    
    % Remove target radius to match format with distractors
    visualsRadii = visualsRadii(2:end); 
end


% Determine pixel-location of distractors.
distractorsX = baseLocsX(trialStruct.distractorLocation);
distractorsY = baseLocsY(trialStruct.distractorLocation);

% Concatenate target and distractor rects into a matrix of ALL visuals rects.
for i = 1:numel(trialStruct.distractorLocation)
    distRect = [distractorsX(i)-visualsRadii(i),...
                distractorsY(i)-visualsRadii(i),...
                distractorsX(i)+visualsRadii(i),...
                distractorsY(i)+visualsRadii(i)];
    visualsRect = [visualsRect; distRect];
end

end


% A general-purpose function which extracts the esthetic characteristic
% stored in column `vtIdx' of VisualTypes for the appropriate type of
% visual (as listed in trialStruct.distractorTypes), including for the
% target.
function arr = MakeArrForVisuals(trialStruct, VisualTypes, vtIdx)

% Returns true if a target should be shown in this trial.
targetIsPresent = trialStruct.targetLocation ~= 0;

% Determine target elem, if a target should appear.
targetElem = [];
if targetIsPresent
    targetElem = VisualTypes{1}{vtIdx};
end

% Determine elem of distractors according to their types.
distractorsElem = [];
for distractorType = trialStruct.distractorTypes
    distractorsElem = [distractorsElem; VisualTypes{distractorType}{vtIdx}];
end

% Concatenate list to produce array of elems for ALL visuals.
arr = [targetElem; distractorsElem];

end


% Returns the maximum dot radius from all the visual types. This is used to
% speed-up the circle drawing function.
function maxDotRad = GetMaxDotRad(visualTypes)

maxDotRad = visualTypes{1}{3};
for i = 2:numel(visualTypes)
    if visualTypes{i}{3} > maxDotRad
        maxDotRad = visualTypes{i}{3};
    end
end

end