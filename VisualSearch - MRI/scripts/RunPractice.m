% Runs the visuals and interaction with Psychtoolbox for the given practice
% structure pracStruct. Unlike RunTrial(), this version waits for keys
% between each set of visuals to proceed.
function RunPractice(window, pracStruct, visualTypes)

global EP;

%%%%%%%%%%%%%%%%%%%%% Psychtoolbox Setup %%%%%%%%%%%%%%%%%%%%

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Let PTB know what keys to monitor.
KbName('UnifyKeyNames');
responseKey = KbName('z'); % left key 
repeatKey = KbName('r');
escKey = KbName('ESCAPE'); %
spaceKey = KbName('Space');
activeKeys = [responseKey, repeatKey, escKey, spaceKey];
RestrictKeysForKbCheck(activeKeys);

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


numTrials = numel(pracStruct);


% Maximum priority level
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);


%%%%%%%%%%%%%%%%%%% Run the experiment for each trial  %%%%%%%%%%%%%%%%%%%%%


for trialNum = 1:numTrials
    % Note: Per-trial preprocessing (this line until the while loop) takes
    % less than 10ms on average.
    
    % Nx4 matrix containing Rect pixel locations of each of N
    % visuals. The first visual is always the target, if one should be shown.                                   
    visualsRect = GetVisualsRect(pracStruct(trialNum), baseLocationsX, ...
                                 baseLocationsY, visualTypes);
                                     
    % Nx3 matrix containing RGB (rows) color values for each of N visuals.
    visualsTexture = MakeArrForVisuals(pracStruct(trialNum), visualTypes, 8);

    % Nx2 matrix containing direction of each of N visuals. Eg [0,1]=vertical.
    motionRestrictor = MakeArrForVisuals(pracStruct(trialNum), visualTypes, 4);
    motionRestrictor = [motionRestrictor, motionRestrictor];
    
    
    numVisuals = size(visualsRect,1);
    
    % Get maximum diameter that we will draw.
    % maxDotDiam = 2 * GetMaxDotRad(visualTypes);
    
    
    % Dots will oscilate with a sine wave function in vert/horz motion. The
    % following parameters define the range of that motion and the frequency.
    motionAmpl = MakeArrForVisuals(pracStruct(trialNum), visualTypes, 7);
    angFreq = 2 * pi * MakeArrForVisuals(pracStruct(trialNum), visualTypes, 5);
    startPhase = MakeArrForVisuals(pracStruct(trialNum), visualTypes, 6);
    time = 0;

    % If the key is already down at the start of the trial, the first
    % release should be ignored. If the key is pressed before the visuals
    % come up, `keyReleaseIsValid' is marked FALSE. Until this flag is TRUE
    % (by releasing the key) a release will not lead to a stored result.
    [keyIsDown, ~, keyCode] = KbCheck;
    keyReleaseIsValid = true;
    if keyIsDown && keyCode(responseKey)     % If user key is down.
        keyReleaseIsValid = false;
    end
    
    % Stores whether the key has been pressed and released respectively,
    % in the current trial (not prior).
    keyWasPressed = false;
    keyReleased = false;

    
    % Loop the animation until timeout occurs or until key is released.
    while ~keyReleased
        
        %%%%%%%%%%%%%%   Monitoring of keypresses   %%%%%%%%%%%%%%
        
        % Check if Escape or User keys are pressed
        [keyIsDown, ~, keyCode] = KbCheck;
        if(keyIsDown && keyCode(escKey)), break; end % If escape pressed exit loop
        
        % Record event if user key is down (to later check for release). Only
        % mark keyWasPressed if user is not holding key from previous trial.
        if keyIsDown && keyCode(responseKey) && keyReleaseIsValid
            keyWasPressed = true;
        end
        
        % If user WAS holding key from previous trial and how now released
        % it, make any later key-presses (and releases) be valid.
        if ~keyCode(responseKey) && ~keyReleaseIsValid
            keyReleaseIsValid = true;
        end
        
        [keyIsDown, ~, keyCode] = KbCheck; 
        % If key was down at some point (in trial) but is no longer, store
        % response time.
        if ~keyCode(responseKey) && keyWasPressed
            keyReleased = true;
        end
        
        
        %%%%%%%%%%%%%%   Performing visual animation   %%%%%%%%%%%%%%
        
        % Calculate the OFFSET between each dot's base position and where
        % they must appear due to their oscillatory motion.
        motionDelta = motionAmpl .* repmat([sin(angFreq * time + startPhase), ...
                                            sin(angFreq * time + startPhase)], [1,2]) ...
                                 .* motionRestrictor;
   
        % Draw 'texture' of nemo on screen.
        Screen('DrawTexture', window, textPtr, [],[],[], 0);

        % motionedRect = visualsRect + motionDelta;
        % for i=1:numVisuals
        %     Screen('FillOval', window, visualsColor(i,:), motionedRect(i,:), maxDotDiam);
        % end

        motionedRect = visualsRect + motionDelta;
        for i = 1:numVisuals
            Screen('DrawTexture', window, visualsTexture(i), [], motionedRect(i,:),[], 0);
        end

        
    
        % Flip to the newly-described screen.
        Screen('Flip', window);

        time = time + ifi;  % Increment the time for dot-movement.
    end

    % If escape was pressed BREAK to the end of the trials (ie stop
    % experiment).
    if(keyIsDown && keyCode(escKey)), break; end
    
    
    
    
    %%%%%%%%%  Start waiting period of just nemo until keystroke  %%%%%%%%%
    if trialNum == numTrials
        Screen('TextFont', window, 'Helvetica');
        Screen('TextSize', window, 50);
        Screen('TextStyle', window, 1);
        
        line1 = 'Try some practice!';
        [nx, ny, ~] = DrawFormattedText(window, line1, 'center', 'center', 1);
        Screen('DrawText', window, '', nx, ny, [255, 0, 0, 255]);
        %Screen('DrawText', window, 'Ready for real practice?', xCenter-400, yCenter-100, [255, 255, 255]);
    else
        % Draw nemo and flip to screen
        Screen('DrawTexture', window, textPtr, [],[],[], 0);
    end
    Screen('Flip', window);
    KbWait;
end

% Disable real-time mode.
Priority(0);

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