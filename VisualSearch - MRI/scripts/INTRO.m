function  [window] = INTRO(window, searchType, featureType)

%   Keys pressed by the subject often show up in the Matlab command window as
%   well, cluttering that window with useless character junk. You can prevent
%   this from happening by disabling keyboard input to Matlab: Add a
%   ListenChar(2); command at the beginning of your script and a
%   ListenChar(0); to the end of your script to enable/disable transmission of
%   keypresses to Matlab. If your script should abort and your keyboard is
%   dead, press CTRL+C to reenable keyboard input -- It is the same as
%   ListenChar(0). See 'help ListenChar' for more info.
ListenChar(2);
%HideCursor;
% 
% %% Compute and load global parameters from exptParams()
global EP;
global CP;
global CV;

if strcmp(searchType, 'conjunction') && strcmp(featureType, 'luminance')
    
    % Translate the human-readable esthetics structure into something easier
    % for PsychToolBox to operate with.
    visualTypes = GetVisualTypes(CV.luminance, window);
    
    BaseStructParams = CP.luminance;
    
elseif strcmp(searchType, 'conjunction') && strcmp(featureType, 'color')
    
    visualTypes = GetVisualTypes(CV.color, window);
    
    BaseStructParams = CP.color;

elseif strcmp(searchType, 'conjunction') && strcmp(featureType, 'shape')
    
    visualTypes = GetVisualTypes(CV.shape, window);
    
    BaseStructParams = CP.shape;
    
else 
    error('ERROR: invalid searchType and/or featureType')
    
end

% make background color specific to feature manipulation 
Screen('FillRect', window, BaseStructParams.BackgroundColor);


%% Let PTB know what keys to monitor.
KbName('UnifyKeyNames');
responseKey = KbName('z'); % left key 
repeatKey = KbName('r');
escKey = KbName('ESCAPE'); %
spaceKey = KbName('Space');
activeKeys = [responseKey, repeatKey, escKey, spaceKey];
RestrictKeysForKbCheck(activeKeys);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% introduction
introFiles = dir(strcat(pwd, filesep,...
    'instructions', filesep,...
    'intro', filesep,...
    featureType, '*Instructions*png'));

for i = 1:size(introFiles,1)
        % specifiy and read instruction images
        introName = [introFiles(i).folder filesep introFiles(i).name];
        [introImg, ~, ~] = imread(introName, 'png');
        imgPtr = Screen('MakeTexture', window, introImg);

        [dstRect] = scaleImgToScreen(introImg, 1, EP.W_WIDTH/2, EP.W_HEIGHT/2);

        Screen('DrawTexture', window, imgPtr, [], dstRect);

        % display image
        Screen('Flip', window);
        WaitSecs(.25);

            % wait for spaceKey to advance to next image
            wait=true;
            while wait
                [~, ~, keyCode] = KbCheck; 

                if keyCode(spaceKey)
                    wait=false;

                end

            end % while wait
end % for introFiles

%% bring cursor back to command window
ListenChar(0);
ShowCursor;

end
