function DONE(window)


%   Keys pressed by the subject often show up in the Matlab command window as
%   well, cluttering that window with useless character junk. You can prevent
%   this from happening by disabling keyboard input to Matlab: Add a
%   ListenChar(2); command at the beginning of your script and a
%   ListenChar(0); to the end of your script to enable/disable transmission of
%   keypresses to Matlab. If your script should abort and your keyboard is
%   dead, press CTRL+C to reenable keyboard input -- It is the same as
%   ListenChar(0). See 'help ListenChar' for more info.
ListenChar(2);
HideCursor;

%% Compute and load global parameters from exptParams()
global EP; % experiment parameters


%% Let PTB know what keys to monitor.
KbName('UnifyKeyNames');
leftKey = KbName('z'); % left key 
rightKey = KbName('m'); % right key
repeatKey = KbName('r');
escKey = KbName('ESCAPE'); %
spaceKey = KbName('Space');
activeKeys = [leftKey, rightKey, repeatKey, escKey, spaceKey];
RestrictKeysForKbCheck(activeKeys);

%% Display Completion screen

finishName = [pwd filesep 'media' filesep 'ColorFoundDory.png'];
[finishImg, ~, ~] = imread(finishName);
imgPtr = Screen('MakeTexture', window, finishImg);

[dstRect] = scaleImgToScreen(finishImg, 1, EP.W_WIDTH/2, EP.W_HEIGHT/2);

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
    
% Clear screen
clearvars
Screen('CloseAll')
end