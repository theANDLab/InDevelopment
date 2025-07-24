function [] = takeBreak(window)
    
    % Let PTB know what keys to monitor.
    KbName('UnifyKeyNames');
    leftKey = KbName('z'); % left key 
    rightKey = KbName('m'); % right key
    repeatKey = KbName('r');
    escKey = KbName('ESCAPE'); %
    spaceKey = KbName('Space');
    activeKeys = [leftKey, rightKey, repeatKey, escKey, spaceKey];
    RestrictKeysForKbCheck(activeKeys);

    line1 = 'Take a';
    line2 = '\n QUICK BREAK';
    [nx, ny, ~] = DrawFormattedText(window, [line1 line2], 'center', 'center', 1);
    Screen('DrawText', window, '', nx, ny, [255, 0, 0, 255]);
    Screen('Flip', window);
    WaitSecs(.25);

    wait=true;
    while wait
        
        % Check if Escape or User keys are pressed
        [keyIsDown, ~, keyCode] = KbCheck;
        if(keyIsDown && keyCode(escKey)), break; end % If escape pressed exit loop
        
        [~, ~, keyCode] = KbCheck; 
        if keyCode(spaceKey)
            wait=false;

        end
        
    end % while wait
    
    
    %% Display fixation screen.
    
    % Load nemo image
    [nemoIm,~,alpha] = imread('media/nemo.png', 'png');
    nemoIm(:,:,4) = alpha;
    
    % Make 'texture' of nemo, to be later repeatedly drawn onto screen.
    textPtr = Screen('MakeTexture', window, nemoIm);
    
    % Draw 'texture' of nemo on screen.
    Screen('DrawTexture', window, textPtr, [],[],[], 0);
    Screen('Flip', window);
    WaitSecs(2);

end 