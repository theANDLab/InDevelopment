%% Display fixation screen.

% Load nemo image
[nemoIm,~,alpha] = imread('media/nemo.png', 'png');
nemoIm(:,:,4) = alpha;

% Make 'texture' of nemo, to be later repeatedly drawn onto screen.
textPtr = Screen('MakeTexture', window, nemoIm);

% Draw 'texture' of nemo on screen.
Screen('DrawTexture', window, textPtr, [],[],[], 0);
Screen('Flip', window);
% Send message that fixation ITI is presented
WaitSecs(2);