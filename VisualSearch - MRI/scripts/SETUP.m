function [window, outputDirectory] = SETUP(SubID)

global EP; % Experiment parameters
global CP; % Conjunction search parameters
global CV; % Conjunction search visual parameters
[EP, CP, CV] = params();

%% Master settings

% This script calls Psychtoolbox commands available only in OpenGL-based 
% versions of the Psychtoolbox. The Psychtoolbox command AssertPsychOpenGL will issue
% an error message if someone tries to execute this script on a computer without
% an OpenGL Psychtoolbox
AssertOpenGL;

% Set the screen number to the external secondary monitor if there is one
% connected
screenNumber = max(Screen('Screens'));
% screenNumber = 1;

% Skip sync tests for demo purposes only
% Screen('Preference', 'SkipSyncTests', 2);
% If you are testing this code on a computer with a sub-par video card, it may fail psychtoolbox's
% screen tests. If this happens, you can set the PTB to skip these tests by uncommenting the above
% line of code. The experiment will run adequately on your computer, but you should not conduct
% actual experiments this way (timing issues).

rng('shuffle');
% This call will reset the random seed based on clock time. This is necessary because MATLAB resets
% its seed to the same thing every time it boots up; thus, each participant would see the same
% random order.

% Setup keyboard
% Let PTB know what keys to monitor.
KbName('UnifyKeyNames');
responseKey = KbName('z'); % left key 
repeatKey = KbName('r');
escKey = KbName('ESCAPE'); %
spaceKey = KbName('Space');
activeKeys = [responseKey, repeatKey, escKey, spaceKey];
RestrictKeysForKbCheck(activeKeys);


%% Create output directory for participant
% Data is saved in the output directory identified below. A folder for the participant ID will be
% automatically created if it does not exist.

% Put the current directory on the path and build the output directory
addpath(cd);
outputDirectory = [pwd filesep 'data' filesep ...
    num2str(SubID) filesep];    % eg ..../data/1001/

% Create the output directory if it doesn't exist.
if ~exist(outputDirectory, 'dir'), mkdir(outputDirectory), end

%% Open an on screen window

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'EnableNative10BitFramebuffer');

[window, ~] = PsychImaging('OpenWindow', screenNumber, EP.BackgroundColor, [0 0 EP.W_WIDTH EP.W_HEIGHT], 32, 2,...
[], [],  kPsychNeed32BPCFloat);

[oldRange, ~ , ~] = Screen('ColorRange',window,1,0,1);

% Setup the text type for the window
Screen('TextFont', window, 'Arial');
Screen('TextSize', window, 36);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

end