
%%%%%%%%%%%%% MAIN SCRIPT TO BE CALLED FOR EACH NEW USER %%%%%%%%%%%%%%%%%

% Clear the workspace and the screen, start from a fresh interpreter
sca;
clear all; 
close all;

addpath(strcat(pwd,filesep,'scripts'))

%% Experimenter input

% --------------------
%       SUBID
% --------------------
inputCheck = 0;
while inputCheck == 0
    SubID = input('Please enter SubID: ', 's');
    if ~isempty(SubID)
        inputCheck = 1;
    end
end

% --------------------------------------------------------------------------------
%                        COUNTERBALANCE TABLE
%	CB number	% Feature Order     
%	1           LUMINANCE, COLOR, SHAPE
%   2           LUMINANCE, SHAPE, COLOR
%   3           COLOR, LUMINANCE, SHAPE
%   4           COLOR, SHAPE, LUMINANCE
%   5           SHAPE, COLOR, LUMINANCE
%   6           SHAPE, LUMINANCE, COLOR
% --------------------------------------------------------------------------------
inputCheck = 0;
while inputCheck == 0
    CB = input('Please enter COUNTERBALANCE number: ', 's');
    if ~isempty(CB)
        inputCheck = 1;
    end
end

% Set up parameters and screen
[window, outputDirectory] = SETUP(SubID);


%% Run tasks
switch CB
    
    case '1'
    %	1           LUMINANCE, COLOR, SHAPE

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = INTRO(window,'conjunction','luminance');
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','luminance',outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','luminance', window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','color', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','color',window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'shape', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','shape',window);
        DONE(window);
               
    case '2'
    %	2           LUMINANCE, SHAPE, COLOR

            connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = INTRO(window,'conjunction','luminance');
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','luminance',outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','luminance', window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','shape', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','shape',window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'color', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','color',window);
        DONE(window);

        % window = INTRO(window, 'luminance');
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'luminance', outputDirectory);
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'shape', outputDirectory);
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'color', outputDirectory);
        % disconnectEyetracker % write this script
        % DONE(window);
        
        
    case '3'
    %	3           COLOR, LUMINANCE, SHAPE

            connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = INTRO(window,'conjunction','color');
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','color',outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','color', window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','luminance', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','luminance',window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'shape', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','shape',window);
        DONE(window);

        % window = INTRO(window, 'color');
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'color', outputDirectory);
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'luminance', outputDirectory);
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'shape', outputDirectory);
        % disconnectEyetracker % write this script
        % DONE(window);

    case '4'
    %	4           COLOR, SHAPE, LUMINANCE

                connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = INTRO(window,'conjunction','color');
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','color',outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','color', window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','shape', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','shape',window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'luminance', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','luminance',window);
        DONE(window);

        % window = INTRO(window, 'color');
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'color', outputDirectory);
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'shape', outputDirectory);
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'luminance', outputDirectory);
        % disconnectEyetracker % write this script
        % DONE(window);
        % 
    case '5'
    %	5           SHAPE, COLOR, LUMINANCE

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = INTRO(window,'conjunction','shape');
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','shape',outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','shape', window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','color', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','color',window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'luminance', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','luminance',window);
        DONE(window);

        % window = INTRO(window, 'shape');
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'shape', outputDirectory);
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'color', outputDirectory);
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'luminance', outputDirectory);
        % disconnectEyetracker % write this script
        % DONE(window);
        
    case '6'
    %	6           SHAPE, LUMINANCE, COLOR

                connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = INTRO(window,'conjunction','shape');
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','shape',outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','shape', window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction','luminance', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','luminance',window);

        connectEyetracker(SubID)
        [window] = calibrateEyetracker(window);
        window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'color', outputDirectory);
        window = disconnectEyetracker(SubID,'conjunction','color',window);
        DONE(window);

        % window = INTRO(window, 'shape');
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'shape', outputDirectory);
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'luminance', outputDirectory);
        % window = EXPERIMENT(window, SubID, CB, 1, 'conjunction', 'color', outputDirectory);
        % disconnectEyetracker % write this script
        % DONE(window);
        % 
    otherwise
        error(['ERROR: Invalid Counterbalance specification.', ...
            '\n CB number	% Feature Order', ...
            '\n 1           LUMINANCE, COLOR, SHAPE', ...
            '\n 2           LUMINANCE, SHAPE, COLOR', ...
            '\n 3           COLOR, LUMINANCE, SHAPE', ...
            '\n 4           COLOR, SHAPE, LUMINANCE', ...
            '\n 5           SHAPE, COLOR, LUMINANCE', ...
            '\n 6           SHAPE, LUMINANCE, COLOR'], ...
            class(length('ERROR: Invalid Counterbalance specification.')))

end


DONE(window);