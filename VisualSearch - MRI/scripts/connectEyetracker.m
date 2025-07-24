function [] = connectEyetracker(SubID)

%% set up connection with eyetracker

% Initialize EyeLink connection (dummymode = 0) or run in "Dummy Mode" without an EyeLink connection (dummymode = 1);
dummymode = 0;

% Optional: Set IP address of eyelink tracker computer to connect to.
% Call this before initializing an EyeLink connection if you want to use a non-default IP address for the Host PC.
%Eyelink('SetAddress', '10.10.10.240');

EyelinkInit(dummymode); % Initialize EyeLink connection
status = Eyelink('IsConnected');
if status < 1 % If EyeLink not connected
    dummymode = 1; 
end
   
% Name EyeLink Data file name entry. File name up to 8 characters
 edfFile = SubID; % Save file name to a variable 
% Print some text in Matlab's Command Window if a file name has not been entered
if  isempty(edfFile)
    fprintf('Session cancelled by user\n')
    error('Session cancelled by user'); % Abort experiment (see cleanup function below)
end     
% Print some text in Matlab's Command Window if file name is longer than 8 characters
if length(edfFile) > 8
    fprintf('Filename needs to be no more than 8 characters long (letters, numbers and underscores only)\n');
    error('Filename needs to be no more than 8 characters long (letters, numbers and underscores only)');
end

% Open an EDF file and name it
failOpen = Eyelink('OpenFile', edfFile);
if failOpen ~= 0 % Abort if it fails to open
    fprintf('Cannot create EDF file %s', edfFile); % Print some text in Matlab's Command Window
    error('Cannot create EDF file %s', edfFile); % Print some text in Matlab's Command Window
end

% Get EyeLink tracker and software version
% <ver> returns 0 if not connected
% <versionstring> returns 'EYELINK I', 'EYELINK II x.xx', 'EYELINK CL x.xx' where 'x.xx' is the software version
ELsoftwareVersion = 0; % Default EyeLink version in dummy mode
[ver, versionstring] = Eyelink('GetTrackerVersion');
if dummymode == 0 % If connected to EyeLink
    % Extract software version number. 
    [~, vnumcell] = regexp(versionstring,'.*?(\d)\.\d*?','Match','Tokens'); % Extract EL version before decimal point
    ELsoftwareVersion = str2double(vnumcell{1}{1}); % Returns 1 for EyeLink I, 2 for EyeLink II, 3/4 for EyeLink 1K, 5 for EyeLink 1KPlus, 6 for Portable Duo         
    % Print some text in Matlab's Command Window
    fprintf('Running experiment on %s version %d\n', versionstring, ver );
end
% Add a line of text in the EDF file to identify the current experimemt name and session. This is optional.
% If your text starts with "RECORDED BY " it will be available in DataViewer's Inspector window by clicking
% the EDF session node in the top panel and looking for the "Recorded By:" field in the bottom panel of the Inspector.
preambleText = sprintf('RECORDED BY Psychtoolbox demo %s session name: %s', mfilename, edfFile);
Eyelink('Command', 'add_file_preamble_text "%s"', preambleText);


%% STEP 2: SELECT AVAILABLE SAMPLE/EVENT DATA
% See EyeLinkProgrammers Guide manual > Useful EyeLink Commands > File Data Control & Link Data Control

% Select which events are saved in the EDF file. Include everything just in case
Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
% Select which events are available online for gaze-contingent experiments. Include everything just in case
Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,BUTTON,FIXUPDATE,INPUT');
% Select which sample data is saved in EDF file or available online. Include everything just in case
if ELsoftwareVersion > 3  % Check tracker version and include 'HTARGET' to save head target sticker data for supported eye trackers
    Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,RAW,AREA,HTARGET,GAZERES,BUTTON,STATUS,INPUT');
    Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT');
else
    Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,RAW,AREA,GAZERES,BUTTON,STATUS,INPUT');
    Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
end