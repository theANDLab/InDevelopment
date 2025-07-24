function [window] = disconnectEyetracker(SubID,SearchType,Condition,window)
%% STEP 6: CLOSE EDF FILE. TRANSFER EDF COPY TO DISPLAY PC. CLOSE EYELINK CONNECTION. FINISH UP

% Put tracker in idle/offline mode before closing file. Eyelink('SetOfflineMode') is recommended.
% However if Eyelink('Command', 'set_idle_mode') is used, allow 50ms before closing the file as shown in the commented code:
% Eyelink('Command', 'set_idle_mode');% Put tracker in idle/offline mode
% WaitSecs(0.05); % Allow some time for transition    
Eyelink('SetOfflineMode'); % Put tracker in idle/offline mode
Eyelink('Command', 'clear_screen 0'); % Clear Host PC backdrop graphics at the end of the experiment
WaitSecs(0.5); % Allow some time before closing and transferring file
Eyelink('CloseFile'); % Close EDF file on Host PC
% Transfer a copy of the EDF file to Display PC
window = transferFile(SubID, window);

movefile([SubID '.edf'],['C:\Users\ANDLab\Documents\Experiments\VisualSearch - Emma Thesis - Eyetracking\task\data\' SubID '\' SubID '_' SearchType '_' Condition '.edf'])