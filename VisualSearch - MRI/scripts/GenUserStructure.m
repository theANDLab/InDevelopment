% Build the user's structure using the provided base paramters.
%   Structure:
%       .targetLocation      Index of the target, 0 if absent.
%       .distractorLocation  Array of clock-form of distractor locations (12 = top).
%       .distractorType      Array of visual-types for distractor at
%                                  associated location.

function userStruct = GenUserStructure(StructParams)

global EP

% Initialize an empty structure.
userStruct = struct([]);

% Increments from 1 to 4 and back around. Determines whether a target is
% visible and whether the adjacent distractor is clockwise of target. There
% are equal numbers of trials where each compound condition is true (for
% trialsPerSetSize being multiples of 4).
variator = 1;

% Stores that previous type that we assigned to a visual. The next type
% will be prevType+1 (or 1 if prevType is the largest type number).
prevType = StructParams.NumDistractorTypes;

trialNum = 1;

% For each SetSize...
for setSize = EP.SetSizes
    % For trial for that SetSize...
    for i = 1:StructParams.TrialsPerSet
        
        % True if adjacent distractor should be clockwise from target.
        % False if it should be counter-clockwise.
        adjDistractorIsCW = false;


        % True if target should be absent in current trial. False otherwise.
        targetIsVisible = false;
        
        % Variator = {1,2}: CW;     Variator = {3,4}: CCW.
        if variator <= 2
            adjDistractorIsCW = true;
        end
        % Variator = {2,4}: Visible;     Variator = {1,3}: Absent.
        if mod(variator, 2) == 0
            targetIsVisible = true;
        end
        
        
        
        % Randomly place target somewhere between locations 1 and 12.
        targetLocation = randi(4);

          
        % Dont include an adjacentDist if no distractors should be shown at
        % all!
        if setSize == 1
            if targetIsVisible
                distractorLocations = [];
            else
                distractorLocations = [targetLocation];
            end
        else
            % If adjacent distractor should be clockwise from target...
            if adjDistractorIsCW
                adjacentDistLocation = mod(targetLocation,4)+1;
            % ... otherwise set ccw from target.
            else
                adjacentDistLocation = mod(targetLocation-2,4)+1;
            end
       
        
            % Randomly assign locations for the other distractors except for
            % at locations of target and adjacent distractor.
            numOtherDistractors = setSize - 2;
            p = randperm(4);
            if size(adjacentDistLocation,1) == 0    % If adjacent dist is []
                distractorLocations = p(p ~= targetLocation);
            else
                distractorLocations = p(p ~= targetLocation & p ~= adjacentDistLocation);
            end
            distractorLocations = distractorLocations(1:numOtherDistractors);
            
            % Add adjacent distractor to list of distractors.
            distractorLocations = [adjacentDistLocation, distractorLocations];
            
             % If target is absent add it to the list of distractors.
            if ~targetIsVisible
                distractorLocations = [targetLocation, distractorLocations];
            end
        end

        
        % Create distractor types in a continuous sequence from previous
        % for each distractor we have a location for.
        distractorTypes = [];
        for d = 1: numel(distractorLocations)
            nextType = 1 + mod(prevType, StructParams.NumDistractorTypes);
            distractorTypes = [distractorTypes, nextType];
            prevType = nextType;
        end
        % Shift our (array of) types up by one because target is type 1.
        if numel(distractorTypes) ~= 0
            distractorTypes = distractorTypes + 1;
        end
        
        
        % Save the local variables to the structure.
        if targetIsVisible
            userStruct(trialNum).targetLocation = targetLocation;
        else
            userStruct(trialNum).targetLocation = 0;
        end
        userStruct(trialNum).distractorLocation = distractorLocations;
        userStruct(trialNum).distractorTypes = distractorTypes;
        
        
        % Increase counters by 1.
        variator = mod(variator, 4) + 1;
        trialNum = trialNum + 1;
    end
end

numTrials = numel(userStruct);

% Initialize to higher-than-desired values so we can enter the loop.
maxSetSize = EP.MaxSetSizesInRow + 1;
maxPresent = EP.MaxAbsOrPresInRow + 1;

% Generate a random permutation to shuffle the trial order, and apply this
% shuffling. Repeat this process until no more than N absent or present
% trials appear consecutively, AND no more than M trials of the same set
% size appear consecutivly.
while maxSetSize > EP.MaxSetSizesInRow || maxPresent > EP.MaxAbsOrPresInRow
    p = randperm(numTrials);
    userStruct = userStruct(p);
    
    % Array where present(i)=1 if trial number i has target present.
    present = [userStruct.targetLocation] ~= 0;

    % Array where setSize(i) setSize of trial number i.
    setSize = zeros(1, numTrials);
    for i=1:numTrials
        setSize(i) = numel(userStruct(i).distractorLocation) + present(i);
    end
    
    % Find number of times present/absent trials are in a row.
    i = find(diff(present));
    n = [i numel(present)] - [0 i];
    c = arrayfun(@(X) X-1:-1:0, n , 'un',0);
    presentY = cat(2,c{:});
    
    % Maximum value from the array is the worst number of repeats in the struct.
    maxPresent = max(presentY) + 1;
    
    
    % Find number of times the same setSize appears in a row.
    i = find(diff(setSize));
    n = [i numel(setSize)] - [0 i];
    c = arrayfun(@(X) X-1:-1:0, n , 'un',0);
    setSizeY = cat(2,c{:});
    
    % Maximum value from the array is the worst number of repeats in the struct.
    maxSetSize = max(setSizeY) + 1; 
end

end