% Build the user's structure using the provided base paramters.
%   Structure:
%       .targetLocation      Index of the target, 0 if absent.
%       .distractorLocation  Array of clock-form of distractor locations (12 = top).
%       .distractorType      Array of visual-types for distractor at
%                                  associated location.
%       .groupDirection      'horz', 'vert', or 'none'. If set this should
%                                   override all of the target/distractors
%                                   in the trial.

function userStruct = GenUserStructureMixed(StructParams)

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
prevType = StructParams.NumDynamicDistractorTypes;

trialNum = 1;

trialsPerSet = 16;   % Ignore value in StructParams.TrialsPerSet. If lower
                     % than 16, we will truncate after shuffling later.

% For each SetSize...
for setSize = EP.SetSizes
    % For trial for that SetSize...
    for i = 1:trialsPerSet
        
        % True if adjacent distractor should be clockwise from target.
        % False if it should be counter-clockwise.
        adjDistractorIsCW = false;


        % True if target should be absent in current trial. False otherwise.
        targetIsVisible = false;
        
        numDistractorTypes = 1;
        
        groupDir = 'vert';
        
        if ismember(variator, [1:4, 9:12])
            groupDir = 'horz';
        end
        % Variator = {1..8}: homogenous     Variator = {9..16}: heterogenous.
        if variator > 8
            numDistractorTypes = 2;
            groupDir = 'none';
        end
        
        % Variator = {2,4,6,8}: Visible;     Variator = {1,3,5,7}: Absent.
        if ismember(variator, [1,2,5,6,9,10,13,14])
            adjDistractorIsCW = true;
        end
        
        if mod(variator, 2) == 0
            targetIsVisible = true;
        end
        
        userStruct(trialNum).groupDir = groupDir;
        
        
        
        
        % Randomly place target somewhere between locations 1 and 12.
        targetLocation = randi(12);

          
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
                adjacentDistLocation = mod(targetLocation,12)+1;
            % ... otherwise set ccw from target.
            else
                adjacentDistLocation = mod(targetLocation-2,12)+1;
            end
       
        
            % Randomly assign locations for the other distractors except for
            % at locations of target and adjacent distractor.
            numOtherDistractors = setSize - 2;
            p = randperm(12);
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
            nextType = 1 + mod(prevType, numDistractorTypes);
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
        variator = mod(variator, 16) + 1;
        trialNum = trialNum + 1;
    end
    
    if StructParams.TrialsPerSet < 16  
        % Copy the userStruct up until this set.
        userStructPre = userStruct(1: trialNum - 17);
        
        % For N=4, include additional restrictions to which trials we
        % select: 2 absent, 2 present, 2 homogenous, 2 hetereogenous.
        if StructParams.TrialsPerSet == 4
%             odds = 1:2:15;  % target absent
%             evens = 2:2:16;  % target present
            
            userStrucSub = userStruct(trialNum - 16:length(userStruct));
            tgtLoc = [userStrucSub.targetLocation]; % comma separated list expansion 
            
            tgtAbs = tgtLoc == 0;
            tgtAbsent = find(tgtAbs); 
            
            tgtPres = tgtLoc > 0;
            tgtPresent = find(tgtPres); 

            
            requirementsAreMet = false;
            
            while requirementsAreMet == false
                absIdxs = tgtAbsent(randperm(8,2)); % Randomly select two absent trials
                presIdxs = tgtPresent(randperm(8,2));    % Randomly select two present trials
                
                                
                if sum(absIdxs <= 8) == 1 && ... % if exactly 1 homogeneous absent
                    sum(absIdxs > 8) == 1 && ... % and exactly 1 heterogeneous absent
                    sum(presIdxs <= 8) == 1 && ... % and exactly 1 homogeneous present 
                    sum(presIdxs > 8) ==1 % and exactly 1 heterogeneous present
                
                    requirementsAreMet = true;
                    
                    idxsToKeep = [absIdxs, presIdxs];
                    
                end
%                 if sum(idxsToKeep <= 8) == 2 && ... % if exactly 2 homogenous
%                    sum(idxsToKeep > 8) == 2         % and exactly 2 heterogenous
%                     
%                     requirementsAreMet = true;
%                 end
            end
            % Shuffle the four trials
            idxsToKeep = idxsToKeep(randperm(4));
        else
            % If we don't need 4 trialsPerSet, the restrictions don't
            % apply. Shuffle the 16 and choose the first four.
            idxsToKeep = randperm(16, StructParams.TrialsPerSet);
        end
            
        % Copy only the selected trials from the current set.
        userStructPost = userStruct(idxsToKeep + trialNum - 17);
        
        % Form a new userStruct based on the two structures formed above.
        userStruct = [userStructPre, userStructPost];
        
        % Reset trialNum now that we have shrunk the structure.
        trialNum = numel(userStruct) + 1;
    end
end

numTrials = numel(userStruct);

% Initialize to higher-than-desired values so we can enter the loop.
maxSetSize = EP.MaxSetSizesInRow + 1;
maxPresent = EP.MaxAbsOrPresInRow + 1;
maxMixedType = StructParams.MaxMixedTypeInRow + 1;
maxGroupDir = StructParams.MaxGroupDirInRow + 1;

% Generate a random permutation to shuffle the trial order, and apply this
% shuffling. Repeat this process until no more than N absent or present
% trials appear consecutively, AND no more than M trials of the same set
% size appear consecutivly.
while maxSetSize > EP.MaxSetSizesInRow ||    ...
      maxPresent > EP.MaxAbsOrPresInRow ||   ...
      maxMixedType > StructParams.MaxMixedTypeInRow ||   ...
      maxGroupDir > StructParams.MaxGroupDirInRow
  
    p = randperm(numTrials);
    userStruct = userStruct(p);
    
    % Array where present(i)=1 if trial number i has target present.
    present = [userStruct.targetLocation] ~= 0;

    % Array where setSize(i) setSize of trial number i.
    setSize = zeros(1, numTrials);
    for i=1:numTrials
        setSize(i) = numel(userStruct(i).distractorLocation) + present(i);
    end
    
    % Array where groupDir(i) is a symbol representing group direction of trial number i.
    groupDir = zeros(1,numTrials);
    for i=1:numTrials
        switch(userStruct(i).groupDir)
            case 'vert'
                groupDir(i) = 1;
            case 'horz'
                groupDir(i) = 2;
            case 'none'
                groupDir(i) = 3;
        end
    end
    
    % Array where mixedType(i) is True if trial i is heterogenous, else False.
    mixedType = groupDir == 3;
    
    
    
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
    
    
    
    % Find number of times the same groupDir appears in a row.
    i = find(diff(groupDir));
    n = [i numel(groupDir)] - [0 i];
    c = arrayfun(@(X) X-1:-1:0, n , 'un',0);
    groupDirY = cat(2,c{:});
    
    % Maximum value from the array is the worst number of repeats in the struct.
    maxGroupDir = max(groupDirY) + 1; 
    
    
    
    % Find number of times the same setSize appears in a row.
    i = find(diff(mixedType));
    n = [i numel(mixedType)] - [0 i];
    c = arrayfun(@(X) X-1:-1:0, n , 'un',0);
    mixedTypeY = cat(2,c{:});
    
    % Maximum value from the array is the worst number of repeats in the struct.
    maxMixedType = max(mixedTypeY) + 1; 
end

end