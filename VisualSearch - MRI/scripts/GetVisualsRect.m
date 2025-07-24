% Produce an Nx4 matrix of the bounds on position for each of N visuals.
%   Bounds vector for each visual is [xMin, yMin, xMax, yMax], denoting the
%   corner edge points of each visual.
function visualsRect = GetVisualsRect(trialStruct, baseLocsX, baseLocsY,...
                                      visualTypes)

% Nx1 matrix containing the radius of each of N visuals we want to display.
visualsRadii = MakeArrForVisuals(trialStruct, visualTypes, 3);

% Returns true if a target should be shown in this trial.
targetIsPresent = trialStruct.targetLocation ~= 0;

% Determine pixel-location of target, if any should appear.
visualsRect = [];
if targetIsPresent
    targetX = baseLocsX(trialStruct.targetLocation);
    targetY = baseLocsY(trialStruct.targetLocation);
    visualsRect = [targetX-visualsRadii(1), targetY-visualsRadii(1),...
                  targetX+visualsRadii(1), targetY+visualsRadii(1)];
    
    % Remove target radius to match format with distractors
    visualsRadii = visualsRadii(2:end); 
end


% Determine pixel-location of distractors.
distractorsX = baseLocsX(trialStruct.distractorLocation);
distractorsY = baseLocsY(trialStruct.distractorLocation);

% Concatenate target and distractor rects into a matrix of ALL visuals rects.
for i = 1:numel(trialStruct.distractorLocation)
    distRect = [distractorsX(i)-visualsRadii(i),...
                distractorsY(i)-visualsRadii(i),...
                distractorsX(i)+visualsRadii(i),...
                distractorsY(i)+visualsRadii(i)];
    visualsRect = [visualsRect; distRect];
end

end


% A general-purpose function which extracts the esthetic characteristic
% stored in column `vtIdx' of VisualTypes for the appropriate type of
% visual (as listed in trialStruct.distractorTypes), including for the
% target.
function arr = MakeArrForVisuals(trialStruct, VisualTypes, vtIdx)

% Returns true if a target should be shown in this trial.
targetIsPresent = trialStruct.targetLocation ~= 0;

% Determine target elem, if a target should appear.
targetElem = [];
if targetIsPresent
    targetElem = VisualTypes{1}{vtIdx};
end

% Determine elem of distractors according to their types.
distractorsElem = [];
for distractorType = trialStruct.distractorTypes
    distractorsElem = [distractorsElem; VisualTypes{distractorType}{vtIdx}];
end

% Concatenate list to produce array of elems for ALL visuals.
arr = [targetElem; distractorsElem];

end