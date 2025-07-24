% Produces two Trial Structures (1 and 2) that can be passed through
% RunPractice() and RunTrial() respectively. The first is a tailored
% sequence of simple tasks, and the second is a subset of real trials the
% user could see.

% Note: Distractor types have been hardcoded (types 4 and 5). If you change
% the order of VisualTypesLegible, change these references too!
function [PracInstruct, PracTrials, PracMixedTrials] = GenPracticeStructure(...
    StructParams, numTrueTrials, SearchType)

pracStruct1 = struct([]);

% 1. Two trials should show a target and one distractor1 (e.g, red & horizontal
% movement), testing whether the subject can distinguish motion directions.
pracStruct1(1).targetLocation = 1;
pracStruct1(1).distractorLocation = 2;
pracStruct1(1).distractorTypes = 5;
pracStruct1(2).targetLocation = 3;
pracStruct1(2).distractorLocation = 2;
pracStruct1(2).distractorTypes = 5;

% 2. Two trials should show a target and one distractor2 (green & vertical 
% movement), testing whether the subject can distinguish colors.
pracStruct1(3).targetLocation = 2;
pracStruct1(3).distractorLocation = 1;
pracStruct1(3).distractorTypes = 4;
pracStruct1(4).targetLocation = 4;
pracStruct1(4).distractorLocation = 1;
pracStruct1(4).distractorTypes = 4;

% 3. Two trials without a target present.
pracStruct1(5).targetLocation = 0;
pracStruct1(5).distractorLocation = [1,3];
pracStruct1(5).distractorTypes = [4,4];
pracStruct1(6).targetLocation = 0;
pracStruct1(6).distractorLocation = [2,1];
pracStruct1(6).distractorTypes = [4,5];

% 4. 2 and 3 should be randomly ordered.
p =randperm(4);
q = 3:6;
q = [1:2, q(p)];
pracStruct1 = pracStruct1(q);

% 5. Following these 6 trials, X (numTrueTrials) trials should be randomly
% selected (determined the same way the general structure is) as ?true?
% practice trials.
pracStruct2 = GenUserStructure(StructParams);
pracStruct2 = pracStruct2(1:numTrueTrials);

if strcmp(SearchType, 'feature')

pracStruct3 = GenUserStructureMixed(StructParams);
pracStruct3 = pracStruct3(1:numTrueTrials);

elseif strcmp(SearchType, 'conjunction')
    
    pracStruct3 = pracStruct2(1:numTrueTrials);
    
else
    error('ERROR: invalid SearchType')
    
end

PracInstruct = pracStruct1;
PracTrials = pracStruct2;
PracMixedTrials = pracStruct3;

