function [] = calculateVariables(SubID, Directory, SearchType, FeatureType, data)
%data=dataset2table(dataset('file','9999_feature_luminance_output.txt'));
%try
    
condition=strcat(SearchType, '_', FeatureType);

% figure1=figure('visible','off');
% histogram(data.RTms(data.RTms>0,:),5)
% saveas(figure1,strcat(Directory, string(SubID), '_', condition, '_', 'RThistogram.png'))

%% get summary table
ds = groupsummary(data, {'SetSize', ...
    'MotionPresent', 'GroupMotion', ...
    'Accuracy','TargetPresent'}, ...
    'mean', 'RTms', 'IncludeEmptyGroups',true);

ds.Acc = (ds.GroupCount / 16)*100;
    
writetable(ds, strcat(Directory, filesep, ...
    SubID, '_', condition, '_', 'data_summary.txt'), 'Delimiter', '\t')

% find which Set Sizes were tested 
setsize = unique(ds.SetSize)';

% create table with condition labels to later concatenate to data tables 
search = cellstr(repmat(SearchType,[length(setsize),1]));
feature = cellstr(repmat(FeatureType,[length(setsize),1]));
conditionTable = cell2table([search feature]);
conditionTable.Properties.VariableNames = {'SearchType', 'Feature'};
    

% Implementing the log-linear approach detailed here
% From Stanislaw & Todorov (1999):
% A third approach, dubbed loglinear, involves adding 0.5 to both the number of
% HITS and the number of FALSE ALARMS and adding 1 to both the number of
% target present (signal) and absent (noise) trials, before calculating the hit and false-alarm
% rates. This seems to work reasonably well (Hautus, 1995). Advocates of the
% loglinear approach recommend using it regardless of whether or not
% extreme rates are obtained.

    switch SearchType

        case 'feature'           
            %% Extract data 
            for i = 1:length(setsize)

                %------------------------------------------
                %           overall Feature search
                %	(collapsed (mixed) across motion type)
                %------------------------------------------


                % Signal Detection Theory measures
                Feature.h(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.Accuracy==1 ...
                    & ds.TargetPresent==1,:)) + .5;

                Feature.m(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.Accuracy==0 ...
                    & ds.TargetPresent==1,:));

                Feature.c(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.Accuracy==1 ...
                    & ds.TargetPresent==0,:));

                Feature.f(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.Accuracy==0 ...
                    & ds.TargetPresent==0,:)) + .5;


                % signal & noise trials
                Feature.p(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.TargetPresent==1,:)) + 1;

                Feature.a(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.TargetPresent==0,:)) + 1;



                % Reaction Time

                % Hits
                if isempty(mean(ds.mean_RTms(ds.SetSize==setsize(i) ...
                        & ds.Accuracy==1 & ds.TargetPresent==1,:), 'omitnan'))

                    Feature.hRT(i) = -1000;

                else

                    Feature.hRT(i) = mean(ds.mean_RTms(ds.SetSize==setsize(i) ...
                        & ds.Accuracy==1 & ds.TargetPresent==1,:), 'omitnan');

                end

                if isnan(Feature.hRT(i)), Feature.hRT(i) = -1000; end


                % False Alarms
                if isempty(mean(ds.mean_RTms(ds.SetSize==setsize(i) ...
                        & ds.Accuracy==0 & ds.TargetPresent==0,:), 'omitnan'))

                    Feature.fRT(i) = -1000;

                else

                    Feature.fRT(i) = mean(ds.mean_RTms(ds.SetSize==setsize(i) ...
                        & ds.Accuracy==0 & ds.TargetPresent==0,:), 'omitnan');

                end

                if isnan(Feature.fRT(i)), Feature.fRT(i) = -1000; end



                %------------------------------------------
                %       subset Dynamic Feature Search
                % need to do this collapsing across heterogeneous and
                % homogeneous motion in order to compare motion to static
                % then can follow up within the motion condition to compare
                % heterogeneous and homogeneous
                %------------------------------------------
                Motion = ds(ds.SetSize ==setsize(i) ...
                    & strcmp(ds.MotionPresent, 'Motion'), :);

                    %%% Dynamic Feature Search, homogeneous motion 

                    % Signal Detection Theory Measures 
                    DynamicHom.h(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.Accuracy==1 ...
                        & Motion.TargetPresent==1 ...
                        & (strcmp(Motion.GroupMotion, 'horz') ...
                        | strcmp(Motion.GroupMotion, 'vert')),:)) + .5;

                    DynamicHom.m(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.Accuracy==0 ...
                        & Motion.TargetPresent==1 ...
                        & (strcmp(Motion.GroupMotion, 'horz') ...
                        | strcmp(Motion.GroupMotion, 'vert')),:));

                    DynamicHom.c(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.Accuracy==1 ...
                        & Motion.TargetPresent==0 ...
                        & (strcmp(Motion.GroupMotion, 'horz') ...
                        | strcmp(Motion.GroupMotion, 'vert')),:));

                    DynamicHom.f(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.Accuracy==0 ...
                        & Motion.TargetPresent==0 ...
                        & (strcmp(Motion.GroupMotion, 'horz') ...
                        | strcmp(Motion.GroupMotion, 'vert')),:)) + .5;


                    % signal & noise trials
                    DynamicHom.p(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.TargetPresent==1 ...
                        & (strcmp(Motion.GroupMotion, 'horz') ...
                        | strcmp(Motion.GroupMotion, 'vert')),:)) + 1;

                    DynamicHom.a(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.TargetPresent==0 ...
                        & (strcmp(Motion.GroupMotion, 'horz') ...
                        | strcmp(Motion.GroupMotion, 'vert')),:)) + 1;


                    % Reaction Time

                    % Hits
                    if isempty(mean(Motion.mean_RTms(Motion.SetSize==setsize(i) ...
                            & Motion.Accuracy==1 ...
                            & Motion.TargetPresent==1 ...
                            & (strcmp(Motion.GroupMotion, 'horz') ...
                            | strcmp(Motion.GroupMotion, 'vert')),:), 'omitnan'))

                        DynamicHom.hRT(i) = -1000;

                    else

                        DynamicHom.hRT(i) = mean(Motion.mean_RTms(Motion.SetSize==setsize(i) ...
                            & Motion.Accuracy==1 ...
                            & Motion.TargetPresent==1 ...
                            & (strcmp(Motion.GroupMotion, 'horz') ...
                            | strcmp(Motion.GroupMotion, 'vert')),:), 'omitnan');

                    end

                    if isnan(DynamicHom.hRT(i)), DynamicHom.hRT(i) = -1000; end


                    % False Alarms
                    if isempty(mean(Motion.mean_RTms(Motion.SetSize==setsize(i) ...
                            & Motion.Accuracy==0 ...
                            & Motion.TargetPresent==0 ...
                            & (strcmp(Motion.GroupMotion, 'horz') ...
                            | strcmp(Motion.GroupMotion, 'vert')),:), 'omitnan'))

                        DynamicHom.fRT(i) = -1000;

                    else

                        DynamicHom.fRT(i) = mean(Motion.mean_RTms(Motion.SetSize==setsize(i) ...
                            & Motion.Accuracy==0 ...
                            & Motion.TargetPresent==0 ...
                            & (strcmp(Motion.GroupMotion, 'horz') ...
                            | strcmp(Motion.GroupMotion, 'vert')),:), 'omitnan');

                    end

                    if isnan(DynamicHom.fRT(i)), DynamicHom.fRT(i) = -1000; end



                    %%% Dynamic Feature Search, heterogeneous motion

                    % Signal Detection Theory Measures
                    DynamicHet.h(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.Accuracy==1 ...
                        & Motion.TargetPresent==1 ...
                        & strcmp(Motion.GroupMotion, 'none'),:)) + .5;

                    DynamicHet.m(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.Accuracy==0 ...
                        & Motion.TargetPresent==1 ...
                        & strcmp(Motion.GroupMotion, 'none'),:));

                    DynamicHet.c(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.Accuracy==1 ...
                        & Motion.TargetPresent==0 ...
                        & strcmp(Motion.GroupMotion, 'none'),:));

                    DynamicHet.f(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.Accuracy==0 ...
                        & Motion.TargetPresent==0 ...
                        & strcmp(Motion.GroupMotion, 'none'),:)) + .5;


                    % signal & noise trials
                    DynamicHet.p(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.TargetPresent==1 ...
                        & strcmp(Motion.GroupMotion, 'none'),:)) + 1;

                    DynamicHet.a(i) = sum(Motion.GroupCount(Motion.SetSize==setsize(i) ...
                        & Motion.TargetPresent==0 ...
                        & strcmp(Motion.GroupMotion, 'none'),:)) + 1;


                    % Reaction Time

                    % Hits
                    if isempty(mean(Motion.mean_RTms(Motion.SetSize==setsize(i) ...
                            & Motion.Accuracy==1 ...
                            & Motion.TargetPresent==1 ...
                            & strcmp(Motion.GroupMotion, 'none'),:), 'omitnan'))

                        DynamicHet.hRT(i) = -1000;

                    else

                        DynamicHet.hRT(i) = mean(Motion.mean_RTms(Motion.SetSize==setsize(i) ...
                            & Motion.Accuracy==1 ...
                            & Motion.TargetPresent==1 ...
                            & strcmp(Motion.GroupMotion, 'none'),:), 'omitnan');

                    end

                    if isnan(DynamicHet.hRT(i)), DynamicHet.hRT(i) = -1000; end


                    % False Alarms
                    if isempty(mean(Motion.mean_RTms(Motion.SetSize==setsize(i) ...
                            & Motion.Accuracy==0 ...
                            & Motion.TargetPresent==0 ...
                            & strcmp(Motion.GroupMotion, 'none'),:), 'omitnan'))

                        DynamicHet.fRT(i) = -1000;

                    else

                        DynamicHet.fRT(i) = mean(Motion.mean_RTms(Motion.SetSize==setsize(i) ...
                            & Motion.Accuracy==0 ...
                            & Motion.TargetPresent==0 ...
                            & strcmp(Motion.GroupMotion, 'none'),:), 'omitnan');

                    end

                    if isnan(DynamicHet.fRT(i)), DynamicHet.fRT(i) = -1000; end
                    


                %------------------------------------------
                %       subset Static Feature Search
                %------------------------------------------


                NoMotion = ds(ds.SetSize ==setsize(i) ...
                    & strcmp(ds.MotionPresent, 'NoMotion') ...
                    & strcmp(ds.GroupMotion, 'none'),:);

                    % Signal Detection Theory Measures
                    Static.h(i) = sum(NoMotion.GroupCount(NoMotion.SetSize==setsize(i) ...
                        & NoMotion.Accuracy==1 ...
                        & NoMotion.TargetPresent==1,:)) + .5;

                    Static.m(i) = sum(NoMotion.GroupCount(NoMotion.SetSize==setsize(i) ...
                        & NoMotion.Accuracy==0 ...
                        & NoMotion.TargetPresent==1,:));

                    Static.c(i) = sum(NoMotion.GroupCount(NoMotion.SetSize==setsize(i) ...
                        & NoMotion.Accuracy==1 ...
                        & NoMotion.TargetPresent==0,:));

                    Static.f(i) = sum(NoMotion.GroupCount(NoMotion.SetSize==setsize(i) ...
                        & NoMotion.Accuracy==0 ...
                        & NoMotion.TargetPresent==0,:)) + .5;



                    Static.p(i) = sum(NoMotion.GroupCount(NoMotion.SetSize==setsize(i) ...
                        & NoMotion.TargetPresent==1,:)) + 1;

                    Static.a(i) = sum(NoMotion.GroupCount(NoMotion.SetSize==setsize(i) ...
                        & NoMotion.TargetPresent==0,:)) + 1;


                    % Reaction Time

                    % Hits
                    if isempty(mean(NoMotion.mean_RTms(NoMotion.SetSize==setsize(i) ...
                        & NoMotion.Accuracy==1 ...
                        & NoMotion.TargetPresent==1,:), 'omitnan'))

                        Static.hRT(i) = -1000;

                    else

                        Static.hRT(i) = mean(NoMotion.mean_RTms(NoMotion.SetSize==setsize(i) ...
                        & NoMotion.Accuracy==1 ...
                        & NoMotion.TargetPresent==1,:), 'omitnan');

                    end

                    if isnan(Static.hRT(i)), Static.hRT(i) = -1000; end


                    % False Alarms
                    if isempty(mean(NoMotion.mean_RTms(NoMotion.SetSize==setsize(i) ...
                        & NoMotion.Accuracy==0 ...
                        & NoMotion.TargetPresent==0,:), 'omitnan'))

                        Static.fRT(i) = -1000;

                    else

                        Static.fRT(i) = mean(Motion.mean_RTms(Motion.SetSize==setsize(i) ...
                            & Motion.Accuracy==0 ...
                            & Motion.TargetPresent==0 ...
                            & strcmp(Motion.GroupMotion, 'none'),:), 'omitnan');

                    end

                    if isnan(Static.fRT(i)), Static.fRT(i) = -1000; end

            end
            
            %% Organize data 
            
            %------------------------------------------
            %       overall Feature search
            %	(collapsed across motion type)
            %------------------------------------------
            
            Feature.hits.values = Feature.h;
            Feature.hits.varnames = strcat(condition, '_', ...
                'overall', '_', 'Hits', '_', string(setsize));
            
            Feature.misses.values = Feature.m;
            Feature.misses.varnames = strcat(condition, '_', ...
                'overall', '_', 'Misses', '_', string(setsize));
            
            Feature.correctrejections.values = Feature.c;
            Feature.correctrejections.varnames = strcat(condition, '_', ...
                'overall', '_', 'CorrectRejections', '_', string(setsize));
            
            Feature.falsealarms.values = Feature.f;
            Feature.falsealarms.varnames = strcat(condition, '_', ...
                'overall', '_', 'FalseAlarms', '_', string(setsize));
            
            Feature.present.values = Feature.p;
            Feature.present.varnames = strcat(condition, '_', ...
                'overall', '_', 'Present', '_', string(setsize));
            
            Feature.absent.values = Feature.a;
            Feature.absent.varnames = strcat(condition, '_', ...
                'overall', '_', 'Absent', '_', string(setsize));
            
            Feature.hitsRT.values = Feature.hRT;
            Feature.hitsRT.varnames = strcat(condition, '_', ...
                'overall', '_', 'HitsRT', '_', string(setsize));
            
            Feature.falsealarmsRT.values = Feature.fRT;
            Feature.falsealarmsRT.varnames = strcat(condition, '_', ...
                'overall', '_', 'FalseAlarmsRT', '_', string(setsize));

            
            % Rates
            Feature.HR.values = Feature.hits.values ./ Feature.present.values;
            Feature.HR.varnames = strcat(condition, '_', ...
                'overall', '_', 'HitRate', '_', string(setsize));

            Feature.MR.values = Feature.misses.values ./ Feature.present.values; 
            Feature.MR.varnames = strcat(condition, '_', ...
                'overall', '_', 'MissRate', '_', string(setsize));

            Feature.CR.values = Feature.correctrejections.values ./ Feature.absent.values;
            Feature.CR.varnames = strcat(condition, '_', ...
                'overall', '_', 'CorrectRejectionRate', '_', string(setsize));

            Feature.FR.values = Feature.falsealarms.values ./ Feature.absent.values;
            Feature.FR.varnames = strcat(condition, '_', ...
                'overall', '_', 'FalseAlarmRate', '_', string(setsize));
            
            % 
                    
            Feature.long = array2table([setsize', ...
                Feature.hits.values', Feature.misses.values', ...
                Feature.correctrejections.values', Feature.falsealarms.values', ...
                Feature.present.values', Feature.absent.values', ...
                Feature.HR.values', Feature.MR.values', ...
                Feature.CR.values', Feature.FR.values', ...
                Feature.hitsRT.values', Feature.falsealarmsRT.values']);

            Feature.long.Properties.VariableNames = {'SetSize', ...
                'Hits', 'Misses', 'CorrectRejects', 'FalseAlarms', ...
                'Present', 'Absent', ...
                'HitRate', 'MissRate', 'CorrectRejectRate', 'FalseAlarmRate', ...
                'HitsRT', 'FalseAlarmRT'};


            Feature.motion = cell2table(cellstr(repmat('mixed',[length(setsize),1])));

            Feature.motion.Properties.VariableNames = {'MotionType'};

            Feature.long = [conditionTable Feature.motion Feature.long];
            
            writetable(Feature.long, strcat(Directory, filesep, SubID, '_',...
            condition, '_', 'overall', '_', 'summary.txt'), 'Delimiter', '\t')


            % calculate RT slope and error
            Feature.x = Feature.long.SetSize(Feature.long.HitsRT>0,:);
            Feature.y = Feature.long.HitsRT(Feature.long.HitsRT>0,:);
            [Feature.coefs, Feature.errorStruc]= polyfit(Feature.x,Feature.y,1);
            [Feature.yFit, Feature.stdError] = polyval(Feature.coefs,Feature.x,Feature.errorStruc);

            Feature.slopeRT.values = Feature.coefs(1);
            Feature.slopeRT.varnames = strcat(condition, '_', 'overall', '_', 'slopeRT');
            Feature.interceptRT.values = Feature.coefs(2);
            Feature.interceptRT.varnames = strcat(condition, '_', 'overall', '_', 'interceptRT');
            Feature.errorRT.values = Feature.stdError';
            
            for j = 1:length(Feature.stdError')

                Feature.errorRT.varnames(j) = strcat(condition, '_', 'overall', '_', 'RTerror', '_', string(j));

            end
            
            Feature.SubID = SubID;
            Feature.SearchType = SearchType;
            Feature.FeatureType = FeatureType;
            
            save(strcat(Directory, filesep, SubID, '_',... 
            condition, '_', 'overall', '.mat'), 'Feature')

            
            Feature.wide = array2table([Feature.hits.values, Feature.misses.values, ...
                Feature.correctrejections.values, Feature.falsealarms.values, ...
                Feature.present.values, Feature.absent.values, ...
                Feature.HR.values, Feature.MR.values, ...
                Feature.CR.values, Feature.FR.values, ...
                Feature.hitsRT.values, Feature.falsealarmsRT.values, ...
                Feature.slopeRT.values, Feature.interceptRT.values, Feature.errorRT.values]);
            
            Feature.wide.Properties.VariableNames = [...
                Feature.hits.varnames, Feature.misses.varnames, ...
                Feature.correctrejections.varnames, Feature.falsealarms.varnames, ...
                Feature.present.varnames, Feature.absent.varnames, ...
                Feature.HR.varnames, Feature.MR.varnames, ...
                Feature.CR.varnames, Feature.FR.varnames, ...
                Feature.hitsRT.varnames, Feature.falsealarmsRT.varnames, ...
                Feature.slopeRT.varnames, Feature.interceptRT.varnames, Feature.errorRT.varnames ...
                ];

            writetable(Feature.wide, strcat(Directory, filesep, SubID, '_',...
                condition, '_', 'overall', '_', ...
                'variables.txt'), 'Delimiter', '\t')

  
%             figure2=figure('visible','off');
%             plot(Feature.x,Feature.y,'bo')
%             hold on
%             plot(Feature.x,Feature.yFit,'r-')
%             plot(Feature.x,Feature.yFit+2*Feature.stdError,'m--', Feature.x,Feature.yFit-2*Feature.stdError,'m--')
%             title('Linear Fit of Data with 95% Prediction Interval')
%             legend('Data','Linear Fit','95% Prediction Interval')
%             saveas(figure2,strcat(Directory, filesep, ...
%                 string(SubID), '_', condition, '_', ...
%                 'overall', 'RTslope.png'))

            
            %------------------------------------------
            %       subset Dynamic Feature Search
            %               (homogeneous)
            %------------------------------------------
            
            DynamicHom.hits.values = DynamicHom.h;
            DynamicHom.hits.varnames = strcat(condition, '_', ...
                'DynamicHom', '_', 'Hits', '_', string(setsize));
            
            DynamicHom.misses.values = DynamicHom.m;
            DynamicHom.misses.varnames = strcat(condition, '_', ...
                'DynamicHom', '_', 'Misses', '_', string(setsize));
            
            DynamicHom.correctrejections.values = DynamicHom.c;
            DynamicHom.correctrejections.varnames = strcat(condition, '_', ...
                'DynamicHom', '_', 'CorrectRejections', '_', string(setsize));
            
            DynamicHom.falsealarms.values = DynamicHom.f;
            DynamicHom.falsealarms.varnames = strcat(condition, '_', ...
                'DynamicHom', 'FalseAlarms', '_', string(setsize));
            
            DynamicHom.present.values = DynamicHom.p;
            DynamicHom.present.varnames = strcat(condition, '_', ...
                'DynamicHom', '_', 'Present', '_', string(setsize));
            
            DynamicHom.absent.values = DynamicHom.a;
            DynamicHom.absent.varnames = strcat(condition, '_', ...
                'DynamicHom', '_', 'Absent', '_', string(setsize));
            
            DynamicHom.hitsRT.values = DynamicHom.hRT;
            DynamicHom.hitsRT.varnames = strcat(condition, '_', ...
                'DynamicHom', '_', 'HitsRT', '_', string(setsize));
            
            DynamicHom.falsealarmsRT.values = DynamicHom.fRT;
            DynamicHom.falsealarmsRT.varnames = strcat(condition, '_', ...
                'DynamicHom', '_', 'FalseAlarmsRT', '_', string(setsize));

            
            % Rates
            DynamicHom.HR.values = DynamicHom.hits.values ./ DynamicHom.present.values;
            DynamicHom.HR.varnames = strcat(condition, '_', ...
                'DynamicHom', '_', 'HitRate', '_', string(setsize));

            DynamicHom.MR.values = DynamicHom.misses.values ./ DynamicHom.present.values; 
            DynamicHom.MR.varnames = strcat(condition, '_', ...
                'DynamicHom', '_', 'MissRate', '_', string(setsize));

            DynamicHom.CR.values = DynamicHom.correctrejections.values ./ DynamicHom.absent.values;
            DynamicHom.CR.varnames = strcat(condition, '_', ...
                'DynamicHom', '_', 'CorrectRejectionRate', '_', string(setsize));

            DynamicHom.FR.values = DynamicHom.falsealarms.values ./ DynamicHom.absent.values;
            DynamicHom.FR.varnames = strcat(condition, '_', ...
                'DynamicHom', '_', 'FalseAlarmRate', '_', string(setsize));
            
            % 
                    
            DynamicHom.long = array2table([setsize', ...
                DynamicHom.hits.values', DynamicHom.misses.values', ...
                DynamicHom.correctrejections.values', DynamicHom.falsealarms.values', ...
                DynamicHom.present.values', DynamicHom.absent.values', ...
                DynamicHom.HR.values', DynamicHom.MR.values', ...
                DynamicHom.CR.values', DynamicHom.FR.values', ...
                DynamicHom.hitsRT.values', DynamicHom.falsealarmsRT.values']);

            DynamicHom.long.Properties.VariableNames = {'SetSize', ...
                'Hits', 'Misses', 'CorrectRejects', 'FalseAlarms', ...
                'Present', 'Absent', ...
                'HitRate', 'MissRate', 'CorrectRejectRate', 'FalseAlarmRate', ...
                'HitsRT', 'FalseAlarmRT'};

            DynamicHom.motion = cell2table(cellstr(repmat('homogeneous',[length(setsize),1])));

            DynamicHom.motion.Properties.VariableNames = {'MotionType'};

            DynamicHom.long = [conditionTable DynamicHom.motion DynamicHom.long];
            
            writetable(DynamicHom.long, strcat(Directory, filesep, SubID, '_',...
            condition, '_', 'dynamic_homogeneous',...
            '_summary.txt'), 'Delimiter', '\t')


            % calculate RT slope and error
            DynamicHom.x = DynamicHom.long.SetSize(DynamicHom.long.HitsRT>0,:);
            DynamicHom.y = DynamicHom.long.HitsRT(DynamicHom.long.HitsRT>0,:);
            [DynamicHom.coefs, DynamicHom.errorStruc]= polyfit(DynamicHom.x,DynamicHom.y,1);
            [DynamicHom.yFit, DynamicHom.stdError] = polyval(DynamicHom.coefs,DynamicHom.x,DynamicHom.errorStruc);

            DynamicHom.slopeRT.values = DynamicHom.coefs(1);
            DynamicHom.slopeRT.varnames = strcat(condition, '_', 'dynamic_homogeneous', '_', 'slopeRT');
            DynamicHom.interceptRT.values = DynamicHom.coefs(2);
            DynamicHom.interceptRT.varnames = strcat(condition, '_', 'dynamic_homogeneous', '_', 'interceptRT');
            DynamicHom.errorRT.values = DynamicHom.stdError';
            
            for j = 1:length(DynamicHom.stdError')

                DynamicHom.errorRT.varnames(j) = strcat(condition, '_', 'dynamic_homogeneous', '_', 'RTerror', '_', string(j));

            end

            
            DynamicHom.SubID = SubID;
            DynamicHom.SearchType = SearchType;
            DynamicHom.FeatureType = FeatureType;
            
            save(strcat(Directory, filesep, SubID, '_',... 
            condition, '_', 'dynamic_homogeneous', '.mat'), 'DynamicHom')

            DynamicHom.wide = array2table([DynamicHom.hits.values, DynamicHom.misses.values, ...
                DynamicHom.correctrejections.values, DynamicHom.falsealarms.values, ...
                DynamicHom.present.values, DynamicHom.absent.values, ...
                DynamicHom.HR.values, DynamicHom.MR.values, ...
                DynamicHom.CR.values, DynamicHom.FR.values, ...
                DynamicHom.hitsRT.values, DynamicHom.falsealarmsRT.values, ...
                DynamicHom.slopeRT.values, DynamicHom.interceptRT.values, DynamicHom.errorRT.values]);

            
            DynamicHom.wide.Properties.VariableNames = [...
                DynamicHom.hits.varnames, DynamicHom.misses.varnames, ...
                DynamicHom.correctrejections.varnames, DynamicHom.falsealarms.varnames, ...
                DynamicHom.present.varnames, DynamicHom.absent.varnames, ...
                DynamicHom.HR.varnames, DynamicHom.MR.varnames, ...
                DynamicHom.CR.varnames, DynamicHom.FR.varnames, ...
                DynamicHom.hitsRT.varnames, DynamicHom.falsealarmsRT.varnames, ...
                DynamicHom.slopeRT.varnames, DynamicHom.interceptRT.varnames, DynamicHom.errorRT.varnames ...
                ];


            writetable(DynamicHom.wide, strcat(Directory, filesep, SubID, '_',...
                condition, '_', 'dynamic_homogeneous', '_', ...
                'variables.txt'), 'Delimiter', '\t')

%   
%             figure3=figure('visible','off');
%             plot(DynamicHom.x,DynamicHom.y,'bo')
%             hold on
%             plot(DynamicHom.x,DynamicHom.yFit,'r-')
%             plot(DynamicHom.x,DynamicHom.yFit+2*DynamicHom.stdError,'m--', DynamicHom.x,DynamicHom.yFit-2*DynamicHom.stdError,'m--')
%             title('Linear Fit of Data with 95% Prediction Interval')
%             legend('Data','Linear Fit','95% Prediction Interval')
%             saveas(figure3,strcat(Directory, filesep, ...
%                 string(SubID), '_', condition, '_', ...
%                 'dynamic_homogeneous', 'RTslope.png'))
            
            
            %------------------------------------------
            %       subset Dynamic Feature Search
            %               (heterogeneous)
            %------------------------------------------
            DynamicHet.hits.values = DynamicHet.h;
            DynamicHet.hits.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'Hits', '_', string(setsize));
            
            DynamicHet.misses.values = DynamicHet.m;
            DynamicHet.misses.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'Misses', '_', string(setsize));
            
            DynamicHet.correctrejections.values = DynamicHet.c;
            DynamicHet.correctrejections.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'CorrectRejections', '_', string(setsize));
            
            DynamicHet.falsealarms.values = DynamicHet.f;
            DynamicHet.falsealarms.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'FalseAlarms', '_', string(setsize));
            
            DynamicHet.present.values = DynamicHet.p;
            DynamicHet.present.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'Present', '_', string(setsize));
            
            DynamicHet.absent.values = DynamicHet.a;
            DynamicHet.absent.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'Absent', '_', string(setsize));
            
            DynamicHet.hitsRT.values = DynamicHet.hRT;
            DynamicHet.hitsRT.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'HitsRT', '_', string(setsize));
            
            DynamicHet.falsealarmsRT.values = DynamicHet.fRT;
            DynamicHet.falsealarmsRT.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'FalseAlarmsRT', '_', string(setsize));

            
            % Rates
            DynamicHet.HR.values = DynamicHet.hits.values ./ DynamicHet.present.values;
            DynamicHet.HR.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'HitRate', '_', string(setsize));

            DynamicHet.MR.values = DynamicHet.misses.values ./ DynamicHet.present.values; 
            DynamicHet.MR.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'MissRate', '_', string(setsize));

            DynamicHet.CR.values = DynamicHet.correctrejections.values ./ DynamicHet.absent.values;
            DynamicHet.CR.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'CorrectRejectionRate', '_', string(setsize));

            DynamicHet.FR.values = DynamicHet.falsealarms.values ./ DynamicHet.absent.values;
            DynamicHet.FR.varnames = strcat(condition, '_', ...
                'DynamicHet', '_', 'FalseAlarmRate', '_', string(setsize));
            
            % 
                    
            DynamicHet.long = array2table([setsize', ...
                DynamicHet.hits.values', DynamicHet.misses.values', ...
                DynamicHet.correctrejections.values', DynamicHet.falsealarms.values', ...
                DynamicHet.present.values', DynamicHet.absent.values', ...
                DynamicHet.HR.values', DynamicHet.MR.values', ...
                DynamicHet.CR.values', DynamicHet.FR.values', ...
                DynamicHet.hitsRT.values', DynamicHet.falsealarmsRT.values']);

            DynamicHet.long.Properties.VariableNames = {'SetSize', ...
                'Hits', 'Misses', 'CorrectRejects', 'FalseAlarms', ...
                'Present', 'Absent', ...
                'HitRate', 'MissRate', 'CorrectRejectRate', 'FalseAlarmRate', ...
                'HitsRT', 'FalseAlarmRT'};

            DynamicHet.motion = cell2table(cellstr(repmat('heterogeneous',[length(setsize),1])));

            DynamicHet.motion.Properties.VariableNames = {'MotionType'};

            DynamicHet.long = [conditionTable DynamicHet.motion DynamicHet.long];
            
            writetable(DynamicHet.long, strcat(Directory, filesep, SubID, '_',...
            condition, '_', 'dynamic_heterogeneous',...
            '_summary.txt'), 'Delimiter', '\t')


            % calculate RT slope and error
            DynamicHet.x = DynamicHet.long.SetSize(DynamicHet.long.HitsRT>0,:);
            DynamicHet.y = DynamicHet.long.HitsRT(DynamicHet.long.HitsRT>0,:);
            [DynamicHet.coefs, DynamicHet.errorStruc]= polyfit(DynamicHet.x,DynamicHet.y,1);
            [DynamicHet.yFit, DynamicHet.stdError] = polyval(DynamicHet.coefs,DynamicHet.x,DynamicHet.errorStruc);

            DynamicHet.slopeRT.values = DynamicHet.coefs(1);
            DynamicHet.slopeRT.varnames = strcat(condition, '_', 'dynamic_heterogeneous', '_', 'slopeRT');
            DynamicHet.interceptRT.values = DynamicHet.coefs(2);
            DynamicHet.interceptRT.varnames = strcat(condition, '_', 'dynamic_heterogeneous', '_', 'interceptRT');
            DynamicHet.errorRT.values = DynamicHet.stdError';
            
            for j = 1:length(DynamicHet.stdError')

                DynamicHet.errorRT.varnames(j) = strcat(condition, '_', 'dynamic_heterogeneous', '_', 'RTerror', '_', string(j));

            end

            DynamicHet.SubID = SubID;
            DynamicHet.SearchType = SearchType;
            DynamicHet.FeatureType = FeatureType;
            
            save(strcat(Directory, filesep, SubID, '_',... 
            condition, '_', 'dynamic_heterogeneous', '.mat'), 'DynamicHet')

            DynamicHet.wide = array2table([DynamicHet.hits.values, DynamicHet.misses.values, ...
                DynamicHet.correctrejections.values, DynamicHet.falsealarms.values, ...
                DynamicHet.present.values, DynamicHet.absent.values, ...
                DynamicHet.HR.values, DynamicHet.MR.values, ...
                DynamicHet.CR.values, DynamicHet.FR.values, ...
                DynamicHet.hitsRT.values, DynamicHet.falsealarmsRT.values, ...
                DynamicHet.slopeRT.values, DynamicHet.interceptRT.values, DynamicHet.errorRT.values]);

            
            DynamicHet.wide.Properties.VariableNames = [...
                DynamicHet.hits.varnames, DynamicHet.misses.varnames, ...
                DynamicHet.correctrejections.varnames, DynamicHet.falsealarms.varnames, ...
                DynamicHet.present.varnames, DynamicHet.absent.varnames, ...
                DynamicHet.HR.varnames, DynamicHet.MR.varnames, ...
                DynamicHet.CR.varnames, DynamicHet.FR.varnames, ...
                DynamicHet.hitsRT.varnames, DynamicHet.falsealarmsRT.varnames, ...
                DynamicHet.slopeRT.varnames, DynamicHet.interceptRT.varnames, DynamicHet.errorRT.varnames ...
                ];


            writetable(DynamicHet.wide, strcat(Directory, filesep, SubID, '_',...
                condition, '_', 'dynamic_heterogeneous', '_', ...
                'variables.txt'), 'Delimiter', '\t')

  
%             figure4=figure('visible','off');
%             plot(DynamicHet.x,DynamicHet.y,'bo')
%             hold on
%             plot(DynamicHet.x,DynamicHet.yFit,'r-')
%             plot(DynamicHet.x,DynamicHet.yFit+2*DynamicHet.stdError,'m--', DynamicHet.x,DynamicHet.yFit-2*DynamicHet.stdError,'m--')
%             title('Linear Fit of Data with 95% Prediction Interval')
%             legend('Data','Linear Fit','95% Prediction Interval')
%             saveas(figure4,strcat(Directory, filesep, ...
%                 string(SubID), '_', condition, '_', ...
%                 'dynamic_heterogeneous', 'RTslope.png'))
            
            
            %------------------------------------------
            %       subset Static Feature Search
            %------------------------------------------
            Static.hits.values = Static.h;
            Static.hits.varnames = strcat(condition, '_', ...
                'Static', '_', 'Hits', '_', string(setsize));
            
            Static.misses.values = Static.m;
            Static.misses.varnames = strcat(condition, '_', ...
                'Static','_', 'Misses', '_', string(setsize));
            
            Static.correctrejections.values = Static.c;
            Static.correctrejections.varnames = strcat(condition, '_', ...
                'Static', '_', 'CorrectRejections', '_', string(setsize));
            
            Static.falsealarms.values = Static.f;
            Static.falsealarms.varnames = strcat(condition, '_', ...
                'Static', '_', 'FalseAlarms', '_', string(setsize));
            
            Static.present.values = Static.p;
            Static.present.varnames = strcat(condition, '_', ...
                'Static', '_', 'Present', '_', string(setsize));
            
            Static.absent.values = Static.a;
            Static.absent.varnames = strcat(condition, '_', ...
                'Static', '_', 'Absent', '_', string(setsize));
            
            Static.hitsRT.values = Static.hRT;
            Static.hitsRT.varnames = strcat(condition, '_', ...
                'Static', '_', 'HitsRT', '_', string(setsize));
            
            Static.falsealarmsRT.values = Static.fRT;
            Static.falsealarmsRT.varnames = strcat(condition, '_', ...
                'Static', '_', 'FalseAlarmsRT', '_', string(setsize));

            
            % Rates
            Static.HR.values = Static.hits.values ./ Static.present.values;
            Static.HR.varnames = strcat(condition, '_', ...
                'Static', '_', 'HitRate', '_', string(setsize));

            Static.MR.values = Static.misses.values ./ Static.present.values; 
            Static.MR.varnames = strcat(condition, '_', ...
                'Static', '_', 'MissRate', '_', string(setsize));

            Static.CR.values = Static.correctrejections.values ./ Static.absent.values;
            Static.CR.varnames = strcat(condition, '_', ...
                'Static', '_', 'CorrectRejectionRate', '_', string(setsize));

            Static.FR.values = Static.falsealarms.values ./ Static.absent.values;
            Static.FR.varnames = strcat(condition, '_', ...
                'Static', '_', 'FalseAlarmRate', '_', string(setsize));
            
            % 
                    
            Static.long = array2table([setsize', ...
                Static.hits.values', Static.misses.values', ...
                Static.correctrejections.values', Static.falsealarms.values', ...
                Static.present.values', Static.absent.values', ...
                Static.HR.values', Static.MR.values', ...
                Static.CR.values', Static.FR.values', ...
                Static.hitsRT.values', Static.falsealarmsRT.values']);

            Static.long.Properties.VariableNames = {'SetSize', ...
                'Hits', 'Misses', 'CorrectRejects', 'FalseAlarms', ...
                'Present', 'Absent', ...
                'HitRate', 'MissRate', 'CorrectRejectRate', 'FalseAlarmRate', ...
                'HitsRT', 'FalseAlarmRT'};

            Static.long = [conditionTable Static.long];
            
            writetable(Static.long, strcat(Directory, filesep, SubID, '_',...
            condition, '_', 'static',...
            '_summary.txt'), 'Delimiter', '\t')


            % calculate RT slope and error
            Static.x = Static.long.SetSize(Static.long.HitsRT>0,:);
            Static.y = Static.long.HitsRT(Static.long.HitsRT>0,:);
            [Static.coefs, Static.errorStruc]= polyfit(Static.x,Static.y,1);
            [Static.yFit, Static.stdError] = polyval(Static.coefs,Static.x,Static.errorStruc);

            Static.slopeRT.values = Static.coefs(1);
            Static.slopeRT.varnames = strcat(condition, '_', 'static', '_', 'slopeRT');
            Static.interceptRT.values = Static.coefs(2);
            Static.interceptRT.varnames = strcat(condition, '_', 'static', '_', 'interceptRT');
            Static.errorRT.values = Static.stdError';
            
            for j = 1:length(Static.stdError')

                Static.errorRT.varnames(j) = strcat(condition, '_', 'static', '_', 'RTerror', '_', string(j));

            end

            Static.SubID = SubID;
            Static.SearchType = SearchType;
            Static.FeatureType = FeatureType;
            
            save(strcat(Directory, filesep, SubID, '_',... 
            condition, '_', 'static', '.mat'), 'Static')

            Static.wide = array2table([Static.hits.values, Static.misses.values, ...
                Static.correctrejections.values, Static.falsealarms.values, ...
                Static.present.values, Static.absent.values, ...
                Static.HR.values, Static.MR.values, ...
                Static.CR.values, Static.FR.values, ...
                Static.hitsRT.values, Static.falsealarmsRT.values, ...
                Static.slopeRT.values, Static.interceptRT.values, Static.errorRT.values]);

            
            Static.wide.Properties.VariableNames = [...
                Static.hits.varnames, Static.misses.varnames, ...
                Static.correctrejections.varnames, Static.falsealarms.varnames, ...
                Static.present.varnames, Static.absent.varnames, ...
                Static.HR.varnames, Static.MR.varnames, ...
                Static.CR.varnames, Static.FR.varnames, ...
                Static.hitsRT.varnames, Static.falsealarmsRT.varnames, ...
                Static.slopeRT.varnames, Static.interceptRT.varnames, Static.errorRT.varnames ...
                ];


            writetable(Static.wide, strcat(Directory, filesep, SubID, '_',...
                condition, '_', 'static', '_', ...
                'variables.txt'), 'Delimiter', '\t')

  
%             figure5=figure('visible','off');
%             plot(Static.x,Static.y,'bo')
%             hold on
%             plot(Static.x,Static.yFit,'r-')
%             plot(Static.x,Static.yFit+2*Static.stdError,'m--', Static.x,Static.yFit-2*Static.stdError,'m--')
%             title('Linear Fit of Data with 95% Prediction Interval')
%             legend('Data','Linear Fit','95% Prediction Interval')
%             saveas(figure5,strcat(Directory, filesep, ...
%                 string(SubID), '_', condition, '_', ...
%                 'static', 'RTslope.png'))



        case 'conjunction'
            
            for i = 1:length(setsize)

                %------------------------------------------
                %       overall Conjunction search
                %------------------------------------------

                % Signal Detection Theory measures
                Conjunction.h(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.Accuracy==1 ...
                    & ds.TargetPresent==1,:)) + .5;

                Conjunction.m(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.Accuracy==0 ...
                    & ds.TargetPresent==1,:));

                Conjunction.c(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.Accuracy==1 ...
                    & ds.TargetPresent==0,:));

                Conjunction.f(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.Accuracy==0 ...
                    & ds.TargetPresent==0,:)) + .5;


                % signal & noise trials
                Conjunction.p(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.TargetPresent==1,:)) + 1;

                Conjunction.a(i) = sum(ds.GroupCount(ds.SetSize==setsize(i) ...
                    & ds.TargetPresent==0,:)) + 1;



                % Reaction Time

                % Hits
                if isempty(mean(ds.mean_RTms(ds.SetSize==setsize(i) ...
                        & ds.Accuracy==1 & ds.TargetPresent==1,:), 'omitnan'))

                    Conjunction.hRT(i) = -1000;

                else

                    Conjunction.hRT(i) = mean(ds.mean_RTms(ds.SetSize==setsize(i) ...
                        & ds.Accuracy==1 & ds.TargetPresent==1,:), 'omitnan');

                end

                if isnan(Conjunction.hRT(i)), Conjunction.hRT(i) = -1000; end


                % False Alarms
                if isempty(mean(ds.mean_RTms(ds.SetSize==setsize(i) ...
                        & ds.Accuracy==0 & ds.TargetPresent==0,:), 'omitnan'))

                    Conjunction.fRT(i) = -1000;

                else

                    Conjunction.fRT(i) = mean(ds.mean_RTms(ds.SetSize==setsize(i) ...
                        & ds.Accuracy==0 & ds.TargetPresent==0,:), 'omitnan');

                end

                if isnan(Conjunction.fRT(i)), Conjunction.fRT(i) = -1000; end
                

            end

            
            %% Organize data 
            
            %------------------------------------------
            %       overall Conjunction search
            %	(collapsed across motion type)
            %------------------------------------------
            
            Conjunction.hits.values = Conjunction.h;
            Conjunction.hits.varnames = strcat(condition, '_', ...
                'overall', '_', 'Hits', '_', string(setsize));
            
            Conjunction.misses.values = Conjunction.m;
            Conjunction.misses.varnames = strcat(condition, '_', ...
                'overall', '_', 'Misses', '_', string(setsize));
            
            Conjunction.correctrejections.values = Conjunction.c;
            Conjunction.correctrejections.varnames = strcat(condition, '_', ...
                'overall', '_', 'CorrectRejections', '_', string(setsize));
            
            Conjunction.falsealarms.values = Conjunction.f;
            Conjunction.falsealarms.varnames = strcat(condition, '_', ...
                'overall', '_', 'FalseAlarms', '_', string(setsize));
            
            Conjunction.present.values = Conjunction.p;
            Conjunction.present.varnames = strcat(condition, '_', ...
                'overall', '_', 'Present', '_', string(setsize));
            
            Conjunction.absent.values = Conjunction.a;
            Conjunction.absent.varnames = strcat(condition, '_', ...
                'overall', '_', 'Absent', '_', string(setsize));
            
            Conjunction.hitsRT.values = Conjunction.hRT;
            Conjunction.hitsRT.varnames = strcat(condition, '_', ...
                'overall', '_', 'HitsRT', '_', string(setsize));
            
            Conjunction.falsealarmsRT.values = Conjunction.fRT;
            Conjunction.falsealarmsRT.varnames = strcat(condition, '_', ...
                'overall', '_', 'FalseAlarmsRT', '_', string(setsize));

            
            % Rates
            Conjunction.HR.values = Conjunction.hits.values ./ Conjunction.present.values;
            Conjunction.HR.varnames = strcat(condition, '_', ...
                'overall', '_', 'HitRate', '_', string(setsize));

            Conjunction.MR.values = Conjunction.misses.values ./ Conjunction.present.values; 
            Conjunction.MR.varnames = strcat(condition, '_', ...
                'overall', '_', 'MissRate', '_', string(setsize));

            Conjunction.CR.values = Conjunction.correctrejections.values ./ Conjunction.absent.values;
            Conjunction.CR.varnames = strcat(condition, '_', ...
                'overall', '_', 'CorrectRejectionRate', '_', string(setsize));

            Conjunction.FR.values = Conjunction.falsealarms.values ./ Conjunction.absent.values;
            Conjunction.FR.varnames = strcat(condition, '_', ...
                'overall', '_', 'FalseAlarmRate', '_', string(setsize));
            
            % 
                    
            Conjunction.long = array2table([setsize', ...
                Conjunction.hits.values', Conjunction.misses.values', ...
                Conjunction.correctrejections.values', Conjunction.falsealarms.values', ...
                Conjunction.present.values', Conjunction.absent.values', ...
                Conjunction.HR.values', Conjunction.MR.values', ...
                Conjunction.CR.values', Conjunction.FR.values', ...
                Conjunction.hitsRT.values', Conjunction.falsealarmsRT.values']);

            Conjunction.long.Properties.VariableNames = {'SetSize', ...
                'Hits', 'Misses', 'CorrectRejects', 'FalseAlarms', ...
                'Present', 'Absent', ...
                'HitRate', 'MissRate', 'CorrectRejectRate', 'FalseAlarmRate', ...
                'HitsRT', 'FalseAlarmRT'};


            Conjunction.motion = cell2table(cellstr(repmat('mixed',[length(setsize),1])));

            Conjunction.motion.Properties.VariableNames = {'MotionType'};

            Conjunction.long = [conditionTable Conjunction.motion Conjunction.long];
            
            writetable(Conjunction.long, strcat(Directory, filesep, SubID, '_',...
            condition, '_', 'overall', '_', 'summary.txt'), 'Delimiter', '\t')


            % calculate RT slope and error
            Conjunction.x = Conjunction.long.SetSize(Conjunction.long.HitsRT>0,:);
            Conjunction.y = Conjunction.long.HitsRT(Conjunction.long.HitsRT>0,:);
            [Conjunction.coefs, Conjunction.errorStruc]= polyfit(Conjunction.x,Conjunction.y,1);
            [Conjunction.yFit, Conjunction.stdError] = polyval(Conjunction.coefs,Conjunction.x,Conjunction.errorStruc);

            Conjunction.slopeRT.values = Conjunction.coefs(1);
            Conjunction.slopeRT.varnames = strcat(condition, '_', 'overall', '_', 'slopeRT');
            Conjunction.interceptRT.values = Conjunction.coefs(2);
            Conjunction.interceptRT.varnames = strcat(condition, '_', 'overall', '_', 'interceptRT');
            Conjunction.errorRT.values = Conjunction.stdError';
            
            for j = 1:length(Conjunction.stdError')

                Conjunction.errorRT.varnames(j) = strcat(condition, '_', 'overall', '_', 'RTerror', '_', string(j));

            end
            
            Conjunction.SubID = SubID;
            Conjunction.SearchType = SearchType;
            Conjunction.ConjunctionType = FeatureType;
            
            save(strcat(Directory, filesep, SubID, '_',... 
            condition, '_', 'overall', '.mat'), 'Conjunction')

            
            Conjunction.wide = array2table([Conjunction.hits.values, Conjunction.misses.values, ...
                Conjunction.correctrejections.values, Conjunction.falsealarms.values, ...
                Conjunction.present.values, Conjunction.absent.values, ...
                Conjunction.HR.values, Conjunction.MR.values, ...
                Conjunction.CR.values, Conjunction.FR.values, ...
                Conjunction.hitsRT.values, Conjunction.falsealarmsRT.values, ...
                Conjunction.slopeRT.values, Conjunction.interceptRT.values, Conjunction.errorRT.values]);
            
            Conjunction.wide.Properties.VariableNames = [...
                Conjunction.hits.varnames, Conjunction.misses.varnames, ...
                Conjunction.correctrejections.varnames, Conjunction.falsealarms.varnames, ...
                Conjunction.present.varnames, Conjunction.absent.varnames, ...
                Conjunction.HR.varnames, Conjunction.MR.varnames, ...
                Conjunction.CR.varnames, Conjunction.FR.varnames, ...
                Conjunction.hitsRT.varnames, Conjunction.falsealarmsRT.varnames, ...
                Conjunction.slopeRT.varnames, Conjunction.interceptRT.varnames, Conjunction.errorRT.varnames ...
                ];

            writetable(Conjunction.wide, strcat(Directory, filesep, SubID, '_',...
                condition, '_', 'overall', '_', ...
                'variables.txt'), 'Delimiter', '\t')

  
%             figure6=figure('visible','off');
%             plot(Conjunction.x,Conjunction.y,'bo')
%             hold on
%             plot(Conjunction.x,Conjunction.yFit,'r-')
%             plot(Conjunction.x,Conjunction.yFit+2*Conjunction.stdError,'m--', Conjunction.x,Conjunction.yFit-2*Conjunction.stdError,'m--')
%             title('Linear Fit of Data with 95% Prediction Interval')
%             legend('Data','Linear Fit','95% Prediction Interval')
%             saveas(figure2,strcat(Directory, filesep, ...
%                 string(SubID), '_', condition, '_', ...
%                 'overall', 'RTslope.png'))


    end
    
% catch
%     
%     error('Error running calculateVariables. Moving on.')
%     
% end 
