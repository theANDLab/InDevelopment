clear left

left.t = []; % time stamp
left.X = []; % X coord of eye gaze
left.Y = []; % Y coord of eye gaze
left.P = []; % pupil diam (mm)
left.d = []; % eye distance from camera (mm)

right.t = []; % time stamp
right.X = []; % X coord of eye gaze
right.Y = []; % Y coord of eye gaze
right.P = []; % pupil diam (mm)
right.d = []; % eye distance from camera (mm)

    for i = 1:5
        response = 0;
        
        left(i).t = [];
        left(i).X = [];
        left(i).Y = [];
        left(i).P = [];
        left(i).d = [];
        
        right(i).t = [];
        right(i).X = [];
        right(i).Y = [];
        right(i).P = [];
        right(i).d = [];
        
        time = 0;
        tstart = GetSecs;
        
        while time < 3

            ret_sam = iView.iV_GetSample(pSampleData);
            
            if (ret_sam == 1)
                
                Smp = libstruct('SampleStruct', pSampleData);
                msg =	[int2str(i) '  ' int2str(Smp.timestamp) ' - GazeX: ' int2str(Smp.leftEye.gazeX) ' - GazeY: ' int2str(Smp.leftEye.gazeY)];
                disp(msg);
                
                left(i).t = [left(i).t, Smp.timestamp];
                left(i).X = [left(i).X, Smp.leftEye.gazeX];
                left(i).Y = [left(i).Y, Smp.leftEye.gazeY];
                left(i).P = [left(i).P, Smp.leftEye.diam];
                left(i).d = [left(i).d, Smp.leftEye.eyePositionZ];
                
                right(i).t = [right(i).t, Smp.timestamp];
                right(i).X = [right(i).X, Smp.rightEye.gazeX];
                right(i).Y = [right(i).Y, Smp.rightEye.gazeY];
                right(i).P = [right(i).P, Smp.rightEye.diam];
                right(i).d = [right(i).d, Smp.rightEye.eyePositionZ];
                
            else
                
                msg = 'Unable to get gaze samples';
                disp(msg)
                
            end
            
            pause(0.01);
            
            tend = GetSecs;
            time = tend - tstart;
            
        end
        
        for q = 2:length(left(i).t)
            left(i).tcorrected(q) = (left(i).t(q) - left(i).t(q-1))/1000;
        end
        
    end