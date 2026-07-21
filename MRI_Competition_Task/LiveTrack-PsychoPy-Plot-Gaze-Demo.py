# LiveTrack-PsychoPy-Plot-Gaze-Demo.py
# version 1.0
# 05/2018 JT

import LiveTrack
from psychopy import visual, core

# initialise LiveTrack
LiveTrack.Init()

# Request calibrated data
LiveTrack.SetResultsTypeCalibrated()

# open a PsychoPy window
win = visual.Window(monitor='testMonitor', size=[1920, 1080], units="pix", fullscr=True, screen=1)

# run a loop to update the window with a dot at the position of the gaze
for x in range(0, 300):
    data = LiveTrack.GetLastResult()
    # plot a circle at the gaze position. Assume calibration in pixels 
    dot = visual.Circle(win,units='pix',radius=5,pos=[data.GazeX,data.GazeY],fillColor=[1,0,0],lineColor=[1,0,0]) #left eye is a red dot
    dot.draw()
    dot = visual.Circle(win,units='pix',radius=5,pos=[data.GazeXRight,data.GazeYRight],fillColor=[0,1,0],lineColor=[0,1,0]) #right eye is a green dot
    dot.draw()
    win.flip()

# close the PsychoPy window
win.close()

# close LiveTrack
LiveTrack.Close()
