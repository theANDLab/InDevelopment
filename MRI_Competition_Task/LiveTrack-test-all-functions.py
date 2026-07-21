#LiveTrack-test-all-functions.py
# version 1.1
import LiveTrack    # used for LiveTrack communication
import time         # used for time.sleep

# initialise LiveTrack
deviceType = LiveTrack.Init()
if deviceType==-100:
    raise Exception('LiveTrack: Error - LiveTrack not connected')

# get firmware version
fwversion = LiveTrack.GetFirmwareVersion()
print 'Firmware version is: ',str(fwversion)

# get library version
libversion = LiveTrack.GetLibraryVersion()
print 'Library version is: ',str(libversion)

# get serial number
serialno = LiveTrack.GetSerialNumber()
print 'Device serial number is: ',serialno

# get video capture configuration
width, height, rate, offsetX, offsetY = LiveTrack.GetCaptureConfig()
print 'Video image width: ',str(width)
print 'Video image height: ',str(height)
print 'Video capture rate: ',str(rate)
print 'video image offset X: ',str(offsetX)
print 'Video image offset Y: ',str(offsetY)

# Find out which eye to get fixation data from
[trackLeftEye, trackRightEye] = LiveTrack.GetTracking()

# Set which eye should be tracked
LiveTrack.SetTracking(trackLeftEye, trackRightEye)

# set pupil calibration
diameter = 5
pixels = 20
LiveTrack.SetPupilCalibration(diameter, pixels)

# get pupil calibration
pixelsToMM = LiveTrack.GetPupilCalibration()
print 'Pupil calibration is: ',str(pixelsToMM),' mm/pixels'

# set calibration
eye = 0 # left eye
cal = [1.1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1] # calibration matrix
viewDist = 500 # viewing distance in mm
xGlintMedian = 0 # median glint position during calibration (for head movement compensation)
yGlintMedian = 0 # median glint position during calibration (for head movement compensation)
LiveTrack.SetCalibration(eye, cal, viewDist, xGlintMedian, yGlintMedian)

# get calibration
cal, viewDist, xGlintMedian, yGlintMedian = LiveTrack.GetCalibration(eye)
print 'Calibration matrix is: ',str(cal)
print 'Viewing distance is set to: ',str(viewDist),' mm'
print 'yGlintMedian is set to: ',str(yGlintMedian)
print 'xGlintMedian is set to: ',str(xGlintMedian)

# Set result type to raw data
LiveTrack.SetResultsTypeRaw()

# Set result type to calibrated data
LiveTrack.SetResultsTypeCalibrated()

# calibrate with some fixation data
eye = 0 # left eye
numberOfFixationTargets = 9
targetsX = [-87.2753246410879,0,87.2753246410879,-87.2753246410879,0,87.2753246410879,-87.2753246410879,0,87.2753246410879] 
targetsY = [-87.2753246410879,-87.2753246410879,-87.2753246410879,0,0,0,87.2753246410879,87.2753246410879,87.2753246410879 ] 
vectX = [15.3326416015625,7.11785888671875,-0.772888183593750,14.8192138671875,6.92382812500000,-0.515563964843750,15.0279541015625,7.54492187500000,-0.210876464843750]
vectY = [-29.0014648437500,-28.1158447265625,-27.4469604492188,-20.2912597656250,-20.1955566406250,-19.9519042968750,-12.6470947265625,-12.0270996093750,-11.9754028320313 ]
viewDist = 500 # viewing distance in mm
LiveTrack.CalibrateDevice(eye, numberOfFixationTargets, targetsX, targetsY, vectX, vectY, viewDist, xGlintMedian, yGlintMedian)

# read back the new calibration
cal, viewDist, xGlintMedian, yGlintMedian = LiveTrack.GetCalibration(eye)
print 'Calibration matrix is: ',str(cal)
print 'Viewing distance is set to: ',str(viewDist),' mm'
print 'yGlintMedian is set to: ',str(yGlintMedian)
print 'xGlintMedian is set to: ',str(xGlintMedian)

# save the calibration to a file on the hard drive
LiveTrack.SaveCalibration("test.ltc")

# load the calibration from a file on the hard drive
LiveTrack.LoadCalibration("test.ltc")

# start buffuring data
LiveTrack.StartTracking()

# wait for one second
time.sleep(1)

# stop buffering data
LiveTrack.StopTracking()

# get the number of samples in the buffer and print it
count = LiveTrack.GetResultsCount()
print 'Number of samples in buffer is: '+str(count)

# get last data sample and print the time stamp
data = LiveTrack.GetLastResult()
print 'Timestamp of last data sample: '+str(data.Timestamp)

# get all the data in the buffer and print the first two time stamps
data = LiveTrack.GetBufferedEyePositions()
print 'Number of samples retreived from buffer: '+str(len(data))
print 'Timestamp of first sample: '+str(data[0].Timestamp)
print 'Timestamp of second sample: '+str(data[1].Timestamp)

# clear all data in the buffer
LiveTrack.ClearDataBuffer()

# get the number of samples in the buffer and print it
count = LiveTrack.GetResultsCount()
print 'Number of samples in buffer is: '+str(count)

# close LiveTrack
LiveTrack.Close()