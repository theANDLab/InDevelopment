import LiveTrackGS
import time

LiveTrackGS.VideoInit(0)
result= LiveTrackGS.VideoStart()

if result == 0:
    time.sleep(5)

    result= LiveTrackGS.VideoStartRecord("test.avi")
    time.sleep(5)

    LiveTrackGS.VideoStopRecord()
    time.sleep(5)

    LiveTrackGS.VideoStop()







