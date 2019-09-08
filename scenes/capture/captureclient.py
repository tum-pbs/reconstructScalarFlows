#!/usr/bin/python
import os
import io
import sys
import picamera
import socket
import datetime
import math
import time
import picamera.array
import numpy as np
from time import sleep

# camera object
camera = None

# my client id
clientID = 0

# buffer containing the incomplete commands
recvBuffer = str()

# awbgains
#awb_rg, awb_bg = (0,0)

# sets the camera parameter
def setCamParams(shutterTime, awb, gain):
    print "st: ", shutterTime
    print "awb: ", awb
    print "gain: ", gain

    # is set to false if the gain adjustment times out
    gainOK = True

    # max deviation of gain from desired value
    maxDiff = 0.1

    currentGain = camera.analog_gain * camera.digital_gain
    if (gain * (1.0 + maxDiff) > currentGain or
       gain * (1.0 - maxDiff) < currentGain):
        # gain not correct -> adjust it

        camera.exposure_mode = 'off'
        camera.shutter_speed = 0
        camera.exposure_compensation = 0

        counter = 0
        phase = 'adjust'
        while (True):
            if phase == 'compare':
                # no check because this is used for adjustment
                currentGain = camera.analog_gain * camera.digital_gain
                print "gain adjustment: target=%s, current = %s, ec=%d"\
                    % (gain, currentGain, camera.exposure_compensation)

                if (gain * (1.0 - maxDiff) < currentGain and
                   gain * (1.0 + maxDiff) > currentGain):
                    # target reached
                    break
                # not yet reached -> adjust

                if currentGain == 0:
                    currentGain = 0.01  # avoid div zero

                if gain > currentGain:

                    # these are just experimental values, but they reduce the
                    # convergence time, compared to a linear correction
                    val = math.ceil(math.log(gain / currentGain, 2) * 3)
                    #val = math.ceil(gain / currentGain * 6)

                    # need to get gain higher
                    camera.exposure_compensation =\
                        min(camera.exposure_compensation + int(val), 25)
                else:
                    val = math.ceil(math.log(currentGain / gain, 2) * 3)
                    # val = math.ceil(currentGain / gain * 6)

                    # need to get gain lower
                    camera.exposure_compensation =\
                        max(camera.exposure_compensation - int(val), -25)

                camera.exposure_mode = 'auto'
                phase = 'adjust'
                sleep(0.9)

            else:
                camera.exposure_mode = 'off'
                phase = 'compare'
                sleep(0.1)

            counter += 1
            if(counter >= 100):
                gainOK = False
                print "Failed to reach desired gain value: target=%s; is=%s"\
                    % (gain, currentGain)
                break

    # camera.iso = max(min(int(gain * 100), 1600), 100)
    # gain is in range -> lock it
    camera.exposure_mode = 'off'
    # camera.exposure_mode = 'auto'
    camera.exposure_compensation = 0

    # other params can be set directly
    camera.awb_mode = 'off'
    camera.shutter_speed = shutterTime
    camera.awb_gains = awb
    sleep(0.3)
    print "awb is now: ", camera.awb_gains
    print "shutter_speed is now: ", camera.shutter_speed
    print "gain is now: ", (camera.analog_gain * camera.digital_gain)

    return gainOK


# gets the camera parameters
# format = (shutterTime, (awb), gain)
def getCamParams():
    return (camera.exposure_speed, camera.awb_gains,
            camera.analog_gain * camera.digital_gain)


# enables / disables automatic parameter control
def enableCamAutoMode(enable):
    if(enable):
        camera.exposure_mode = 'auto'
        camera.shutter_speed = 0
        camera.awb_mode = 'auto'

    else:
        camera.exposure_mode = 'off'
        camera.shutter_speed = camera.exposure_speed
        tmp = camera.awb_gains
        camera.awb_mode = 'off'
        camera.awb_gains = tmp

        # need to wait here a bit because the values are set delayed
        # (meaning that reading awb_gains would result in (0, 0))
        sleep(0.3)


# all my socket messages will follow the scheme: "<Control code>|<data>~"
def sendMsg(s, msg):
    s.sendall("%s~" % msg)


# waits until a full message is received
def getMsg(s):
    global recvBuffer

    while True:

        # receive until full message
        delim = recvBuffer.find("~")
        if(delim != -1):

            # full message -> extract it and remove from buffer
            result = recvBuffer[0:delim]
            recvBuffer = recvBuffer[delim + 1:]
            return result

        try:
            currentRecv = s.recv(4096, 0)

        except KeyboardInterrupt:
            print "Keyborad interrupt -> EXIT"
            camera.close()
            s.close()
            sys.exit(0)

        except:
            return ""

        if(len(currentRecv) == 0):
            # this means a empty string was received -> this should not happen
            return ''

        print "recv: %s" % currentRecv
        recvBuffer = recvBuffer + currentRecv


# checks if a number can be casted to a float
def isFloat(str):
    return str.replace(".", "", 1).isdigit()


def main(argv):

    global camera

    if len(argv) <= 1:
        # it requires one argument (the host ip)
        print "Missing arguments!\nUsage: flowcapture.py <control host>"
        return

    # read client id
    cidfile = open("/home/pi/flow.df/clientid.conf")
    clientID = int(cidfile.read())
    print "Client ID: %d" % clientID

    s = socket.socket()
    host = socket.gethostbyname(argv[1])

    try:
        # connect
        s.connect((host, 54321))

        # sen HI message with client id
        sendMsg(s, "HI|%d" % clientID)

        # wait for answer...
        m = getMsg(s)

        # ... and check if answer is expected
        if(m != ("CON|%d" % clientID)):
            print "Invalid answer from control host: %s" % m
            return

    except:
        print "Failed to connect to control host"
        return

    # get the camera object
    camera = picamera.PiCamera(resolution=(1920,1080),sensor_mode=1)

    # camera.start_preview()
    # camera.annotate_frame_num = True

    # currently recording?
    recording = False

    #currently streaming?
    streaming = False

    #socket handles for streaming
    streaming_socket = None
    connection = None

    # current manual parameters
    currentParams = (0, (0, 0), 0)

    # currentParams valid?
    paramsManual = False

    # main loop
    try:
        while True:

            # get a command
            msg = getMsg(s)

            # split command
            delim = msg.find("|")

            if (msg == "" or delim == -1):
                # command invalid
                print "Connection terminated or received invalid command"
                camera.close()
                s.close()
                sys.exit(0)

            # cmd  ~ command
            # data ~ data for command
            cmd = msg[0:delim]
            data = msg[delim + 1:]

            print "CMD: %s" % cmd

            if(cmd == "EXIT"):
                # end program
                s.close()
                sys.exit(0)

            elif(cmd == "CFP"):
                # Lock parameters

                print "Lock camera parameters"
                enableCamAutoMode(False)
                paramsManual = True
                currentParams = getCamParams()
                sendMsg(s, "OK|CFP")

            elif(cmd == "CGP"):
                # Get parameteres
                
                #old code for awb calculation
                #camera.awb_mode = 'off'
                #rg, bg = (0.5, 0.5)
                #camera.awb_gains = (rg, bg)
                #with picamera.array.PiRGBArray(camera, size=(128, 72)) as output:
                #        for i in range(50):
                #                camera.capture(output, format='rgb', resize=(128, 72), use_video_port=True)
                #                r, g, b = (np.mean(output.array[..., i]) for i in range(3))
                #                if abs(r - g) > 2:
                #                        if r > g:
                #                                rg -= 0.05
                #                        else:
                #                                rg += 0.05
                #                if abs(b - g) > 1:
                #                        if b > g:
                #                                bg -= 0.05
                #                        else:
                #                                bg += 0.05
                #                camera.awb_gains = (rg, bg)
                #                output.seek(0)
                #                output.truncate()
                ss, (awb_1, awb_2), g = getCamParams()
                # ss, (awb_1, awb_2), g = (1000, (0.5, 0.5), 1.5)
                print "Get camera params: shutter=%s, awb=(%s, %s), gain=%s"\
                    % (ss, awb_1, awb_2, g)

        
                # Host has default params, so this is no longer required

                #=========Setting parameters manually here (definitely not the best way to do this)========
                #superfluous code should be removed, but due to time constraints can't be done currently
                # 'g' is the parameter for exposure correction and '(awb_1, awb_2)' are the AWB parameters
                #(awb_1, awb_2) = (1.5,1.5)
                #g =6.0

                sendMsg(s, "CP|%d:%f:%f:%f" % (ss, awb_1, awb_2, g))

            elif(cmd == "CSP"):
                # Set parameters

                d = data.split(":")
                if(len(d) != 4 or (not isFloat(d[0])) or (not isFloat(d[1])) or
                   (not isFloat(d[2])) or (not isFloat(d[3]))):
                    # invalid data
                    print "Invalid CSP message: %s" % msg
                    sendMsg(s, "ERR|CSP")
                    continue

                print "Set cam params to: shutter=%s, awb=(%s, %s), gain=%s"\
                    % (d[0], d[1], d[2], d[3])

                gainOK = setCamParams(int(d[0]), (float(d[1]),
                                      float(d[2])), float(d[3]))

                enableCamAutoMode(False)
                paramsManual = True

                currentParams = (int(d[0]), (float(d[1]),
                                 float(d[2])), float(d[3]))

                # report either OK or ERR
                if gainOK:
                    sendMsg(s, "OK|CSP")
                else:
                    sendMsg(s, "ERR|CSP")

            elif(cmd == "RES"):
                # set resolution / framerate

                d = data.split(":")

                if(len(d) != 3 or (not d[0].isdigit()) or
                   (not d[1].isdigit()) or (not isFloat(d[2]))):
                    print "Invalid RES message: %s" % msg
                    sendMsg(s, "ERR|RES")
                    continue

		if(0):
		        print ("Set camera resolution to: w=%s, h=%s, fps=%s" % (d[0], d[1], d[2]))
		        camera.resolution = (int(d[0]), int(d[1]))
		        camera.framerate = float(d[2])
		else:
			a = 1920
			b = 1080
			print ("Set camera resolution to: w=%s, h=%s, fps=%s" % (a, b, d[2]))
		        camera.resolution = (int(a), int(b))
		        camera.framerate = float(d[2])

                sendMsg(s, "OK|RES")

            elif(cmd == "CAP"):
                # Automatic parameters

                print "Enable auto parameters"
                enableCamAutoMode(True)
                paramsManual = False
                sendMsg(s, "OK|CAP")

            elif(cmd == "SAM"):
                # Sample (single image)

                print "Taking sample picture: %s" % (data)

                sleep(0.5)
                camera.capture(data)

                sendMsg(s, "OK|SAM")

            elif(cmd == "REC"):
                # start recording
                
                if 0:
                    if recording:
                        print "Error: Still recording"
                        continue

                    no = int(sys.argv[2]);
                    print no

                    print "Start recording at: %s"\
                        % (datetime.datetime.now().time().isoformat())
                    recording = True
                    time.sleep(2)
                    # Set up 40 in-memory streams
                    outputs = [io.BytesIO() for i in range(no)]
                    start = time.time()
                    camera.capture_sequence(outputs, 'jpeg', use_video_port=True)
                    finish = time.time()
                    # How fast were we?
                    print('Captured '+str(no)+' images at %.2ffps' % (no / (finish - start)))
                    

                    for i in range(no):
                        with open('/home/pi/flow.df/streams/stream_%03d.jpg'%i, 'wb') as f:
                            f.write(outputs[i].getvalue())
                        f.close()
                
                
                
                # old code
                if 1:
                    if recording:
                        print "Error: Still recording"
                        continue

                    print "Start recording at: %s"\
                        % (datetime.datetime.now().time().isoformat())
                    recording = True

                    print "Recording clip: %s" % (data)
                    #camera.start_recording('/home/pi/flow.df/video.data', format='yuv')
                    #==========================modify bitrate here==============================
                    camera.start_recording('/home/pi/flow.df/video.h264', format='h264', bitrate=10000000)
                    #camera.annotate_frame_num = True
                    #(awb_test1, awb_test2) = camera.awb_gains
                    #camera.annotate_text = "(%f, %f)" %(awb_test1, awb_test2)
                    print "recording..."

                    if(paramsManual):
                        sleep(0.3)
                        ss, awb, g = currentParams
                        # re-run the set params function,
                        # camera forgets manual setting when video starts
                        setCamParams(ss, awb, g)

            elif(cmd == "STP"):
                # stop recording

                if not recording:
                    print "Error: Not recording"
                    continue

                print "Stop recording"
                recording = False
                camera.stop_recording()


            elif(cmd == "STR"):
                streaming = True
                print "recording for stream"
                streaming_socket = socket.socket()
                streaming_socket.connect((host, 8000))
                connection = streaming_socket.makefile('wb')
                try:
                    # Start a preview and let the camera warm up for 2 seconds
                    camera.start_preview()
                    sleep(2)
                    # Start recording, sending the output to the connection for 60
                    # seconds, then stop
                    camera.start_recording(connection, format='h264')
                    print "stream recording initialized"
                    #camera.wait_recording(60)
                    #camera.stop_recording()
                except socket.error:
                    print "streaming error, closing sockets"             
                    connection.close()
                    streaming_socket.close()
                    streaming = False

            elif(cmd == "STSTP"):
                # stop streaming

                if not streaming:
                    print "Error: Not streaming"
                    continue

                print "Stop streaming"
                #camera.stop_recording()
                #sleep(5)
                try:
                    print "closing sockets"
                    connection.close()
                    streaming_socket.close()
                    streaming = False
                    camera.stop_recording()
                except:
                    print "could not close straming sockets"

    except KeyboardInterrupt:
        s.close()
        camera.close()
        sys.exit(0)


if __name__ == "__main__":
    main(sys.argv)

