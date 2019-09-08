import time, sys
import RPi.GPIO as GPIO

# todo remove
print('in controlSmoke 1.')

### setup
GPIO.setwarnings(False)
GPIO.setmode(GPIO.BCM)
GPIO.setup(26, GPIO.IN, pull_up_down=GPIO.PUD_DOWN) # green LED
GPIO.setup(19, GPIO.IN, pull_up_down=GPIO.PUD_DOWN) # red LED
GPIO.setup(13, GPIO.OUT) # push button -> red LED is on

# todo remove
print('in controlSmoke 2.')

if len(sys.argv)<2:
	print('controlSmoke: not enough input arguments.')
else:
	smokeSeconds = int(sys.argv[1])
	print('Attempt to push smoke button for %d seconds. '%smokeSeconds)

	### push button
	didPush = False
	tooMuchTimePassed = False
	start_time = time.time()
	while not didPush and not tooMuchTimePassed:
		# hack! remove True once green LED works again
		if True or GPIO.input(26):
			GPIO.output(13, GPIO.HIGH)
			time.sleep(smokeSeconds)
			GPIO.output(13, GPIO.LOW)
			didPush = True
			print('Pushed the button for %d seconds!'%smokeSeconds)
		else: 
			time.sleep(5)
		tooMuchTimePassed = (time.time() - start_time) > 300
	if tooMuchTimePassed: print('Too much time has passed, green LED didn\'t switch on for a while.')
	if not didPush: print('Did not push the button.')

