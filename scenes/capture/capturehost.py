#!/usr/bin/env python3

import argparse
import socket
import asyncio
import subprocess
import time, os, sys
from PIL import Image

# The default camera parameters
# None marks them as undefined
# if defined they are (shutter, awb_1, awb_2, gain)
defaultCamParams = (33164, 1.5, 1.5, 4.8)

# ID -> ClienConnection
# Ordinary clients with camera for capturing
camera_clients = {}

# Additional clients for room periphrals (like servos)
aux_clients = {}


class ClientMsgError(Exception):
	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return "ClientMsgError: {}".format(self.msg)

class ClientResultError(ClientMsgError):
	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return "ClientResultError: {}".format(self.msg)

class ClientNotFoundError(Exception):
	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return "ClientNotFoundError: {}".format(self.msg)

class ClientSocketError(Exception):
	pass

# An instance of this class handles a connection to a single client.
# It offsers async functions to execute commands on the client.
# These functions wait until the command finishes.
# It may be a good idea to split this class into subclasses for the different
# client types. For now, all functions are in one class.
class ClientConnection:
	def __init__(self, reader, writer):
		self._reader = reader
		self._writer = writer
		self._buffer = str()
		self._cid = None

	# This function should be called first, it gets the ID of the client by reading the HI message form the socket
	@asyncio.coroutine
	def client_id(self):
		if self._cid != None:
			return self._cid

		print("Waitng for client message")
		msg = yield from self._get_next_message()
		print(msg[0])
		d = msg.split('|')
		if len(d) != 2 or d[0] != 'HI':
			print("Unexpected message from client\
				   (Expected HI message):\n%s" % msg)
			self.remove_client()
			raise ClientMsgError(msg)

		if d[1].isdigit():
			self._cid = int(d[1])
			print("Camera client connected, id=%d" % self._cid)
		else:
			self._cid = d[1]
			print("Auxiliary client connected, name=%s" % self._cid)

		# confirm connection
		yield from self._send("CON|%s~" % d[1])

		return self._cid

	def remove_client(self):
		global camera_clients
		global aux_clients

		if self._cid in camera_clients.keys():
			del camera_clients[self._cid]
		if self._cid in aux_clients.keys():
			del aux_clients[self._cid]
		self._writer.close()


	# Motor client only: Set the position of the slide. pos is from 0 to 1
	@asyncio.coroutine
	def aux_set_slide(self, pos):
		yield from self._send("SET|%f~" % pos)
		yield from self._wait_for_ok("SET")

	# Set the camera resolution
	# mode = (w, h, fps)
	@asyncio.coroutine
	def set_resolution(self, mode):
		yield from self._send("RES|%d:%d:%d~" % mode)
		yield from self._wait_for_ok("RES")

	# Take a single image (and store it as <name>)
	@asyncio.coroutine
	def take_image(self, name):
		yield from self._send("SAM|%s~" % name)
		yield from self._wait_for_ok("SAM")

	# Starts the recording
	@asyncio.coroutine
	def start_video(self, name):
		yield from self._send("REC|%s~" % name)

	# Stops the recording
	@asyncio.coroutine
	def stop_video(self):
		yield from self._send("STP|~")

	# Starts the stream
	@asyncio.coroutine
	def start_stream(self):
		yield from self._send("STR|%s~")

	# Stops the stream
	@asyncio.coroutine
	def stop_stream(self):
		yield from self._send("STSTP|~")

	# Automatic camera parameters
	@asyncio.coroutine
	def param_auto(self):
		yield from self._send("CAP|~")
		yield from self._wait_for_ok("CAP")

	# Lock camera parameters to the current values
	@asyncio.coroutine
	def param_lock(self):
		yield from self._send("CFP|~")
		yield from self._wait_for_ok("CFP")

	# Set the camera parameters
	@asyncio.coroutine
	def param_set(self, params):
		yield from self._send("CSP|%i:%f:%f:%f~" % params)
		yield from self._wait_for_ok("CSP")

	# Get the current parameters form a camera
	@asyncio.coroutine
	def get_cam_params(self, mode):
		yield from self.param_auto()
		yield from self.set_resolution(mode)
		print("Resolution set! Wait 5 sec for exposure values to adapt")
		yield from asyncio.sleep(5)
		yield from self.param_lock()
		yield from self._send("CGP|~")
		status, data = yield from self._get_resposne()

		if status != "CP":
			print("Unexpected response to parameter query: %s" % status)
			raise ClientResultError(status)

		d = data.split(":")
		if (len(d) != 4 or (not isFloat(d[0])) or (not isFloat(d[1])) or
		  (not isFloat(d[2])) or (not isFloat(d[3]))):
			print("Invalid CP message: %s" % data)
			raise ClientMsgError(data)

		params = (int(d[0]), float(d[1]), float(d[2]), float(d[3]))
		return params

	# Internal: Read next message form the socket
	@asyncio.coroutine
	def _get_next_message(self):
		while True:
			delim = self._buffer.find('~')
			if delim != -1:
				res = self._buffer[:delim]
				self._buffer = self._buffer[delim + 1:]
				return res.strip('\r\n')
			try:
				data = yield from self._reader.read(100)
				self._buffer += data.decode('ascii')
			except Exception as e:
				print(e)
				self.remove_client()
				raise ClientSocketError()

	# Internal: Gets the response for a command
	@asyncio.coroutine
	def _get_resposne(self):
		msg = yield from self._get_next_message()
		delim = msg.find("|")
		if (msg == "" or delim == -1):
			# format invalid
			print("Received invalid response: %s" % msg)
			raise ClientMsgError(msg)

		# status ~ code of the message
		# data   ~ data of the message
		status = msg[0:delim]
		data = msg[delim + 1:]

		return (status, data)

	# Wait for a 'OK' result
	@asyncio.coroutine
	def _wait_for_ok(self, command):
		status, data = yield from self._get_resposne()
		if status != "OK" or data != command:
			raise ClientResultError(command)

	# Send a command
	@asyncio.coroutine
	def _send(self, msg):
		self._writer.write(msg.encode('ascii'))
		yield from self._writer.drain()


# Calls the given function with the given paramerters on all clients and waits for the results
@asyncio.coroutine
def command_all_clients(function, *args):
	coros = []
	for client in camera_clients.values():
		if client != None:
			coros += [function(client, *args)]
	yield from asyncio.gather(*coros)


# Take an image with all cameras
@asyncio.coroutine
def take_client_images(filename, resolution=None, params=None):
	if resolution != None:
		print("Set resolution")
		yield from command_all_clients(ClientConnection.set_resolution, resolution)

	if params != None:
		print("Set params")
		yield from command_all_clients(ClientConnection.param_set, params)

	print("Take image")
	yield from command_all_clients(ClientConnection.take_image, filename)


# Starts the capture on all cameras
@asyncio.coroutine
def start_capture(filename, resolution, params):
	#print("Set resolution")
	yield from command_all_clients(ClientConnection.set_resolution, resolution)

	#print("Set params")
	yield from command_all_clients(ClientConnection.param_set, params)

	#print("Starting video")
	yield from command_all_clients(ClientConnection.start_video, filename)


# Stops the capture
@asyncio.coroutine
def stop_capture():
	#print("Stop video")
	yield from command_all_clients(ClientConnection.stop_video)
	#print("--> DONE <--")


# Move the marker to a given position and wait for a given time
@asyncio.coroutine
def move_marker(pos, wait=5):
	print("Move marker to %f" % pos)

	if not "CS" in aux_clients.keys():
		raise ClientNotFoundError("Callibration silde client is not connected!")

	yield from aux_clients["CS"].aux_set_slide(pos)
	if wait > 0:
		print("OK, wait %i sec for marker to stop wobbling" % wait)
		yield from asyncio.sleep(wait)


# Execute the calibration
# Move the marker step-by-step and take images
@asyncio.coroutine
def do_calibration(resolution, steps, params):
	print("Start calibartion run")

	#params = (1000, 1, 1, 1.5);
	#defaultCamParams = (33164, 1.5, 1.5, 4.8)
	print("Set params")
	yield from command_all_clients(ClientConnection.param_set, params)

	print("Set resolution")
	yield from command_all_clients(ClientConnection.set_resolution, resolution)

	print("Now start moving marker")
	for p in range(steps + 1):
		cur = p / float(steps)
		yield from move_marker(cur)

		yield from take_client_images("/home/pi/flow.df/calib_%02d.jpg" % p)

	print("Returning marker to home position")
	yield from move_marker(0, 0)


# Starts the stream on a given client (UNTESTED)
@asyncio.coroutine
def do_streaming(client):
	#open socket
	listening_sock = socket.socket()
	listening_sock.bind(('0.0.0.0', 8000))
	listening_sock.listen(0)

	#send msg to cam to start stream
	yield from client.start_stream()

	#accept connection from camera
	connection = listening_sock.accept()[0].makefile('rb')

	print("streaming connection opened")
	try:
		#open player
		print("open player")
		cmdline = ['mplayer', '-fps', '90', '-cache', '1024', '-']
		player = subprocess.Popen(cmdline, stdin=subprocess.PIPE)
		print("player opened")
		while True:
			data = connection.read(1024)
			#print "data received"
			if not data:
				print("no data")
			break
		player.stdin.write(data)
	except:
		print("data reading and writing to mplayer failed")

	connection.close()
	listening_sock.close()
	player.terminate()
	client.stop_stream()
	print("Stream eneded")

# This is the 'main' function of this.
# In here all the CLI input processing is done and the commands are called.
@asyncio.coroutine
def handle_user_input(stdin, master_id, exposure_correction, cam_mode, loop=None):
	global camera_clients

	# the current camera parameters
	currentParamData = defaultCamParams

	while True:
		line = yield from stdin.readline()
		data = line.decode('ascii').strip('\r\n\t ')
		if data == 'h':
			# help
			print("c   Start capture")
			print("e   Get exposure values")
			print("ec  Set exposure correction value")
			print("f   Close / open valves to fill box with smoke")
			print("p   Take position image")
			print("s   Take sample image")
			print("q   quit")
			print("l   live video (requires mplayer)")
			print("cal start calibration")

		elif data.startswith('ec'):
			p = data.split(' ')
			if len(p) == 1:
				# no argument -> print current value
				print("exposure_correction = %f" % exposure_correction)

			else:
				# else: has argument -> check if valid
				if not isFloat(p[1]) or float(p[1]) <= 0:
					print("illegal value: %s" % p[1])
					continue

				# set new value and print
				exposure_correction = float(p[1])
				print("OK, exposure_correction = %f" % exposure_correction)

				# print resulting gain (if parameters are set)
				if currentParamData != None:
					(s, a1, a2, g) = currentParamData
					g = min(max(g * exposure_correction, 1), 12)
					print("gain is now = %f" % g)
					currentParamData = (s, a1, a2, g)

		elif data == 'p':
			# position image
			try:
				yield from take_client_images("/home/pi/flow.df/pos.jpg", (2592, 1944, 1))
			except (ClientMsgError, ClientSocketError) as err:
				print(err)

			print("--> DONE <--")

		elif data == 's':
			print("Take sample image")

			try:
				yield from take_client_images("/home/pi/flow.df/sample.jpg", cam_mode, currentParamData)
			except (ClientMsgError, ClientSocketError) as err:
				print(err)

			sampleFilename = '/home/student/frejek/rpi/%d/sample.jpg'
			images = []

			for i in range(0,5):
				im = Image.open(sampleFilename%(i+1))
				im = im.transpose(Image.ROTATE_90)
				#im.show()
				images.append(im)

			new_im = Image.new('L', (5*1080, 1920))

			x_offset = 0
			for im in images:
			  new_im.paste(im, (x_offset,0))
			  x_offset += im.size[0]

			new_im.save('/home/student/frejek/samples.png')
			size = 1.6*1080, 1.6*384
			new_im.thumbnail(size, Image.ANTIALIAS)
			new_im.show()

			print("--> DONE <--")

		elif data == 'e':
			print("Get exposure values")

			if (not master_id in camera_clients.keys()) or camera_clients[master_id] == None:
				print(camera_clients.keys())
				print("No client: %d (reference client id)" % master_id)
				continue

			try:
				currentParamData = yield from camera_clients[master_id].get_cam_params(cam_mode)
			except (ClientMsgError, ClientSocketError) as err:
				print(err)
				continue

			(s, a1, a2, g) = currentParamData
			g2 = min(max(g * exposure_correction, 1), 12)
			currentParamData = (s, a1, a2, g2)

			print("Camera params: shutter=%s, awb=(%s, %s), gain=%s, gain(corrected)=%f" % (s, a1, a2, g, g2))

		elif data.startswith('c'):
			arg = data.split(" ")
			if arg[0] == 'c':
				if len(arg) < 6: 
					print('C is missing additional arguments (videoFolder, percentageLower, numOfCaptures, smokeSecondsInit, smokeSecondsLoop)')
				else:
					
					print('')
					print('---------------------------------------------------------')
					print('Make sure there is enough smoke fluid in smoke machine!!!')
					print('---------------------------------------------------------')
					print('')

					#print("Release smoke and record videos.")
					if currentParamData == None:
						# capture without parameters set is not allowed
						print("Can't start capture without exposure data")
						continue

					print("Close valves.")
					subprocess.call(['./servocontrol.sh', 'fill', '0'])

					videoFolder = arg[1]
					percentageLower = arg[2]
					numOfCaptures = int(arg[3])
					smokeSecondsInit = arg[4]
					smokeSecondsLoop = arg[5]
					
					print("Fill box with smoke, init.")
					subprocess.call(['./fillBox.sh', smokeSecondsInit])
					print("Wait 150s for smoke to disappear and people leaving the room.")
					time.sleep(150)

					# loop for each capture process
					numberNaming = 0
					for capture in range(numOfCaptures):
						try:
							print("Start capture.")   
							yield from start_capture("rec.h264", cam_mode, currentParamData)
							time.sleep(5)
							print("Open valves.")
							subprocess.call(['./servocontrol.sh', 'plume', percentageLower])
						except (ClientMsgError, ClientSocketError) as err:
							print(err)

						print('Capture in progress.')
						time.sleep(18)

						try: 
							print("Stopping capture.")
							yield from stop_capture()
							print("Closing valves.")
							subprocess.call(['./servocontrol.sh', 'fill', '0'])
							print("Download videos.")
							while os.path.isdir(videoFolder+'_%s_%04d' % (percentageLower, numberNaming)):
								numberNaming = numberNaming+1
								if numberNaming>9999: 
									print('Folder exists: %s'%(videoFolder+'_%s_%04d' % (percentageLower, numberNaming)))
									break
							subprocess.call(['./grabVideos.sh', videoFolder+'_%s_%04d' % (percentageLower, numberNaming)])
						except (ClientMsgError, ClientSocketError) as err:
							print(err)

						print("Done with capture %04d."%capture)
						if capture+1<numOfCaptures:
							time.sleep(60)
							print("Close valves to fill box with smoke.")
							subprocess.call(['./servocontrol.sh', 'fill', '0'])
							print("Fill box with smoke, wait for smoke to disappear, loop.")
							subprocess.call(['./fillBox.sh', smokeSecondsLoop])
							time.sleep(180)

					print("Done with all captures 'c'.")
					print("Close valves.")
					subprocess.call(['./servocontrol.sh', 'fill', '0'])

			elif arg[0] == 'cal':
				if len(arg) < 2: 
					print('CAL is missing additional argument (calibFolder)')
				else:
					calibFolder = arg[1]
					try:
						yield from do_calibration((2592, 1944, 1), 20, currentParamData)
						print("--> DONE CAL <--")
						print("--> gather calibs <--")
						subprocess.call(['./gather_calib.sh', calibFolder])
						print("--> process images <--")
						subprocess.Popen(['./calibration/process_images.sh', calibFolder, '--ed=1000'])
						print("--> DONE process images <--")
					except (ClientMsgError, ClientSocketError, ClientNotFoundError) as err:
						print(err)					

			else:
				print("Unknown command: ", data)

		elif data.startswith('l'):
			arg = data.split(' ')
			if len(arg) == 1:
				# no argument -> print current value
				print("please add which client should be streamed")
				continue
			else:
				# else: has argument -> check if valid
				if not arg[1].isdigit():
					print("Not a decimal number: %s" % arg[1])
					continue

			client_num = int(arg[1])
			if not client_num in camera_clients.keys() or camera_clients[client_num] == None:
				print("No client with id: %d (reference client id)" % client_num)
				continue

			try:
				yield from do_streaming(camera_clients[client_num])
			except (ClientMsgError, ClientSocketError) as err:
				print(err)

		elif data == 'q':
			# q is same as CRTL+c
			raise KeyboardInterrupt
		else:
			print("Unknown command: ", data)


# checks if a value can be interpreted as a float
# fails with exponential values, but they do not occure here
def isFloat(str):
	return str.replace(".", "", 1).isdigit()


# Handles a connection from a new client
@asyncio.coroutine
def handle_connection(reader, writer):
	global camera_clients
	global aux_clients

	try:
		c = ClientConnection(reader, writer)
		new_client_id = yield from c.client_id()
		if type(new_client_id) is int:
			camera_clients[new_client_id] = c
		else:
			aux_clients[new_client_id] = c

	except (ClientMsgError, ClientSocketError) as e:
		print(e)
		print("Failed to establish connection")
		return

# Returns async reader for stdio
@asyncio.coroutine
def setup_stdio_reader(loop=None):
	if not loop:
		loop = asyncio.get_event_loop()
	reader = asyncio.StreamReader()
	reader_protocol = asyncio.StreamReaderProtocol(reader)

	yield from loop.connect_read_pipe(lambda: reader_protocol, sys.stdin)
	return reader


def main(argv):
	parser = argparse.ArgumentParser(
				  description='Host program for smoke capture camera system')

	parser.add_argument('-r', '--reference', type=int, default=1,
						help='Client id of the reference client (used for\
							  parameter synchronization)')

	parser.add_argument('-W', '--width', type=int, default=1296,
						help='Width (resolution) of the recorded video')

	parser.add_argument('-H', '--height', type=int, default=972,
						help='Height (resolution) of the recorded video')

	parser.add_argument('-F', '--fps', type=int, default=30,
						help='Frame rate of the recorded video')
	#parser.add_argument('-b', '--bind', default='131.159.40.51',
	#					help='Bind address used for the listening socket')
	parser.add_argument('-e', '--exposure_correction', default='0.8',
						type=float,
						help='Defines the initial value for the gain\
							  correction (this value is multilied with the\
							  current gain)')

	args = parser.parse_args()

	master_id = args.reference
	default_ec = args.exposure_correction
	cam_mode = (args.width, args.height, args.fps)

	loop = asyncio.get_event_loop()
	stdin_reader = loop.run_until_complete(setup_stdio_reader(loop=loop))
	input_coro = handle_user_input(stdin_reader, master_id, default_ec, cam_mode, loop=loop)
	server_coro = asyncio.start_server(handle_connection, '0.0.0.0', 54321, loop=loop)
	server = loop.run_until_complete(server_coro)

	try:
		loop.run_until_complete(input_coro)
	except KeyboardInterrupt:
		pass

	# Close the server
	server.close()
	loop.run_until_complete(server.wait_closed())
	loop.close()

if __name__ == "__main__":
	main(sys.argv)

