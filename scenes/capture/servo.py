#!/usr/bin/env python3


'''
This is a simple control application for the servos on the box.
The application itself writes the control ommands for the servoblaster on stdout so it
can be directly forwarded by netcat.

The servo IDs (Channel IDs of servoblaster are hardcoded here)
'''

import sys
from time import sleep
import tkinter as tk


# Channel IDs for the servos
channel_upper = 0
channel_lower = 2
channel_side = 4

idle_timers = {}

# Tuples with (last, now)
servo_values = {}

# Helper function to print on stderr
def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)


def idle_timer_func(id):
	# The idea is to turn the servo back by 1% to stop it from making noise.
	l, c = servo_values[id]
	n = c
	if l > c:
		n = c + 2
	elif l < c:
		n = c - 2
	set_channel(id, n, False)
	eprint("IDLE ", id)


def end_idle_timer(id):
	if not id in idle_timers:
		return
	t = idle_timers[id]
	if t is not None:
		root.after_cancel(t)


def start_idle_timer(id):
	end_idle_timer(id)
	idle_timers[id] = root.after(500, idle_timer_func, id)


def set_channel(id, percent, startIdleTimer=True):
	last_val = 0
	if id in servo_values:
		_, last_val = servo_values[id]

	d = ((float(percent) * (218 - 80)) / 100 + 80)#+205-80
	print("%d=%d" % (id, d))
	sys.stdout.flush()
	eprint("set %d to %d" %(id, d))

	servo_values[id] = (last_val, percent)
	if (startIdleTimer):
		start_idle_timer(id)


class App:
	def __init__(self, master, initValues=True, useGUI=False):
		# Add a delay here to avoid all servos statrting at the same time
		if initValues:
			  self.update_upper(0, useGUI)
			  sleep(0.5)
			  self.update_lower(0, useGUI)
			  sleep(0.5)

		if useGUI:
			frame = tk.Frame(master)
			frame.pack()
			tk.Label(frame, text="<---- CLOSE ---------------------------------------- OPEN ---->").grid(row=0, column=1)
			upper = tk.Scale(frame, from_=0, to=100, length=300,
				  orient=tk.HORIZONTAL, command=self.update_upper)
			upper.grid(row=1, column=1)
			upper.set(0)
			tk.Label(frame, text="\nTop").grid(row=1)

			lower = tk.Scale(frame, from_=0, to=100, length=300,
				  orient=tk.HORIZONTAL, command=self.update_lower)
			lower.grid(row=2, column=1)
			lower.set(0)
			tk.Label(frame, text="\nBottom").grid(row=2)

			#smoke = tk.Scale(frame, from_=0, to=100, length=300,
			#	  orient=tk.HORIZONTAL, command=self.update_smoke)
			#smoke.grid(row=3, column=1)
			#smoke.set(50)
			#tk.Label(frame, text="\nSide").grid(row=3)


	def update_upper(self, percent, useGUI):
		set_channel(channel_upper, int(percent), useGUI)
		#set_channel(channel_upper, int(percent) * 0.92 + 8)

	def update_lower(self, percent, useGUI):
		set_channel(channel_lower, int(percent), useGUI)

	def update_smoke(self, percent, useGUI):
		set_channel(channel_side, (100 - int(percent) * 0.8), useGUI)

	def fill(self, useGUI):
		#set_channel(channel_side, int(100))
		set_channel(channel_upper, int(0), useGUI)
		set_channel(channel_lower, int(0), useGUI)

	def plume(self, percentLower, useGUI):
		#set_channel(channel_side, int(0))
		set_channel(channel_lower, int(percentLower), useGUI)
		set_channel(channel_upper, int(70), useGUI)
		

#eprint("Start servo client")
useGUI = False
if useGUI:
	root = tk.Tk()
	root.wm_title('Servo Control')
else:
	root = 0
app = App(root, len(sys.argv) < 2, useGUI)
if useGUI: root.geometry("380x155+0+0")
if len(sys.argv) >= 2:
	if sys.argv[1] == 'fill':
		app.fill(useGUI)
	elif sys.argv[1] == 'plume':
		if len(sys.argv) < 3: print('not enough arguments passed to servo.py with option \'plume\': %s'%sys.argv)
		app.plume(sys.argv[2],useGUI)
	else: 
		print('argument passed to servo.py invalid: %s'%sys.argv[1])
else:
	print('not enough arguments passed to servo.py: %s, running gui.'%sys.argv)
	if useGUI: root.mainloop()
#eprint('Number of arguments:', len(sys.argv), 'arguments.')
#eprint('Argument List:', str(sys.argv) )
#root.mainloop()
