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

    d = ((float(percent) * (218 - 80)) / 100 + 80)
    print("%d=%d" % (id, d))
    sys.stdout.flush()
    eprint("set %d to %d" %(id, d))

    servo_values[id] = (last_val, percent)
    if (startIdleTimer):
        start_idle_timer(id)


class App:
    def __init__(self, master):
        # Add a delay here to avoid all servos statrting at the same time
        self.update_upper(0)
        sleep(0.5)
        self.update_lower(0)
        sleep(0.5)
        self.update_smoke(100)
        sleep(0.5)

        frame = tk.Frame(master)
        frame.pack()
        tk.Label(frame, text="<---- CLOSE ---------------------------------------- OPEN ---->").grid(row=0, column=1)
        upper = tk.Scale(frame, from_=0, to=100, length=300,
              orient=tk.HORIZONTAL, command=self.update_upper)
        upper.grid(row=1, column=1)
        upper.set(50)
        tk.Label(frame, text="\nTop").grid(row=1)

        lower = tk.Scale(frame, from_=0, to=100, length=300,
              orient=tk.HORIZONTAL, command=self.update_lower)
        lower.grid(row=2, column=1)
        lower.set(50)
        tk.Label(frame, text="\nBottom").grid(row=2)

        smoke = tk.Scale(frame, from_=0, to=100, length=300,
              orient=tk.HORIZONTAL, command=self.update_smoke)
        smoke.grid(row=3, column=1)
        smoke.set(50)
        tk.Label(frame, text="\nSide").grid(row=3)


    def update_upper(self, percent):
        set_channel(channel_upper, int(percent) * 0.92 + 8)

    def update_lower(self, percent):
        set_channel(channel_lower, int(percent))

    def update_smoke(self, percent):
        set_channel(channel_side, (100 - int(percent) * 0.8), False)


eprint("Start servo client")

root = tk.Tk()
root.wm_title('Servo Control')
app = App(root)
root.geometry("380x155+0+0")
#root.mainloop()
