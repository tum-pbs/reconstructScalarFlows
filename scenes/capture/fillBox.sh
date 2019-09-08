#!/bin/bash

# first copy controlSmoke file to raspi (raspi must be mounted)
cp ~/frejek/controlSmoke.py ~/frejek/rpi/smoky/scripts/

# execute controlSmoke script on raspi
ssh pi@rpi-smoky.ge.in.tum.de -i ~/frejek/ssh/raspi "python ~/scripts/controlSmoke.py $1"
