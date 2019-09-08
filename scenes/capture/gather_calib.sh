#!/bin/sh

IDS="1 2 3 4 5"
MOUNT_PATH="/home/student/frejek/rpi"

for i in $IDS
do
  mkdir $1
  mkdir $1/$i
  cp $MOUNT_PATH/$i/calib_*.jpg $1/$i&
done

for i in $IDS
do
  wait
done
