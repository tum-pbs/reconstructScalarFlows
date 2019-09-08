#!/bin/bash
fps=60
MOUNT_PATH="/home/student/frejek/rpi"

#if mounting point changes path for the video has to be changed
#similiarly if more cameras are used the loop has to be extended
mkdir $1
for i in 1 2 3 4 5
do
  echo $i
  MP4Box -fps $fps -add $MOUNT_PATH/$i/video.h264 $1/cam$i.mp4&
done

for i in 1 2 3 4 5
do
	wait
done
