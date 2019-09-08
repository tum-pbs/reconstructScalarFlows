#!/bin/bash

# Configuration for the marker
MARKER="-d=7 -h=48 -w=29 --ml=0.01228125 --sl=0.02046875" 

# z-distance between two images in meters
ZDIST="0.02025"

if [ "$#" -lt 1 ]; then
    echo "USAGE: process_images.sh <IMG> [EXTRA AGRS]"
    echo "IMG is the path to the directory containing the images for all cameras. The images must be named calib_xx.jpg and placed subfolders named from 1 to 5 according to the cameras."
    echo "EXTRA ARGS are additional arguments that are passed to process_image_stack"
    exit
fi


#!/bin/sh

# Path to the "process_image_stack" executeable
# This is assumed in the same dir as this script
PIS=$(dirname $0)/process_image_stack

# Path to the folder containing the images
IMG=$1

# Consume the first arg
shift

PARAM="$MARKER -z=$ZDIST $@"

echo "args: $PARAM"

# $PIS $PARAM -v=$IMG/perfectMarker.png -o=1&
$PIS $PARAM -v=$IMG/2/calib_%02d.jpg -o=2&
$PIS $PARAM -v=$IMG/3/calib_%02d.jpg -o=3&

# Wait for a processes to finish before starting the next
# I don't want more than three running at the same time
wait

$PIS $PARAM -v=$IMG/1/calib_%02d.jpg -o=1&
$PIS $PARAM -v=$IMG/5/calib_%02d.jpg -o=5&
$PIS $PARAM -v=$IMG/4/calib_%02d.jpg -o=4&

# Wait for the rest
wait
wait
wait
