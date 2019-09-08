#!/bin/bash
k=1
IFS="="

while read -r name value
do
  if [[ $name == *"clientip"* ]]; then
    mkdir -p ~/frejek/rpi/$k
    echo "mounting client $value as id $k"
    sshfs pi@$value:/home/pi/flow.df ~/frejek/rpi/$k -o IdentityFile=~/frejek/ssh/raspi 
    ((k++))
  fi
  if [[ $name == *"smoky"* ]]; then
    mkdir -p ~/frejek/rpi/smoky
    echo "mounting client $value as id smoky"
    sshfs pi@$value:/home/pi ~/frejek/rpi/smoky -o IdentityFile=~/frejek/ssh/raspi 
    ((k++))
  fi
done < global.cfg

