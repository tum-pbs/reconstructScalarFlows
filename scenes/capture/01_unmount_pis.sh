#!/bin/bash
k=1
IFS="="

while read -r name value
do
  if [[ $name == *"clientip"* ]]; then
    mkdir -p ~/frejek/rpi/$k
    echo "unmounting client $value as id $k"
    fusermount -u ~/frejek/rpi/$k 
    rm -rf ~/frejek/rpi/$k 
    ((k++))
  fi
  if [[ $name == *"smoky"* ]]; then
    mkdir -p ~/frejek/rpi/smoky
    echo "unmounting client $value as id smoky"
    fusermount -u ~/frejek/rpi/smoky
    rm -rf ~/frejek/rpi/smoky
    ((k++))
  fi
done < global.cfg
