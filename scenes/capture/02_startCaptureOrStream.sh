#!/bin/bash
k=1
IFS="="
capture=$1 # set to 1 if capture, to 0 if stream desired
no=2 # number of frames for capture, time for streaming
sleepno=$((5+$no))

# first copy all important files
while read -r name value
do
  if [[ $name == *"clientip"* ]]; then
    cp ~/frejek/captureclient.py ~/frejek/rpi/$k/
    cp ~/frejek/streamingclient.py ~/frejek/rpi/$k/
    cp ~/frejek/global.cfg ~/frejek/rpi/$k/
    ((k++))
  fi
done < global.cfg

# call captureclient or streamingclient
if [ $capture -eq 1 ] 
then
  #gnome-terminal -e "python /home/frejek/capturinghost.py" --title CaptureHost
  host=$(
  while read -r name value
  do
    if [[ $name == *"hostip"* ]]; then
      echo -n "$value "
    fi
  done < $"global.cfg")
  ./99_pi_cmd.sh "python ~/flow.df/captureclient.py $host $no"
else
  # the streaming process must be done sequentially
  IFS="="

  while read -r name value
  do
    if [[ $name == *"clientip"* ]]; then
      echo $name
      gnome-terminal -e "python /home/student/frejek/streaminghost.py" --title StreamingHost
      gnome-terminal -e "ssh pi@$value -o IdentityFile=~/frejek/ssh/raspi python /home/student/frejek/streamingclient.py $no" --title "StreamingClient $value"
      sleep $sleepno
    fi
  done < $"global.cfg"
fi

