#!/bin/bash

# This starts the "servo.py" script and a netcat instance to forward the data
# from the servo control GUI to the MotorPI where it is directly forwarded to
# the servoblaster software.


IFS="="

IP=$(
  while read -r name value
  do
    if [[ $name == "motorpi"* ]]; then
      echo -n "$value "
    fi
  done < $"global.cfg")

# 'fill' or 'plume', if 'plume', specify percentage of lower opening
./servo.py $1 $2 | nc $IP 1234
