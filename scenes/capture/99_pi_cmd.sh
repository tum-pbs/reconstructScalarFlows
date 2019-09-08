#!/bin/bash
IFS="="

clients=$(
echo -n "-H "
while read -r name value
do
  if [[ $name == *"clientip"* ]]; then
    echo -n "pi@$value "
  fi
done < $"global.cfg")

parallel-ssh $clients -t 7200 -x -2 -x -i -x ~/frejek/ssh/raspi -o ~/pssh_result -e ~/pssh_err $1

# pass command to execute it on each raspi
# -x -o -x IdentityFile=~/capture2/ssh/raspi

# fyi:
# -t timeout in seconds
# -x: extra ssh command line arguments
#    -i: display standard output and error as each host completes
#    -2: ssh protocol version 2 only
#    pass private key
# -o standard output directory
# -e error directory (same form as filenames for -o)
