cd $(dirname "$0")
while [ 1 ]
do
	chrt --rr 99 python2 motorclient.py 131.159.40.51
	echo Client exited with code $?
	# Retry after 5 sec
	sleep 5
done
