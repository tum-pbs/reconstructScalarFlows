gcc -std=c99 -o stepper stepper.c -lwiringPi
gcc -std=c99 -shared -fPIC -Wl,-soname,stepper.so -o stepper.so -DBUILD_LIB stepper.c -lwiringPi
