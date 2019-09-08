#define _XOPEN_SOURCE 600

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <wiringPi.h>

#define PIN_PULSE     16
#define PIN_DIRECTION 20
#define PIN_ENDSWITCH 21

#define DIRECTION_FORWARD 0
#define DIRECTION_BACKWARD 1

#define MICROSTEP 8
#define SLIDE_LENGTH 20000
#define RAMP_STEPS (300 * MICROSTEP)
#define RAMP_OFFSET (10 * MICROSTEP)

// 10 kHz max
const uint32_t DRIVE_DELAY = 50;

// 1 kHz
const uint32_t HOME_DELAY = 500;

// Reduced speed in the last part
// This avoids wobbling
const uint32_t DRIVE_DELAY_SLOW_ADD = 150;
const long SLOW_AREA_START = 0;//(16000 * MICROSTEP);

/*
 * Returns nonzero if the switch is pressed.
 */
static inline uint8_t endswitch_pressed()
{
	return digitalRead(PIN_ENDSWITCH) != 0;
}

/*
 * microsecond delay.
 * Waits until the current time is time + us
 * At call time should be the current time,
 * after return time is the current time again (after the delay).
 */
static inline void us_delay(uint32_t us, struct timespec *time)
{
	time_t end_sec = time->tv_sec;
	time_t end_ns = time->tv_nsec;
	end_ns += us * 1000;

	// Check / fix overflow
	while (end_ns >= 1000000000)
	{
		end_ns -= 1000000000;
		end_sec += 1;
	}

	// Loop until specified time elapsed.
	do
	{
		clock_gettime(CLOCK_MONOTONIC, time);
	}
	while (time->tv_sec < end_sec || (time->tv_sec == end_sec && time->tv_nsec <= end_ns));
}

/*
 * Performs a single step with a delay of delay_us.
 */
static inline void singlestep(uint32_t delay_us, struct timespec *time)
{
	// set pin
	digitalWrite(PIN_PULSE, 1);
	us_delay(delay_us, time);

	// clear pin
	digitalWrite(PIN_PULSE, 0);
	us_delay(delay_us, time);
}

/*
 * Performs <distance> steps into the direction <direction>
 * min_delay specifies the minimum delay when stepping.
 * The speed at start and end is slower.
 * The start pos is only used to reduce the speed at the end of the slide.
 */
static uint8_t drive(long distance, uint32_t min_delay, long start_pos)
{
	int8_t direction = 0;
	if (distance > 0)
	{
		digitalWrite(PIN_DIRECTION, DIRECTION_FORWARD);
		direction = 1;
	}
	else
	{
		digitalWrite(PIN_DIRECTION, DIRECTION_BACKWARD);
		distance = -distance;
		direction = -1;
	}

	struct timespec tp;
	clock_gettime(CLOCK_MONOTONIC, &tp);

	for (long d = 0; d < distance; d++)
	{
		if ((direction == -1) && endswitch_pressed())
		{
			puts("EMERGENCY STOP, END SWITCH PRESSED");
			return 0;
		}

		long end_distance = (distance - d) - 1;
		if (end_distance > d)
		{
			end_distance = d;
		}

		uint32_t wait_us = min_delay;

		if (end_distance < RAMP_STEPS)
		{
			wait_us = min_delay * (RAMP_STEPS + RAMP_OFFSET) / (end_distance + RAMP_OFFSET);
		}

		start_pos += direction;
		if (start_pos > SLOW_AREA_START)
		{
			if (start_pos > SLOW_AREA_START + RAMP_STEPS)
			{
				wait_us += DRIVE_DELAY_SLOW_ADD;
			}
			else
			{
				wait_us += DRIVE_DELAY_SLOW_ADD * (start_pos - SLOW_AREA_START) / RAMP_STEPS;
			}
		}

		//printf("dist: %ld, w: %u\n", end_distance, wait_us);

		singlestep(wait_us, &tp);

	}
	return 1;
}


/*
 * Returns the slide to the home position and returns
 * the number of steps required to reach the switch.
 */
static long home(uint32_t delay)
{
	if (endswitch_pressed())
	{
		// at end -> make some steps to release the switch
		printf("Already at end position -> try to release the switch\n");
		drive(1000, DRIVE_DELAY, 0);
	}

	if (endswitch_pressed())
	{
		printf("Endswitch failed!\n");
		return -1;
	}

	// Set direction to backwards
	digitalWrite(PIN_DIRECTION, DIRECTION_BACKWARD);

	long counter = 0;

	// Time for delay function
	struct timespec tp;
	clock_gettime(CLOCK_MONOTONIC, &tp);

	// Step until switch reached
	while (!endswitch_pressed())
	{
		singlestep(delay, &tp);
		counter++;
	}

	printf("Reached home in %ld steps\n", counter);

	return counter;
}


static uint8_t setup(uint8_t drive_home)
{
	if (wiringPiSetupGpio() == -1)
	{
		printf("Setup failed!\n");
		return 0;
	}

	// Setup pins
	pinMode(PIN_PULSE, OUTPUT);
	pinMode(PIN_DIRECTION, OUTPUT);
	pinMode(PIN_ENDSWITCH, INPUT);
	pullUpDnControl(PIN_ENDSWITCH, PUD_UP);

	if (drive_home)
	{
		printf("Init done -> Return to home position\n");
		// First return to home position
		if (home(HOME_DELAY) < 0)
		{
			return 0;
		}
	}

	return 1;
}

#ifndef BUILD_LIB
// Build application for testing

/*
 * Reads the stdin and returns the parsed position value.
 */
static long get_next_position()
{
	for (;;)
	{
		printf("Enter position: ");
		char command_buffer[32];
		if (fgets(command_buffer, sizeof(command_buffer) - 1, stdin) == NULL)
		{
			puts("fgets() failed");
			continue;
		}
		command_buffer[sizeof(command_buffer) - 1] = 0;

		char *end;
		long pos = strtol(command_buffer, &end, 10);
		if (end == command_buffer)
		{
			printf("Invalid position: %s", command_buffer);
		}
		else if (pos > SLIDE_LENGTH)
		{
			printf("Position value too large: %l\n", pos);
		}
		else
		{
			return pos * MICROSTEP;
		}
	}
}

int main (int argc, char **argv)
{
	if (argc == 2 && strcmp(argv[1], "unstuck") == 0)
	{
		// Unstuck the slide in the case the endswitch failed
		setup(0);

		// Drive 2000 steps forward
		drive(2000, DRIVE_DELAY, 0);

		// Reset the direction to backwards
		digitalWrite(PIN_DIRECTION, DIRECTION_BACKWARD);
		return 0;
	}

	if (!setup(1))
	{
		return 1;
	}

	// This is my current position
	long pos = 0;

	printf("Start command loop\n");
	// Now start the command loop
	for (;;)
	{
		long next = get_next_position();
		if (next < 0)
		{
			printf("Reset postition\n");
			if (home(HOME_DELAY) < 0)
			{
				return 1;
			}
			pos = 0;

			continue;
		}

		long diff = next - pos;
		printf("Drive to %d (diff: %d)\n", next, diff);
		if (!drive(diff, DRIVE_DELAY, pos))
		{
			return 0;
		}
		pos = next;
	}
	return 0;
}

#else
// Build lib for using with python

long current_position = -1;

int init()
{
	if(setup(1))
	{
		current_position = 0;
		return 1;
	}
	return 0;
}

int set_position(long pos)
{
	if (current_position < 0)
	{
		printf("Stepper not initialized or in failsafe!\n");
		return 0;
	}

	if (pos < 0)
	{
		printf("Reset postition\n");
		if (home(HOME_DELAY) < 0)
		{
			current_position = -1;
			return 0;
		}
		current_position = 0;
		return 1;
	}

	if (pos > SLIDE_LENGTH)
	{
		printf("Position out of range!\n");
		return 0;
	}

	pos *= MICROSTEP;

	long diff = pos - current_position;
	printf("Drive to %d (diff: %d)\n", pos, diff);
	if (!drive(diff, DRIVE_DELAY, current_position))
	{
		return 0;
	}
	current_position = pos;
	return 1;
}

#endif
