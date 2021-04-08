#ifndef TIMING_H_
#define TIMING_H_

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef char msg[256];

typedef struct {
	unsigned int num_events;
	double *timings;
	clock_t *events;
	long long int *times_added;
	msg *msgs;
} LR_timer;

void init_timer(LR_timer *timer, unsigned int num_events, msg *msgs);
void get_time(LR_timer *timer, unsigned int event);
void process_times(LR_timer *timer);
void prepare_timer_results(LR_timer *timer);
void divide_given_timing(LR_timer *timer, unsigned int event, double factor);
void divide_timer_results(LR_timer *timer, double factor);
void print_times(LR_timer *timer, FILE *out);
void destroy_timer(LR_timer *timer);

#endif /* TIMING_H_ */
