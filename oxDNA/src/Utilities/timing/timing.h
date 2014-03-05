/**
 * @file    timing.h
 *
 *
 */
#ifndef TIMING_H_
#define TIMING_H_

#define TIMER_NOT_INIT 0
#define TIMER_READY 1

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef char msg[256];

/**
 * @brief Manages timings. Useful to rate performances and for code optimising.
 *
 * Timings are managed in plain c for historical reasons (old code ported to oxDNA)
 */
typedef struct {
	int state;
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
