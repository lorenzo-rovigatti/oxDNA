#include "timing.h"

void init_timer(LR_timer *timer, unsigned int num_events, msg *msgs) {
	unsigned int i;

	timer->num_events = num_events;
	timer->events = (clock_t *) calloc(2*num_events, sizeof(clock_t));
	timer->times_added = (long long int *) calloc(2*num_events, sizeof(long long int));
	timer->timings = (double *) calloc(num_events, sizeof(double));
	timer->msgs = (msg *) malloc(num_events * sizeof(msg));

	for(i = 0; i < num_events; i++) strncpy(timer->msgs[i], msgs[i], sizeof(msg));
}

void get_time(LR_timer *timer, unsigned int event) {
	if(event % 2) timer->timings[(event-1)/2] += (double)clock() - (double)timer->events[event-1];
	//if(event % 2) timer->timings[(event-1)/2] += (double)1 - (double)timer->events[event-1];
	else timer->events[event] = clock();
	timer->times_added[event]++;
}

void process_times(LR_timer *timer) {
	// we want millisec
	unsigned int i;

	for(i = 0; i < timer->num_events; i++) {
//		timer->timings[i] += (double)timer->events[2*i + 1]*ratio - (double)timer->events[2*i]*ratio;
		timer->events[2*i + 1] = timer->events[2*i] = 0;
	}
}

void prepare_timer_results(LR_timer *timer) {
	unsigned int i;
	double ratio = 1000. / CLOCKS_PER_SEC;

	for(i = 0; i < timer->num_events; i++) {
		if(timer->times_added[2*i] > 0) {
			timer->timings[i] *= ratio;
			timer->timings[i] /= (double) timer->times_added[2*i];
		}
	}
}

void divide_given_timing(LR_timer *timer, unsigned int event, double factor) {
	timer->timings[event] /= factor;
}

void divide_timer_results(LR_timer *timer, double factor) {
	unsigned int i;

	for(i = 0; i < timer->num_events; i++) {
			timer->timings[i] /= factor;
	}
}

void print_times(LR_timer *timer, FILE *out) {
	unsigned int i;

	for(i = 0; i < timer->num_events; i++) fprintf(out, "%s: %lf ms\n", timer->msgs[i], timer->timings[i]);
}

void destroy_timer(LR_timer *timer) {
	free(timer->events);
	free(timer->times_added);
	free(timer->timings);
	free(timer->msgs);
}
