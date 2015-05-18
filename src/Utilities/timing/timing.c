#include "timing.h"

void init_timer(LR_timer *timer, unsigned int num_events, msg *msgs) {
	unsigned int i;

	timer->state = TIMER_READY;
	timer->num_events = num_events;
	timer->events = (clock_t *) calloc(2*num_events, sizeof(clock_t));
	timer->times_added = (long long int *) calloc(2*num_events, sizeof(long long int));
	timer->timings = (double *) calloc(num_events, sizeof(double));
	timer->msgs = (msg *) malloc(num_events * sizeof(msg));

	for(i = 0; i < num_events; i++) strncpy(timer->msgs[i], msgs[i], sizeof(msg));
}

void get_time(LR_timer *timer, unsigned int event) {
	if(timer->state != TIMER_READY) return;
	if(event % 2) {
		double diff = (double)clock() - (double)timer->events[event-1];
		// we need this because of the idiocy of the standard regarding the clock_t type. Indeed,
		// the standard does not give any requirement on the size of this typedef, which can be
		// either an int, an unsigned int or a long long int. When clock()'s internal counter
		// (which is of type clock_t) reaches its maximum it resets to 0 (or to a negative value
		// if it is a signed variable). If this happens in between two calls to get_time, the
		// diff variable has a negative value which would screw up the associated timing. In
		// order to avoid this, if diff is negative we just skip ahead and disregard the event.
		// It's dirty but it works.
		if(diff < 0.) return;
		timer->timings[(event-1)/2] += diff;
	}
	else timer->events[event] = clock();
	timer->times_added[event]++;
}

void process_times(LR_timer *timer) {
	if(timer->state != TIMER_READY) return;
	unsigned int i;
	for(i = 0; i < timer->num_events; i++) {
		timer->events[2*i + 1] = timer->events[2*i] = 0;
	}
}

void prepare_timer_results(LR_timer *timer) {
	if(timer->state != TIMER_READY) return;
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
	if(timer->state != TIMER_READY) return;
	timer->timings[event] /= factor;
}

void divide_timer_results(LR_timer *timer, double factor) {
	if(timer->state != TIMER_READY) return;
	unsigned int i;

	for(i = 0; i < timer->num_events; i++) {
			timer->timings[i] /= factor;
	}
}

void print_times(LR_timer *timer, FILE *out) {
	if(timer->state != TIMER_READY) return;
	unsigned int i;

	for(i = 0; i < timer->num_events; i++) fprintf(out, "%s: %lf ms\n", timer->msgs[i], timer->timings[i]);
}

void destroy_timer(LR_timer *timer) {
	if(timer->state != TIMER_READY) return;

	free(timer->events);
	free(timer->times_added);
	free(timer->timings);
	free(timer->msgs);
	timer->state = TIMER_NOT_INIT;
}
