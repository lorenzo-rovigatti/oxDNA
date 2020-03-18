#include "time_scales.h"

void initTimeScale(time_scale * ts, int time_scale_type) {
	if(isTSLoaded(ts)) return;

	switch(time_scale_type) {
	case TS_LOG_LIN:
		ts->type = TS_LOG_LIN;
		break;
	case TS_LIN:
		ts->type = TS_LIN;
		break;
	default:
		fprintf(stderr, "Selected time scale %d not supported\n", time_scale_type);
		exit(-1);
	}

	ts->cycle = 0;
	ts->current_point = 0;
	ts->next_step = 0;
	ts->points_per_cycle = 0;
	ts->interval = 0;
}

void setTSPPC(time_scale *ts, int ppc) {
	if(isTSLoaded(ts)) return;

	ts->points_per_cycle = ppc;
	if(ts->type == TS_LOG_LIN) _computeAlpha(ts);
}

void setTSInterval(time_scale *ts, int interval) {
	if(isTSLoaded(ts)) return;

	ts->interval = interval;
	if(ts->type == TS_LOG_LIN) _computeAlpha(ts);
}

void _computeAlpha(time_scale *ts) {
	if(isTSLoaded(ts)) return;

	ts->alpha = ts->interval * pow(1.3, ts->points_per_cycle-1);
}

void setTSInitialStep(time_scale *ts, long long int initial) {
	if(isTSLoaded(ts)) return;

	if(ts->type == TS_LIN) ts->next_step = initial;
	switch(ts->type) {
	case TS_LOG_LIN:
		while(ts->next_step < initial) setTSNextStep(ts);
		break;
	case TS_LIN:
		ts->next_step = initial;
		break;
	default:
		fprintf(stderr, "You can't be here! Selected time scale %d not supported\n", ts->type);
		exit(1);
	}
}

void setTSNextStep(time_scale *ts) {
	switch(ts->type) {
	case TS_LOG_LIN:
		ts->next_step = (long long int)((ts->interval*pow(1.3, ts->current_point)) + ts->alpha * ts->cycle);
		ts->current_point++;
		if(ts->current_point == ts->points_per_cycle) {
			ts->cycle++;
			ts->current_point = 0;
		}
		break;
	case TS_LIN:
		ts->next_step += (long long int) ts->interval;
		break;
	default:
		fprintf(stderr, "You can't be here! Selected time scale %d not supported\n", ts->type);
		exit(-1);
	}
}

void cleanTimeScale(time_scale *ts) {
	ts->type = TS_NOT_SELECTED;
}

int isTSLoaded(time_scale *ts) {
	return (ts->state == TS_LOADED);
}

void dumpTSState(time_scale *ts, char *filename) {
	FILE *state = fopen(filename, "w");

	if(state == NULL) {
		fprintf(stderr, "Error opening file %s, exiting", filename);
		exit(-1);
	}

	fprintf(state, "%lld %d %d %lf %d %d %d\n", ts->next_step, ts->type, ts->interval, ts->alpha, ts->cycle, ts->current_point, ts->points_per_cycle);

	fclose(state);
}

void loadTSState(time_scale *ts, char *filename) {
	FILE *state = fopen(filename, "r");
	int tmpi;

	if(state == NULL) {
		fprintf(stderr, "File %s not found, exiting", filename);
		exit(-1);
	}

	tmpi = fscanf(state, "%lld %d %d %lf %d %d %d\n", &ts->next_step, &ts->type, &ts->interval, &ts->alpha, &ts->cycle, &ts->current_point, &ts->points_per_cycle);
	if (tmpi < 1) exit(-2);
	ts->state = TS_LOADED;

	fclose(state);
}
