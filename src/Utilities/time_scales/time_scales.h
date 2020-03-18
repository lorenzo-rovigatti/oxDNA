/**
 * @file    time_scales.h
 *
 */

#ifndef TIME_SCALES_H_
#define TIME_SCALES_H_

#define TS_LIN 0
#define TS_LOG_LIN 1
#define TS_NOT_SELECTED -1

#define TS_LOADED 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * @brief Incapsulates a time scale. It is used by oxDNA to manage output. It will probably be slowly integrated into ObservableOutput.
 *
 * Timescales are managed in plain c for historical reasons (old code ported to oxDNA)
 */
typedef struct {
	long long int next_step;
	int type;
	int interval;
	/// used by TS_LOG_LIN
	double alpha;
	int cycle;
	int current_point;
	int points_per_cycle;
	int state;
} time_scale;

void initTimeScale(time_scale *ts, int time_scale_type);
/// used by TS_LOG_LIN
void setTSPPC(time_scale *ts, int ppc);
void _computeAlpha(time_scale *ts);
/// used by TS_LIN and TS_LOG_LIN
void setTSInterval(time_scale *ts, int interval);
void setTSNextStep(time_scale *ts);
void setTSInitialStep(time_scale *ts, long long int);
void cleanTimeScale(time_scale *ts);
int isTSLoaded(time_scale *ts);
void dumpTSState(time_scale *ts, char *filename);
void loadTSState(time_scale *ts, char *filename);

#endif
