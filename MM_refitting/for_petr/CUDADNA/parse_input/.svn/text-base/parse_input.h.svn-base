#ifndef PARSE_INPUT_H_
#define PARSE_INPUT_H_

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define UNPARSED 0
#define PARSED 1
#define ERROR 2
#define OPT_MAX_LENGTH 256
#define KEY_NOT_FOUND -1
#define KEY_FOUND 0

typedef char input_string[OPT_MAX_LENGTH];

typedef struct {
	int state;
	int N_opts;
	int N_alloc;
	int N_unread_keys;
	input_string *keys;
	input_string *values;
	int *reads;
	input_string *unread_keys;
} input_file;

void loadInputFile(input_file *inp, const char *filename);
void loadInput(input_file *inp, FILE *desc);

int getInputKeyIndex(input_file *inp, const char *skey, int mandatory);
int getInputString(input_file *inp, const char *skey, char *dest, int mandatory);
int getInputInt(input_file *inp, const char *skey, int *dest, int mandatory);
int getInputUInt(input_file *inp, const char *skey, unsigned int *dest, int mandatory);
int getInputLLInt(input_file *inp, const char *skey, long long int *dest, int mandatory);
int getInputDouble(input_file *inp, const char *skey, double *dest, int mandatory);
int getInputFloat(input_file *inp, const char *skey, float *dest, int mandatory);
int getInputChar(input_file *inp, const char *skey, char *dest, int mandatory);

void getTrimmedString(const char *src, char *dest);

void setUnreadKeys(input_file *inp);
void cleanInputFile(input_file *inp);

#endif
