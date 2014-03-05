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
#define OPT_MAX_LENGTH 4096
#define KEY_NOT_FOUND -1
#define KEY_FOUND 0
#define KEY_INVALID 1

#define INP_EOF -1
#define KEY_READ 0
#define NOTHING_READ 1

typedef char input_string[OPT_MAX_LENGTH];

/**
 * @brief Incapsulates a parsed input file.
 *
 * The parsing of the input file is carried out in plain c for historical reasons (old code ported to oxDNA)
 */
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

/**
 * @brief Load the keys and values found in the filename into the given input_file
 * This function is a simple wrapper around loadInput so that it is possible to accept a
 * filename instead of a file descriptor as input.
 * @param inp
 * @param filename
 */
void loadInputFile(input_file *inp, const char *filename);

/**
 * @brief Load the keys and values found in the desc file into the given input_file
 * @param inp
 * @param desc
 */
void loadInput(input_file *inp, FILE *desc);

/**
 * @brief Add the keys and values found in the desc file to the given input_file.
 * @param inp target input_file structure
 * @param desc input file to be parsed
 */
void addInput(input_file *inp, FILE *desc);

/**
 * @brief Print out all the keys and values stored in the given input_file structure.
 * @param inp
 * @param filename
 */
void printInput(input_file *inp, char *filename);

/**
 * @brief Add the keys and values stored in argv to the given input_file
 * @param inp target input_file structure
 * @param argc number of options
 * @param argv array of options
 */
void addCommandLineArguments(input_file *inp, int argc, char *argv[]);

/**
 * @brief Parse a line of the input file
 * This function allocates memory for the key and the value.
 * @param inp_file
 * @param key
 * @param value
 * @return INP_EOF if EOF, NOTHING_READ if the line is malformed, empty or commented and KEY_READ otherwise.
 */
int _readLine(FILE *inp_file, char **key, char **value);

int getInputKeyIndex(input_file *inp, const char *skey, int mandatory);
int getInputString(input_file *inp, const char *skey, char *dest, int mandatory);
int getInputInt(input_file *inp, const char *skey, int *dest, int mandatory);
int getInputBoolAsInt(input_file *inp, const char *skey, int *dest, int mandatory);
int getInputUInt(input_file *inp, const char *skey, unsigned int *dest, int mandatory);
int getInputLLInt(input_file *inp, const char *skey, long long int *dest, int mandatory);
int getInputDouble(input_file *inp, const char *skey, double *dest, int mandatory);
int getInputFloat(input_file *inp, const char *skey, float *dest, int mandatory);
int getInputChar(input_file *inp, const char *skey, char *dest, int mandatory);

/**
 * @brief Strip whitespace from the beginning and end of src
 * src is left unchanged and the resulting, trimmed string is stored in dest. This function
 * does not allocate memory for dest.
 * @param src
 * @param dest
 */
void getTrimmedString(const char *src, char *dest);

void setUnreadKeys(input_file *inp);

/**
 * @brief Clean up the given input_file structure
 * @param inp
 */
void cleanInputFile(input_file *inp);

#endif
