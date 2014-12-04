#ifndef PARSE_INPUT_H_
#define PARSE_INPUT_H_

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <set>

#define UNPARSED 0
#define PARSED 1
#define ERROR 2
#define KEY_NOT_FOUND -1
#define KEY_FOUND 0
#define KEY_INVALID 1

#define INP_EOF -1
#define KEY_READ 0
#define NOTHING_READ 1

struct input_value {
	std::string value;
	int read;

	input_value() : value(""), read(0) {}
	input_value(std::string v) : value(v), read(0) {}
};

typedef std::map<std::string, input_value> input_map;
struct input_file {
	input_map keys;
	std::vector<std::string> unread_keys;
	int state;

	std::set<std::string> true_values;
	std::set<std::string> false_values;

	input_file();
};

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
 * @brief Add the keys and values found in the string to the given input_file.
 * @param inp target input_file structure
 * @param s_inp string to be parsed
 */
void addInput(input_file *inp, std::string s_inp);

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
 * @param it iterator begin position
 * @param end iterator end position
 * @param key
 * @param value
 * @return INP_EOF if EOF, NOTHING_READ if the line is malformed, empty or commented and KEY_READ otherwise.
 */
int _readLine(std::vector<std::string>::iterator &it, std::vector<std::string>::iterator &end, std::string &key, std::string &value); 

int getInputString(input_file *inp, const char *skey, std::string &dest, int mandatory);
int getInputString(input_file *inp, const char *skey, char *dest, int mandatory);
int getInputInt(input_file *inp, const char *skey, int *dest, int mandatory);
int getInputBoolAsInt(input_file *inp, const char *skey, int *dest, int mandatory);
int getInputBool(input_file *inp, const char *skey, bool *dest, int mandatory);
int getInputUInt(input_file *inp, const char *skey, unsigned int *dest, int mandatory);
int getInputLLInt(input_file *inp, const char *skey, long long int *dest, int mandatory);
int getInputDouble(input_file *inp, const char *skey, double *dest, int mandatory);
int getInputFloat(input_file *inp, const char *skey, float *dest, int mandatory);
template<typename number> int getInputNumber(input_file *inp, const char *skey, number *dest, int mandatory);
int getInputChar(input_file *inp, const char *skey, char *dest, int mandatory);

/**
 * @brief find all keys that begin with a substring
 *
 * @param inp input file
 * @param begins_with substring. If empty string, the function will return ALL the keys in the input file.
 * @param dest pointer to a vector of strings to be filled
 * @param mandatory whether at least one of these strings should be mandatory
 * @return the number of keys matching the string.
 */
int getInputKeys (input_file *inp, std::string begins_with, std::vector<std::string> * dest, int mandatory);

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
