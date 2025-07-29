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

#ifdef JSON_ENABLED
#include <nlohmann/json_fwd.hpp>
#endif

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
	std::string key, value, expanded_value;
	std::vector<std::string> depends_on;
	int read = 0;

	input_value() : key(""), value("") {}
	input_value(std::string k, std::string v) : key(k), value(v) {}

	bool has_dependencies();
	void expand_value(std::map<std::string, std::string> expanded_dependency_values);
};

typedef std::map<std::string, input_value> input_map;

struct input_file {
	static input_file *main_input;
	static std::set<std::string> true_values;
	static std::set<std::string> false_values;

	input_map keys;
	std::vector<std::string> unread_keys;
	int state;
	bool show_overwrite_warnings = true;

	bool is_main_input = false;

	input_file(bool is_main);
	input_file();
	virtual ~input_file();

	/**
	 * @brief Initialise the structure from the given input file, given as a file descriptor.
	 * @param inp_file
	 */
	void init_from_file(FILE *inp_file);

	/**
	 * @brief Initialise the structure from the given input file, given as a filename.
	 * @param filename
	 */
	void init_from_filename(std::string filename);

	/**
	 * @brief Initialise the structure from the given string.
	 * @param s_inp string to be parsed
	 */
	void init_from_string(std::string s_inp);

#ifdef JSON_ENABLED
	void init_from_json(const nlohmann::json &json);
#endif

	/**
	 * @brief Load keys and values from command line's argc and argv variables
	 *
	 * This function will build the data structures out of argv. The first non-zero element should be
	 * the name of the input file, while the rest of the elements will be interpreted as key=value pairs.
	 * This method can optional skip additional argv elements that come after the input file and that
	 * should not interpreted as key=value pairs.
	 *
	 * @param argc
	 * @param argv
	 * @param args_to_skip optional parameter that sets the number of argv elements that should not be considered key=value pairs (in addition to the first one)
	 */
	void init_from_command_line_args(int argc, char *argv[], int args_to_skip=0);

	/**
	 * @brief Add the keys and values found in the desc file.
	 * @param desc input file to be parsed
	 */
	void add_input_source(FILE *desc);

	/**
	 * @brief Add the keys and values found in the string.
	 * @param s_inp string to be parsed
	 */
	void add_input_source(std::string s_inp);

	/**
	 * @brief Print out all the keys and values stored.
	 * @param filename
	 */
	void print(char *filename);

	void set_unread_keys();

	std::string get_value(std::string key, int mandatory, bool &found);
	void set_value(std::string key, std::string value);
	void unset_value(std::string key);

	std::string to_string() const;
};

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
int getInputKeys(input_file *inp, std::string begins_with, std::vector<std::string> * dest, int mandatory);

#endif
