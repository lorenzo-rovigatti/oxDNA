#include "parse_input.h"
#include "../Utils.h"
#include "../Logger.h"
#include "../oxDNAException.h"
#include <algorithm>

using std::string;
using std::map;
using std::set;
using std::vector;

input_file::input_file() {
	true_values.insert("true");
	true_values.insert("1");
	true_values.insert("yes");
	true_values.insert("yup");
	true_values.insert("of course");

	false_values.insert("false");
	false_values.insert("0");
	false_values.insert("no");
	false_values.insert("nope");
	false_values.insert("are you crazy?");
}

void printInput(input_file *inp, char *filename) {
	OX_LOG(Logger::LOG_INFO, "Printing the input file as used by oxDNA in `%s'", filename);
	FILE *out = fopen(filename, "w");

	for(input_map::iterator it = inp->keys.begin(); it != inp->keys.end(); it++) {
		fprintf(out, "%s = %s\n", it->first.c_str(), it->second.value.c_str());
	}

	fclose(out);
}

void loadInputFile(input_file *inp, const char *filename) {
	FILE *inp_file = fopen(filename, "r");
	if(inp_file == NULL) {
	fprintf(stderr, "Input file '%s' not found\n", filename);
		inp->state = ERROR;
		return;
	}
	loadInput(inp, inp_file);
	fclose(inp_file);
	return;
}

void addCommandLineArguments(input_file *inp, int argc, char *argv[]) {
	string s_inp("");
	for(int i = 0; i < argc; i++) s_inp += string(argv[i]) + string("\n");
	addInput(inp, s_inp);
}

int _readLine(std::vector<string>::iterator &it, std::vector<string>::iterator &end, string &key, string &value) {
	string option(*it);
	option = Utils::trim(option);

	if (option.size() > 0) {
		std::vector<string> words = Utils::split(option, '=');	
		
		if (words.size() == 1) {
			OX_LOG (Logger::LOG_WARNING, "Malformed line '%s' found. Ignoring it", option.c_str());
			return NOTHING_READ;
		}

		if (words.size() > 2) {
			OX_LOG (Logger::LOG_WARNING, "Multiple `=' symbols found in line '%s'. Assuming you know what you are doing.", option.c_str());
			for (unsigned int i = 2; i < words.size(); i ++) words[1] += string("=") + words[i];
		}

		string my_key = Utils::trim(words[0]);
		string my_value = Utils::trim(words[1]);

		if(my_value[0] == '{') {
			// counts the number of open and closed curly brackets 
			size_t open = std::count(my_value.begin(), my_value.end(), '{');
			size_t close = std::count(my_value.begin(), my_value.end(), '}');

			int sum = (int)open - (int)close; 

			if (sum < 0) {
				OX_LOG (Logger::LOG_WARNING, "Malformed line '%s' found. Extra `}'.", option.c_str());
				return NOTHING_READ;
			}

			if (sum > 0) my_value += string("\n");

			while (sum > 0) {
				it++;
				
				if(it == end) throw oxDNAException ("Unclosed `{' at the end of file. Aborting");
				
				string new_line = string(*it);
				new_line = Utils::trim(new_line);

				int n_open = std::count(new_line.begin(), new_line.end(), '{');
				int n_closed = std::count(new_line.begin(), new_line.end(), '}');
			
				sum += n_open;
				sum -= n_closed;

				if(n_closed > 0) {
					int last_pos = new_line.find_last_of("}");
					string after_end = new_line.substr(last_pos + 1);
					if(after_end.size() > 0) throw oxDNAException("Found the string '%s' after a closing curly brace. You should either comment it or remove it. Aborting", after_end.c_str());
				}

				my_value += new_line;
				if(sum != 0) my_value += string("\n");
			}
		}

		key = my_key;
		value = my_value;
	}
	else return NOTHING_READ;

	return KEY_READ;
}

void loadInput(input_file *inp, FILE *inp_file) {
	inp->state = UNPARSED;
	addInput(inp, inp_file);
	inp->state = PARSED;
}

void addInput(input_file *inp, std::string s_inp) {
	vector<string> tot_lines = Utils::split(s_inp, '\n');
	vector<string> lines;

	for(vector<string>::iterator it = tot_lines.begin(); it != tot_lines.end(); it++) {
		// remove in-line comments
		size_t comment_start = it->find('#');
		if(comment_start != string::npos) it->erase (comment_start, it->size() - comment_start);

		// split the string using ; as a delimiter
		std::vector<string> to_add = Utils::split(*it, ';');

		lines.insert(lines.end(), to_add.begin(), to_add.end());
	}

	std::vector<string>::iterator l_end = lines.end();
	for(std::vector<string>::iterator it = lines.begin(); it != lines.end(); it++) {
		string key, value;
		int res = _readLine(it, l_end, key, value);

		if(res == KEY_READ){
			input_value new_value(value);

			input_map::iterator old_val = inp->keys.find(key);
			if(old_val != inp->keys.end()) OX_LOG(Logger::LOG_WARNING, "Overwriting key `%s' (`%s' to `%s')", key.c_str(), old_val->second.value.c_str(), value.c_str());
			inp->keys[key] = value;
		}
	}
}

void addInput(input_file *inp, FILE *inp_file) {
	size_t alloc_size;
	char *c_option = NULL;
	string file_contents("");

	int count = 0;
	while(count != -1) {
		count = getline(&c_option, &alloc_size, inp_file);
		if(count != -1) file_contents += string(c_option);
		free(c_option);
		c_option = NULL;
	}

	addInput(inp, file_contents);
}

input_map::iterator getInputValue(input_file *inp, const char *skey, int mandatory)  {
	std::map<string, input_value>::iterator it = inp->keys.find(string(skey));
	if(it != inp->keys.end()) it->second.read++;
	else if(mandatory) throw oxDNAException("Mandatory key `%s' not found, exiting", skey);

	return it;
}

int getInputString(input_file *inp, const char *skey, std::string &dest, int mandatory) {
	input_map::iterator it = getInputValue(inp, skey, mandatory);
	if(it == inp->keys.end()) return KEY_NOT_FOUND;
		
	dest = it->second.value;

	return KEY_FOUND;
}

int getInputString(input_file *inp, const char *skey, char *dest, int mandatory) {
	string s_dest;
	int res = getInputString(inp, skey, s_dest, mandatory);
	if(res != KEY_FOUND) return res;
	strncpy(dest, s_dest.c_str(), sizeof(char) * (s_dest.size()+1));

	return KEY_FOUND;
}

int getInputInt(input_file *inp, const char *skey, int *dest, int mandatory) {
	input_map::iterator it = getInputValue(inp, skey, mandatory);
	if(it == inp->keys.end()) return KEY_NOT_FOUND;

	*dest = (int) floor(atof(it->second.value.c_str())+0.1);

	return KEY_FOUND;
}

int getInputBool(input_file *inp, const char *skey, bool *dest, int mandatory) {
	input_map::iterator it = getInputValue(inp, skey, mandatory);
	if(it == inp->keys.end()) return KEY_NOT_FOUND;

	// make it lower case
	string val = it->second.value;
	std::transform(val.begin(), val.end(), val.begin(), ::tolower);

	set<string>::iterator res = inp->true_values.find(val);
	if(res != inp->true_values.end()) *dest = true;
	else {
		res = inp->false_values.find(val);
		if(res != inp->false_values.end()) *dest = false;
		else throw oxDNAException("boolean key `%s' is invalid (`%s'), aborting.", skey, val.c_str());
	}

	return KEY_FOUND;
}

int getInputBoolAsInt(input_file *inp, const char *skey, int *dest, int mandatory) {
	bool ret;
	int result = getInputBool(inp, skey, &ret, mandatory);
	if(result == KEY_FOUND) *dest = (int) ret;

	return result;
}

int getInputLLInt(input_file *inp, const char *skey, long long int *dest, int mandatory) {
	input_map::iterator it = getInputValue(inp, skey, mandatory);
	if(it == inp->keys.end()) return KEY_NOT_FOUND;
	*dest = (long long) floor(atof(it->second.value.c_str())+0.1);

	return KEY_FOUND;
}

int getInputUInt(input_file *inp, const char *skey, unsigned int *dest, int mandatory) {
	input_map::iterator it = getInputValue(inp, skey, mandatory);
	if(it == inp->keys.end()) return KEY_NOT_FOUND;

	*dest = (unsigned int) floor (atof(it->second.value.c_str())+0.1);

	return KEY_FOUND;
}

int getInputDouble(input_file *inp, const char *skey, double *dest, int mandatory) {
	input_map::iterator it = getInputValue(inp, skey, mandatory);
	if(it == inp->keys.end()) return KEY_NOT_FOUND;

	*dest = atof(it->second.value.c_str());

	return KEY_FOUND;
}

int getInputFloat(input_file *inp, const char *skey, float *dest, int mandatory) {
	input_map::iterator it = getInputValue(inp, skey, mandatory);
	if(it == inp->keys.end()) return KEY_NOT_FOUND;

	*dest = atof(it->second.value.c_str());

	return KEY_FOUND;
}

int getInputChar(input_file *inp, const char *skey, char *dest, int mandatory) {
	input_map::iterator it = getInputValue(inp, skey, mandatory);
	if(it == inp->keys.end()) return KEY_NOT_FOUND;

	*dest = it->second.value[0];

	return KEY_FOUND;
}

template<typename number>
int getInputNumber(input_file *inp, const char *skey, number *dest, int mandatory) {
	input_map::iterator it = getInputValue(inp, skey, mandatory);
	if(it == inp->keys.end()) return KEY_NOT_FOUND;

	*dest = (number) atof(it->second.value.c_str());

	return KEY_FOUND;
}
template int getInputNumber(input_file *inp, const char *skey, float *dest, int mandatory);
template int getInputNumber(input_file *inp, const char *skey, double *dest, int mandatory);

void setUnreadKeys(input_file *inp) {
	for(input_map::iterator it = inp->keys.begin(); it != inp->keys.end(); it++) {
		if(it->second.read == 0) inp->unread_keys.push_back(it->first);
	}
}

int getInputKeys (input_file *inp, string begins_with, vector<string> * dest, int mandatory) {
	int ret = 0;
	bool add_all = false;
	if (begins_with.size() == 0) add_all = true;
	for (input_map::iterator it = inp->keys.begin(); it != inp->keys.end(); it++) {
		if (add_all || it->first.compare(0, begins_with.length(), begins_with) == 0) {
			dest->push_back (it->first);
			ret ++;
		}
	}
	
	if(mandatory && ret == 0) throw oxDNAException("At least one key starting with `%s' is required. Found none, exiting", begins_with.c_str());

	return ret;
}

void cleanInputFile(input_file *inp) {
	inp->state = UNPARSED;
}

