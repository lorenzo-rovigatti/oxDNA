#include "parse_input.h"

#include "../Utils.h"
#include "../Logger.h"
#include "../oxDNAException.h"
#include "../ConfigInfo.h"

#ifdef JSON_ENABLED
#include <nlohmann/json.hpp>
#endif

#include <exprtk/exprtk.hpp>

#include <regex>
#include <algorithm>

using std::string;
using std::map;
using std::set;
using std::vector;

input_file *input_file::main_input = nullptr;
std::set<std::string> input_file::true_values = {"true", "1", "yes", "yup", "of course"};
std::set<std::string> input_file::false_values = {"false", "0", "no", "nope", "are you crazy?"};

bool input_value::has_dependencies() {
	return !depends_on.empty();
}

void input_value::expand_value(std::map<std::string, std::string> expanded_dependency_values) {
	expanded_value = value;
	for(auto k : depends_on) {
		std::string to_sub = Utils::sformat("$(%s)", k.c_str());

		size_t pos = expanded_value.find(to_sub);
		expanded_value.replace(pos, to_sub.length(), expanded_dependency_values[k]);
	}

	// here we match patterns that are like this: ${some_text}, and we use parentheses to make sure that the second element of the std::smatch is "some_text"
	std::regex pattern("\\$\\{(.*)\\}"); // backslashes have to be escaped or the compiler complains
	std::smatch m;
	std::string to_search = expanded_value;
	exprtk::parser<double> parser;
	while(std::regex_search(to_search, m, pattern)) {
		exprtk::expression<double> expression;

		if(!parser.compile(m[1].str(), expression) || std::isnan(expression.value())) {
			throw oxDNAException("An error occurred during the evaluation of the mathematical expression '%s' required to expand the '%s' key", m[1].str().c_str(), key.c_str());
		}
		std::string expr_value = Utils::sformat("%lf", expression.value());

		size_t pos = expanded_value.find(m[0].str());
		expanded_value.replace(pos, m[0].str().length(), expr_value);

		to_search = m.suffix().str();
	}
}

input_file::input_file(bool is_main) :
				is_main_input(is_main) {
	if(is_main) {
		if(input_file::main_input != nullptr) {
			throw oxDNAException("There is already a main input file: Cannot initialise another one");
		}
		input_file::main_input = this;
	}

	state = UNPARSED;
}

input_file::input_file() :
				input_file(false) {

}

input_file::~input_file() {
	state = UNPARSED;
}

void input_file::init_from_file(FILE *inp_file) {
	state = UNPARSED;
	add_input_source(inp_file);
	state = PARSED;
}

void input_file::init_from_filename(std::string filename) {
	FILE *inp_file = fopen(filename.c_str(), "r");
	if(inp_file == NULL) {
		fprintf(stderr, "Input file '%s' not found\n", filename.c_str());
		state = ERROR;
		return;
	}
	init_from_file(inp_file);
	fclose(inp_file);
}

void input_file::init_from_string(std::string s_inp) {
	state = UNPARSED;
	add_input_source(s_inp);
	state = PARSED;
}

#ifdef JSON_ENABLED
void input_file::init_from_json(const nlohmann::json &json) {
	for(auto &item : json.items()) {
		try {
			std::string key(item.key());
			std::string value(item.value());
			set_value(key, value);
		}
		catch(nlohmann::detail::type_error &e) {
			// here we use a stringstream since if we are here it means that we cannot cast item.key() and/or item.value() to a string
			std::stringstream ss;
			ss << "The JSON line \"" << item.key() << " : " << item.value() << "\" cannot be converted to strings. ";
			ss << "Please make sure that all keys and values are quoted.";
			throw oxDNAException(ss.str());
		}
	}
}
#endif

void input_file::init_from_command_line_args(int argc, char *argv[], int args_to_skip) {
	init_from_filename(argv[1]);
	if(state == ERROR) {
		throw oxDNAException("Caught an error while opening the input file");
	}

	if(argc > (2 + args_to_skip)) {
		string s_inp("");
		for(int i = (2 + args_to_skip); i < argc; i++) {
			s_inp += string(argv[i]) + string("\n");
		}
		add_input_source(s_inp);
	}
}

void input_file::add_input_source(std::string s_inp) {
	vector<string> tot_lines = Utils::split(s_inp, '\n');
	vector<string> lines;

	for(vector<string>::iterator it = tot_lines.begin(); it != tot_lines.end(); it++) {
		// remove in-line comments
		size_t comment_start = it->find('#');
		if(comment_start != string::npos) {
			it->erase(comment_start, it->size() - comment_start);
		}

		// split the string using ; as a delimiter
		std::vector<string> to_add = Utils::split(*it, ';');

		lines.insert(lines.end(), to_add.begin(), to_add.end());
	}

	std::vector<string>::iterator l_end = lines.end();
	for(std::vector<string>::iterator it = lines.begin(); it != lines.end(); it++) {
		string key, value;
		int res = _readLine(it, l_end, key, value);

		if(res == KEY_READ) {
			set_value(key, value);
		}
	}
}

void input_file::add_input_source(FILE *inp_file) {
	size_t alloc_size;
	char *c_option = NULL;
	string file_contents("");

	int count = 0;
	while(count != -1) {
		count = getline(&c_option, &alloc_size, inp_file);
		if(count != -1) {
			file_contents += string(c_option);
		}
		free(c_option);
		c_option = NULL;
	}

	add_input_source(file_contents);
}

void input_file::print(char *filename) {
	OX_LOG(Logger::LOG_INFO, "Printing the input file as used by oxDNA in `%s'", filename);
	FILE *out = fopen(filename, "w");

	for(auto pair : keys) {
		fprintf(out, "%s = %s\n", pair.first.c_str(), pair.second.value.c_str());
	}

	fclose(out);
}

void input_file::set_unread_keys() {
	for(input_map::iterator it = keys.begin(); it != keys.end(); it++) {
		if(it->second.read == 0) {
			unread_keys.push_back(it->first);
		}
	}
}

std::string input_file::get_value(std::string key, int mandatory, bool &found) {
	found = false;
	std::map<string, input_value>::iterator it = keys.find(string(key));
	if(it != keys.end()) {
		input_value &value = it->second;
		value.read++;

		// this lambda recursively checks that there are no circular dependencies and expands all the values of the involved keys
		std::function<void(std::string,input_value &)> expand_value = [&expand_value](std::string root_key, input_value &current_value) {
			std::map<std::string, std::string> expanded_dependency_values;
			for(const auto &k : current_value.depends_on) {
				if(k == root_key) {
					throw oxDNAException("Circular dependency found between keys '%s' and '%s', aborting", root_key.c_str(), current_value.key.c_str());
				}
				auto dep_it = input_file::main_input->keys.find(k);
				if(dep_it == input_file::main_input->keys.end()) {
					throw oxDNAException("Key '%s' (which is expanded by '%s') is not defined", k.c_str(), current_value.key.c_str());
				}
				input_value &dep_value = input_file::main_input->keys[k];

				expand_value(root_key, dep_value);
				expanded_dependency_values[k] = dep_value.expanded_value;
			}
			current_value.expand_value(expanded_dependency_values);
		};

		// here we launch the recursive variable expansion
		expand_value(key, value);

		found = true;
		return value.expanded_value;
	}
	else if(mandatory) {
		throw oxDNAException("Mandatory key `%s' not found", key.c_str());
	}
	else {
		return "";
	}
}

void input_file::set_value(std::string key, std::string value) {
	if(key == "show_overwrite_warnings") {
		auto res = false_values.find(value);
		if(res != false_values.end()) {
			show_overwrite_warnings = false;
		}
	}

	input_map::iterator old_val = keys.find(key);
	if(old_val != keys.end() && show_overwrite_warnings) {
		OX_LOG(Logger::LOG_WARNING, "Overwriting key `%s' (`%s' to `%s')", key.c_str(), old_val->second.value.c_str(), value.c_str());
		keys[key].depends_on.clear();
	}
	keys[key] = input_value(key, value);

	// now we update the dependencies

	// here we match patterns that are like this: $(some_text), and we use parentheses to make sure that the second element of the std::smatch is "some_text"
	std::regex pattern("\\$\\(([\\w\\[\\]]+)\\)"); // backslashes have to be escaped or the compiler complains

	std::smatch m;
	while(std::regex_search(value, m, pattern)) {
		keys[key].depends_on.push_back(m[1].str());
		value = m.suffix().str();
	}
}

void input_file::unset_value(std::string key) {
	keys.erase(key);
}

std::string input_file::to_string() const {
	std::string ret;

	for(auto key : keys) {
		ret += Utils::sformat("%s = %s\n", key.first.c_str(), key.second.value.c_str());
	}

	// remove the last \n
	ret.pop_back();

	return ret;
}

int _readLine(std::vector<string>::iterator &it, std::vector<string>::iterator &end, string &key, string &value) {
	string option(*it);
	option = Utils::trim(option);

	if(option.size() > 0) {
		std::vector<string> words = Utils::split(option, '=');

		if(words.size() == 1) {
			OX_LOG(Logger::LOG_WARNING, "Malformed line '%s' found. Ignoring it", option.c_str());
			return NOTHING_READ;
		}

		if(words.size() > 2) {
			OX_LOG(Logger::LOG_WARNING, "Multiple `=' symbols found in line '%s'. Assuming you know what you are doing.", option.c_str());
			for(unsigned int i = 2; i < words.size(); i++) {
				words[1] += string("=") + words[i];
			}
		}

		string my_key = Utils::trim(words[0]);
		string my_value = Utils::trim(words[1]);

		if(my_value[0] == '{') {
			// counts the number of open and closed curly brackets 
			size_t open = std::count(my_value.begin(), my_value.end(), '{');
			size_t close = std::count(my_value.begin(), my_value.end(), '}');

			int sum = (int) open - (int) close;

			if(sum < 0) {
				OX_LOG(Logger::LOG_WARNING, "Malformed line '%s' found. Extra `}'.", option.c_str());
				return NOTHING_READ;
			}

			if(sum > 0) {
				my_value += string("\n");
			}

			while(sum > 0) {
				it++;

				if(it == end) {
					throw oxDNAException("Unclosed `{' at the end of file. Aborting");
				}

				string new_line = string(*it);
				new_line = Utils::trim(new_line);

				int n_open = std::count(new_line.begin(), new_line.end(), '{');
				int n_closed = std::count(new_line.begin(), new_line.end(), '}');

				sum += n_open;
				sum -= n_closed;

				if(n_closed > 0) {
					int last_pos = new_line.find_last_of("}");
					string after_end = new_line.substr(last_pos + 1);
					if(after_end.size() > 0)
						throw oxDNAException("Found the string '%s' after a closing curly brace. You should either comment it or remove it. Aborting", after_end.c_str());
				}

				my_value += new_line;
				if(sum != 0) {
					my_value += string("\n");
				}
			}
		}

		key = my_key;
		value = my_value;
	}
	else {
		return NOTHING_READ;
	}

	return KEY_READ;
}

int getInputString(input_file *inp, const char *skey, std::string &dest, int mandatory) {
	bool found;
	std::string tmp_value = inp->get_value(skey, mandatory, found);
	if(!found) {
		return KEY_NOT_FOUND;
	}

	// we don't want "dest" to be overwritten if the key was not found in the input file
	dest = tmp_value;

	return KEY_FOUND;
}

int getInputString(input_file *inp, const char *skey, char *dest, int mandatory) {
	string s_dest;
	int res = getInputString(inp, skey, s_dest, mandatory);
	if(res != KEY_FOUND) {
		return KEY_NOT_FOUND;
	}
	strncpy(dest, s_dest.c_str(), sizeof(char) * (s_dest.size() + 1));

	return KEY_FOUND;
}

int getInputInt(input_file *inp, const char *skey, int *dest, int mandatory) {
	bool found;
	std::string value = inp->get_value(skey, mandatory, found);
	if(!found) {
		return KEY_NOT_FOUND;
	}

	*dest = (int) std::floor(std::atof(value.c_str()) + 0.1);

	return KEY_FOUND;
}

int getInputBool(input_file *inp, const char *skey, bool *dest, int mandatory) {
	bool found;
	std::string value = inp->get_value(skey, mandatory, found);
	if(!found) {
		return KEY_NOT_FOUND;
	}

	// make it lower case
	std::transform(value.begin(), value.end(), value.begin(), ::tolower);

	set<string>::iterator res = inp->true_values.find(value);
	if(res != inp->true_values.end()) {
		*dest = true;
	}
	else {
		res = inp->false_values.find(value);
		if(res != inp->false_values.end()) {
			*dest = false;
		}
		else {
			throw oxDNAException("boolean key `%s' is invalid (`%s'), aborting.", skey, value.c_str());
		}
	}

	return KEY_FOUND;
}

int getInputBoolAsInt(input_file *inp, const char *skey, int *dest, int mandatory) {
	bool ret;
	int result = getInputBool(inp, skey, &ret, mandatory);
	if(result == KEY_FOUND) {
		*dest = (int) ret;
	}

	return result;
}

int getInputLLInt(input_file *inp, const char *skey, long long int *dest, int mandatory) {
	bool found;
	std::string value = inp->get_value(skey, mandatory, found);
	if(!found) {
		return KEY_NOT_FOUND;
	}

	*dest = (long long) std::floor(std::atof(value.c_str()) + 0.1);

	return KEY_FOUND;
}

int getInputUInt(input_file *inp, const char *skey, unsigned int *dest, int mandatory) {
	bool found;
	std::string value = inp->get_value(skey, mandatory, found);
	if(!found) {
		return KEY_NOT_FOUND;
	}

	*dest = (unsigned int) floor(atof(value.c_str()) + 0.1);

	return KEY_FOUND;
}

int getInputDouble(input_file *inp, const char *skey, double *dest, int mandatory) {
	bool found;
	std::string value = inp->get_value(skey, mandatory, found);
	if(!found) {
		return KEY_NOT_FOUND;
	}

	*dest = std::atof(value.c_str());

	return KEY_FOUND;
}

int getInputFloat(input_file *inp, const char *skey, float *dest, int mandatory) {
	bool found;
	std::string value = inp->get_value(skey, mandatory, found);
	if(!found) {
		return KEY_NOT_FOUND;
	}

	*dest = std::atof(value.c_str());

	return KEY_FOUND;
}

int getInputChar(input_file *inp, const char *skey, char *dest, int mandatory) {
	bool found;
	std::string value = inp->get_value(skey, mandatory, found);
	if(!found) {
		return KEY_NOT_FOUND;
	}

	*dest = value[0];

	return KEY_FOUND;
}

template<typename number>
int getInputNumber(input_file *inp, const char *skey, number *dest, int mandatory) {
	bool found;
	std::string value = inp->get_value(skey, mandatory, found);
	if(!found) {
		return KEY_NOT_FOUND;
	}

	*dest = (number) std::atof(value.c_str());

	return KEY_FOUND;
}
template int getInputNumber(input_file *inp, const char *skey, float *dest, int mandatory);
template int getInputNumber(input_file *inp, const char *skey, double *dest, int mandatory);

int getInputKeys(input_file *inp, string begins_with, vector<string> *dest, int mandatory) {
	int ret = 0;
	bool add_all = false;
	if(begins_with.size() == 0) {
		add_all = true;
	}
	for(input_map::iterator it = inp->keys.begin(); it != inp->keys.end(); it++) {
		if(add_all || it->first.compare(0, begins_with.length(), begins_with) == 0) {
			dest->push_back(it->first);
			ret++;
		}
	}

	if(mandatory && ret == 0) {
		throw oxDNAException("At least one key starting with `%s' is required. Found none, exiting", begins_with.c_str());
	}

	return ret;
}

