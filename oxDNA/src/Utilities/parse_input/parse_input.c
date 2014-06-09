#include "parse_input.h"

void printInput(input_file *inp, char *filename) {
	printf("%s\n", filename);
	FILE *out = fopen(filename, "w");

	int i;
	for (i = 0; i < inp->N_opts; i++) fprintf(out, "%s = %s\n", inp->keys[i], inp->values[i]);

	fclose(out);
}

void getTrimmedString(const char *src, char *dest) {
	int start = 0;
	int end = strlen(src)-1;
	if(end < 0) {
		dest[0] = '\0';
		return;
	}
	else if(end == 0) {
		dest[0] = src[0];
		dest[1] = '\0';
		return;
	}

	while(isspace(src[start]) && start <= end) start++;
	if (start < end) while(isspace(src[end]) && end >= 0) end--;

	strncpy(dest, src + start, sizeof(char) * (end - start + 1));
	dest[end-start+1] = '\0';
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
	FILE *temp = tmpfile();
	int i;
	for(i = 0; i < argc; i++) fprintf(temp, "%s\n", argv[i]);
	rewind(temp);

	addInput(inp, temp);

	fclose(temp);
}

int _readLine(FILE *inp_file, char **key, char **value) {
	size_t alloc_size;
	char *delim = NULL, *option = NULL;
	char t_option[OPT_MAX_LENGTH];
	int return_value = NOTHING_READ;

	int count = getline(&option, &alloc_size, inp_file);

	if(count == -1) {
		free(option);
		return INP_EOF;
	};

	if(strlen(option) > 0 && option[strlen(option)-1] == '\n') option[strlen(option)-1] = '\0';
	getTrimmedString(option, t_option);

	if(strlen(t_option) > 0 && t_option[0] != '#') {
		delim = strchr(t_option, '=');
		if(delim == NULL) {
			fprintf(stderr, "WARNING: malformed line '%s'\n", option);
			return_value = NOTHING_READ;
		}
		else {
			*delim = '\0';

			// split the option in key and value and trim them both
			*key = malloc(sizeof(char) * (strlen(t_option)+1));
			*value = malloc(sizeof(char) * (strlen(delim+1)+1));
			getTrimmedString(t_option, *key);
			getTrimmedString(delim+1, *value);

			if((*value)[0] == '{') {
				// check whether the parentheses are closed or not
				int i = 0;
				int sum = 0;
				for(i = 0; i < strlen(*value); i++) {
					if((*value)[i] == '{') sum++;
					else if((*value)[i] == '}') sum--;
				}
				if(sum > 0) {
					// we removed the trailing \n before, now we add it back again since we are dealing
					// with a multi-line option
					int len = strlen(*value);
					(*value)[len] = '\n';
					(*value)[len+1] = '\0';
					while(sum != 0) {
						char *new_line = NULL;
						int res = getline(&new_line, &alloc_size, inp_file);
						if(res == -1) {
							fprintf(stderr, "CRITICAL: Option '%s' has no closing curly bracket\n", *key);
							exit(1);
						}

						char * new_trimmed_line = malloc((strlen(new_line)+1)*sizeof(char));
						getTrimmedString(new_line, new_trimmed_line);
						int new_trimmed_len = strlen(new_trimmed_line);
						new_trimmed_line[new_trimmed_len] = '\n';
						new_trimmed_line[new_trimmed_len + 1] = '\0';
						if(new_trimmed_line[0] != '#') {
							for(i = 0; i < strlen(new_line); i++) {
								if(new_line[i] == '{') sum++;
								else if(new_line[i] == '}') sum--;
							}
							char *buffer = malloc((strlen(*value)+1) * sizeof(char));
							strcpy(buffer, *value);
							*value = realloc(*value, (strlen(*value) + strlen(new_trimmed_line) + 1) * sizeof(char));
							sprintf(*value, "%s%s", buffer, new_trimmed_line);
							free(buffer);
						}
						free(new_line);
						free(new_trimmed_line);
					}
				}
			}

			return_value = KEY_READ;
		}
	}
	else return_value = NOTHING_READ;

	free(option);

	return return_value;
}

void loadInput(input_file *inp, FILE *inp_file) {
	inp->keys = NULL;
	inp->values = NULL;
	inp->reads = NULL;
	inp->unread_keys = NULL;
	inp->N_unread_keys = 0;
	inp->state = UNPARSED;

	inp->N_alloc = inp->N_opts = 0;

	addInput(inp, inp_file);

	inp->state = PARSED;
}

void addInput(input_file *inp, FILE *inp_file) {
	int eof = 0;

	while(!eof) {
		char *key, *value;
		int res = _readLine(inp_file, &key, &value);

		if(res == INP_EOF) eof = 1;
		else if(res == KEY_READ){
			// we could need a realloc
			if(inp->N_opts == inp->N_alloc) {
				inp->N_alloc += 10;
				inp->keys = (input_string *) realloc(inp->keys, inp->N_alloc * sizeof(input_string));
				inp->unread_keys = (input_string *) realloc(inp->unread_keys, inp->N_alloc * sizeof(input_string));
				inp->values = (input_string *) realloc(inp->values, inp->N_alloc * sizeof(input_string));
				inp->reads = (int *) realloc(inp->reads, inp->N_alloc * sizeof(int));
				int i;
				for(i = inp->N_alloc-10; i < inp->N_alloc; i++) inp->reads[i] = 0;
			}

			if(strlen(value) < OPT_MAX_LENGTH) {
				int new_index = getInputKeyIndex(inp, key, 0);
				if(new_index == KEY_NOT_FOUND) {
					new_index = inp->N_opts;
					inp->N_opts++;
				}
				sprintf(inp->keys[new_index], "%s", key);
				sprintf(inp->values[new_index], "%s", value);
			}
			else {
				fprintf(stderr, "WARNING: option '%s' is too long (maximum length allowed is %d characters) and it will be ignored\n", key, OPT_MAX_LENGTH);
			}

			free(key);
			free(value);
		}
	}
}

int getInputKeyIndex(input_file *inp, const char *skey, int mandatory)  {
	int i;
	for (i = 0; i < inp->N_opts; i++) {
		if (!strncmp(skey, inp->keys[i], OPT_MAX_LENGTH)) {
			if (strlen(inp->values[i]) > 0) {
				inp->reads[i]++;
				return i;
			}
		}
	}

	if (mandatory) {
		fprintf(stderr, "Mandatory key '%s' not found, exiting\n", skey);
		exit(-1);
	}

	return KEY_NOT_FOUND;
}

int getInputString(input_file *inp, const char *skey, char *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) {
		strncpy(dest, inp->values[key], sizeof(char) * (strlen(inp->values[key])+1));
	}
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputInt(input_file *inp, const char *skey, int *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = (int) floor(atof(inp->values[key])+0.1);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputBoolAsInt(input_file *inp, const char *skey, int *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) {
		int val_len = strlen(inp->values[key]);
		char *lower_value = malloc((val_len+1) * sizeof(char));
		int i;
		for(i = 0; i < val_len; i++) lower_value[i] = tolower(inp->values[key][i]);
		lower_value[val_len] = '\0';

		*dest = -1;
		if(!strcmp(lower_value, "yes")) *dest = 1;
		if(!strcmp(lower_value, "no")) *dest = 0;
		if(!strcmp(lower_value, "true")) *dest = 1;
		if(!strcmp(lower_value, "false")) *dest = 0;
		if(!strcmp(lower_value, "1")) *dest = 1;
		if(!strcmp(lower_value, "0")) *dest = 0;
		free(lower_value);
		if(*dest == -1) return KEY_INVALID;
	}
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputLLInt(input_file *inp, const char *skey, long long int *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = (long long) floor(atof(inp->values[key])+0.1);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputUInt(input_file *inp, const char *skey, unsigned int *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = (unsigned int) floor (atof(inp->values[key])+0.1);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputDouble(input_file *inp, const char *skey, double *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = atof(inp->values[key]);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputFloat(input_file *inp, const char *skey, float *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = atof(inp->values[key]);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

/*
template <typename number>
int getInputNumber(input_file *inp, const char *skey, number *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = strtod(inp->values[key]);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

template <typename number>
int getInputDefaultNumber(input_file *inp, const char *skey, number *dest, number default_value) {
	int found = getInputNumber<number>(inp, skey, dest, 0);
	if(found == KEY_NOT_FOUND) *dest = default_value;

	return found;
}*/

int getInputChar(input_file *inp, const char *skey, char *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) *dest = inp->values[key][0];
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

void setUnreadKeys(input_file *inp) {
	int i;

	inp->N_unread_keys = 0;
	for(i = 0; i < inp->N_opts; i++) {
		if(inp->reads[i] == 0) {
			sprintf(inp->unread_keys[inp->N_unread_keys], "%s", inp->keys[i]);
			inp->N_unread_keys++;
		}
	}
}

void cleanInputFile(input_file *inp) {
	free(inp->keys);
	free(inp->values);
	free(inp->reads);
	free(inp->unread_keys);
	inp->state = UNPARSED;
}
