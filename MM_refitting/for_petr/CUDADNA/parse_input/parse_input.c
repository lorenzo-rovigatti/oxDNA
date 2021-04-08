#include "parse_input.h"

void getTrimmedString(const char *src, char *dest) {
	int start = 0;
	int end = strlen(src)-1;
	if(end < 0) {
		dest[0] = '\0';
		return;
	}

	while(isspace(src[start])) start++;
	while(isspace(src[end])) end--;

	strncpy(dest, src + start, end+1);
	//if(dest[end+1] != '\0') dest[end+1] = '\0';
	dest[end+1] = '\0';
}

void loadInputFile(input_file *inp, const char *filename) {
	FILE *inp_file = fopen(filename, "r");
	if(inp_file == NULL) {
	fprintf(stderr, "Input file '%s' not found\n", filename);
		inp->state = ERROR;
		return;
	}
	loadInput (inp, inp_file);
	fclose (inp_file);
	return;
}

void loadInput(input_file *inp, FILE *inp_file) {
	int eof, i;
	size_t alloc_size;
	char *delim, *option = NULL;
	char t_option[512];

	inp->keys = NULL;
	inp->values = NULL;
	inp->reads = NULL;
	inp->unread_keys = NULL;
	inp->N_unread_keys = 0;
	inp->state = UNPARSED;

	inp->N_alloc = inp->N_opts = eof = 0;
	while(!eof) {
		unsigned int count = getline(&option, &alloc_size, inp_file);
		if(count == -1) {
			free (option);
			option = NULL;
			eof = 1;
			continue;
		};
		if(strlen(option) == 0) continue;

		if(strlen(option) > 0 && option[strlen(option)-1] == '\n') option[strlen(option)-1] = '\0';
		getTrimmedString(option, t_option);

		if(strlen(t_option) > 0 && t_option[0] != '#') {
			delim = strchr(t_option, '=');
			if(delim == NULL) {
				fprintf(stderr, "WARNING: malformed line '%s'\n", option);
				continue;
			}
			*delim = '\0';

			// we could need a realloc
			if(inp->N_opts == inp->N_alloc) {
				inp->N_alloc += 10;
				inp->keys = (input_string *) realloc(inp->keys, inp->N_alloc * sizeof(input_string));
				inp->unread_keys = (input_string *) realloc(inp->unread_keys, inp->N_alloc * sizeof(input_string));
				inp->values = (input_string *) realloc(inp->values, inp->N_alloc * sizeof(input_string));
				inp->reads = (int *) realloc(inp->reads, inp->N_alloc * sizeof(int));
				for(i = inp->N_alloc-10; i < inp->N_alloc; i++) inp->reads[i] = 0;
			}

			// split the option in key and value and trim them both
			getTrimmedString(t_option, inp->keys[inp->N_opts]);
			getTrimmedString(delim+1, inp->values[inp->N_opts]);

//			printf("%d %s %s\n", ret, inp->keys[inp->N_opts], inp->values[inp->N_opts]);

			inp->N_opts++;
		}

		free(option);
		option = NULL;
	}

	inp->state = PARSED;
}

int getInputKeyIndex(input_file *inp, const char *skey, int mandatory)  {
	int i;
	for(i = 0; i < inp->N_opts; i++)
		if(!strncmp(skey, inp->keys[i], OPT_MAX_LENGTH)) {
		   if (strlen(inp->values[i]) > 0) {
			inp->reads[i]++;
			return i;
		}
		}

	if(mandatory) {
		fprintf(stderr, "Mandatory key '%s' not found, exiting\n", skey);
		exit(-1);
	}
	else fprintf(stderr, "INFO: optional key '%s' not found\n", skey);

	return KEY_NOT_FOUND;
}

int getInputString(input_file *inp, const char *skey, char *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	if(key != KEY_NOT_FOUND) strncpy(dest, inp->values[key], OPT_MAX_LENGTH);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputInt(input_file *inp, const char *skey, int *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	//if(key != KEY_NOT_FOUND) *dest = atoi(inp->values[key]);
	if(key != KEY_NOT_FOUND) *dest = (int) floor (atof(inp->values[key])+0.1);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputLLInt(input_file *inp, const char *skey, long long int *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	//if(key != KEY_NOT_FOUND) *dest = atol(inp->values[key]);
	if(key != KEY_NOT_FOUND) *dest = (long long) floor (atof(inp->values[key])+0.1);
	else return KEY_NOT_FOUND;

	return KEY_FOUND;
}

int getInputUInt(input_file *inp, const char *skey, unsigned int *dest, int mandatory) {
	int key = getInputKeyIndex(inp, skey, mandatory);
	//if(key != KEY_NOT_FOUND) *dest = atoi(inp->values[key]);
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
