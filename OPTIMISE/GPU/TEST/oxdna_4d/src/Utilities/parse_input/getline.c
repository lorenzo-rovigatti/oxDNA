/*
 * getline.c
 *
 *  Created on: 31/gen/2012
 *      Author: lorenzo
 */

#ifdef __APPLE__

#include "getline.h"

/* Read up to (and including) a TERMINATOR from STREAM into *LINEPTR
   (and null-terminate it). *LINEPTR is a pointer returned from malloc (or
   NULL), pointing to *N characters of space.  It is realloc'd as
   necessary.  Returns the number of characters read (not including the
   null delimiter), or -1 on error or EOF.  */
ssize_t getdelim(char **linebuf, size_t *linebufsz, int delimiter, FILE *file) {
	static const int GROWBY = 80; /* how large we will grow strings by */

	int ch;
	int idx = 0;

	if (file == NULL || linebuf == NULL || linebufsz == NULL) {
	    errno = EINVAL;
	    return -1;
	}

	if (*linebuf == NULL || *linebufsz < 2) {
		*linebuf = malloc(GROWBY);
		if (!*linebuf) {
			errno = ENOMEM;
			return -1;
		}
		*linebufsz += GROWBY;
	}

	while (1) {
		ch = fgetc(file);
		if (ch == EOF)
			break;
		/* grow the line buffer as necessary */
		while (idx > *linebufsz-2) {
			*linebufsz += GROWBY;
			*linebuf = realloc(*linebuf, *linebufsz);
			if (!*linebuf) {
				errno = ENOMEM;
				return -1;
			}
		}
		(*linebuf)[idx++] = (char)ch;
		if ((char)ch == delimiter)
			break;
	}

	if (idx != 0)
	    (*linebuf)[idx] = 0;
	else if ( ch == EOF )
		return -1;
	return idx;
}

/* Basically getdelim() with the delimiter hard wired to '\n' */
ssize_t getline(char **linebuf, size_t *n, FILE *file) {
	return (getdelim (linebuf, n, '\n', file));
}

#endif
