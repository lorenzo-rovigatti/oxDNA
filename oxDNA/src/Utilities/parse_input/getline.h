/*
 * @file    getline.h
 * @date    31/gen/2012
 * @author  lorenzo
 */

/* @brief   getdelim & getline for uClibc (http://cristi.indefero.net/p/uClibc-cristi/source/tree/0_9_14/libc/stdio)
 *
 * Copyright (C) 2000 by Lineo, inc. and Erik Andersen
 * Copyright (C) 2000,2001 by Erik Andersen <andersen@uclibc.org>
 * Written by Erik Andersen <andersen@uclibc.org>
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Library General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef GETLINE_H_
#define GETLINE_H_

#ifdef __APPLE__

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#define __USE_GNU
#include <sys/types.h>

ssize_t getline(char **linebuf, size_t *n, FILE *file);
ssize_t getdelim(char **linebuf, size_t *linebufsz, int delimiter, FILE *file);

#endif

#endif /* GETLINE_H_ */
