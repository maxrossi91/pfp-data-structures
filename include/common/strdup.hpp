/* pfp-ds - prefix free parsing data structures
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file strdup.hpp
   \brief strdup.hpp contains implementation of strdup to override std. lib one in order to use malloc_count malloc.
   \author Massimiliano Rossi
   \date 12/03/2020
*/

#ifndef _STRDUP_HH
#define _STRDUP_HH

#include <malloc_count.h>
#include <string.h>
#include <stdlib.h>

extern char *strdup(const char *s) throw()
{
    size_t len = strlen (s) + 1;
    char *result = (char*) malloc (len);
    if (result == (char*) 0)
        return (char*) 0;
    return (char*) memcpy (result, s, len);
}

#endif /* end of include guard: _STRDUP_HH */