/*
 * Copyright (C) 2006 Murphy Lab,Carnegie Mellon University
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 * 
 * For additional information visit http://murphylab.web.cmu.edu or
 * send email to murphy@cmu.edu
 */
/* basic/gsort.c --- Quicksort routine: interface to standard qsort(). */

/*--------------------------------------------------------------------------*/

#include "basic.h"

/*---------------------------------------------------------------------------*/

void basic_qsort (table, i, j, compare)
     int  table[];  /* input/output */
     int  i, j;
     int (*compare)();
     /* This routine sorts table[i..j] in place: int (*compare)() is the
        comparison function, which is called with two arguments that point
        to the elements of table[i..j] that are compared; it is assumed to
        return -1, 0, or +1 with the usual meaning.  Cf, man qsort. */
{
  qsort (&(table[i]), j - i + 1, sizeof (int), compare);
}
