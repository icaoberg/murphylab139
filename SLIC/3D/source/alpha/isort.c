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
/* basic/isort.c  --- Inplace sort routines for very small int lists. */

/*--------------------------------------------------------------------------*/

#include "basic.h"

/*--------------------------------------------------------------------------*/

#define swap(A,B) \
do {              \
     aux = A;     \
     A = B;       \
     B = aux;     \
     swaps++;     \
   } once

/*--------------------------------------------------------------------------*/

int basic_isort2 (a, b)
     int *a, *b;  /* input/output */
     /* Sorts (&a,&b) with 1 comparison. */
{
  int swaps = 0, aux;
  if (*a > *b)
    swap (*a, *b);
  return (swaps);
}

/*--------------------------------------------------------------------------*/

int basic_isort3 (a, b, c)
     int *a, *b, *c;  /* input/output */
     /* Sorts (&a,&b,&c) with <= 3 comparisons and <= 3 swaps. */
{
  int swaps = 0, aux;
  if (*a > *b)
    swap (*a, *b);
  if (*b > *c)
    {
      swap (*b, *c);
      if (*a > *b)
        swap (*a, *b);
    }
  return (swaps);
}

/*--------------------------------------------------------------------------*/

int basic_isort4p (a, b, c, d)
     int *a, *b, *c, *d;  /* input/output */
     /* Sorts (&a,&b,&c,&d) with <= 3 comparisons and  <= 3 swaps,
        BUT assuming &a <= &b <= &c to begin with! */

{
  int swaps = 0, aux;
  Assert ((*a <= *b) and (*b <= *c));
  if (*d < *c)
    {
      if (*d < *a)
        {
          swap (*c, *d);
          swap (*b, *c);
          swap (*a, *b);
        }
      else if (*d < *b)
        {
          swap (*c, *d);
          swap (*b, *c);
        }
      else
        swap (*c, *d);
    }
  return (swaps);
}

/*--------------------------------------------------------------------------*/

int basic_isort4 (a, b, c, d)
     int *a, *b, *c, *d;  /* input/output */
     /* Sorts (&a,&b,&c,&d) with <= 6 = 3 + 3 comparisons
        and <= 6 = 3 + 3 swaps (insertion sort). */
{
  int swaps = 0, aux;
  /* step1: isort3 */
  if (*a > *b)
    swap (*a, *b);
  if (*b > *c)
    {
      swap (*b, *c);
      if (*a > *b)
        swap (*a, *b);
    }
  /* step2: isort4p */
  if (*d < *c)
    {
      if (*d < *a)
        {
          swap (*c, *d);
          swap (*b, *c);
          swap (*a, *b);
        }
      else if (*d < *b)
        {
          swap (*c, *d);
          swap (*b, *c);
        }
      else
        swap (*c, *d);
    }
  return (swaps);
}

/*--------------------------------------------------------------------------*/

int basic_isort5p (a, b, c, d, e)
     int *a, *b, *c, *d, *e;
     /* Sorts (&a,&b,&c,&d,&e) with <= 4 comparisons and <= 4 swaps,
        BUT assuming &a <= &b <= &c <= &d to begin with! */
{
  int swaps = 0, aux;
  if (*e < *d)
    {
      if (*e < *a)
        {
          swap (*d, *e);
          swap (*c, *d);
          swap (*b, *c);
          swap (*a, *b);
        }
      else if (*e < *b)
        {
          swap (*d, *e);
          swap (*c, *d);
          swap (*b, *c);
        }
      else if (*e < *c)
        {
          swap (*d, *e);
          swap (*c, *d);
        }
      else
        swap (*d, *e);
    }
  return (swaps);
}
