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
/* basic/uhash.c  ---  Unversal hash functions. */

/*--------------------------------------------------------------------------*/

#include "basic.h"

/*--------------------------------------------------------------------------*/

#define max_a(R)  (MAXINT / powerof2 (bitsof (Basic_byte)) / (R))
        /* To avoid overflow w/i uhash_function(). */

/*--------------------------------------------------------------------------*/

void basic_uhash_new (s, r, m, a)
     int s, r;
     int *m;   /* output */
     int a[];  /* output: [0..r-1] */
     /* Given s, the desired size of the hash table, and r, the number of
        bytes per key, this procedure returns m, a prime number larger
        than s, and a vector a[0..r-1] defining the universal hash function
        It's the user responsibility to set up a hash table of size m, with
        r-bit keys.
        NOTE: Use a = MALLOC (int, r) and FREE (a) to allocate and free the
        memory of a[].  Also, call srandom (seed) to set up the seed for the
        random() function which will be called within this procedure.
        Reference, eg: Thomas H Cormen, Charles E Leiserson, and Ronald L
        Rivest. "Introduction to Algorithms."  MIT Press, 1990, p229ff. */
{
  int i, ma;
  *m = s = basic_prime_successor (s);
  ma = If ((s < max_a (r)), s - 1, max_a (r) - 1);
  upfor (i, 0, r - 1)
    a[i] = 1 + (random () mod ma);
}

/*--------------------------------------------------------------------------*/

int basic_uhash (a, r, m, x)
     int a[];  /* [0..r-1] */
     int r, m;
     Basic_byte x[];  /* [0..r-1] */
     /* Returns value of unversal hash function given by a[0..r-1] 
        for the (bytes of the) key x[0..r-1].  See basic_uhash_new(). */
{
  int i, sum = 0, m_reg = m;
  upfor (i, 0, r - 1)
    sum += (a[i] * x[i]) mod m_reg;
  return (sum mod m_reg);
}
