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
/* basic/math2.c --- Substitutes for missing math funcs PLUS some int math. */

/*--------------------------------------------------------------------------*/

#include "basic.h"

/*--------------------------------------------------------------------------*/

int basic_ipower (x, y)
     int x, y;
     /* Returns (int) x "to the power of" y, assuming y >= 0.
        Time complexity: O (ln(y)). */
{
  int aux;
  if (y == 0)
    return (1);
  else if (Odd (y))
    return (x * basic_ipower (x, y - 1));
  else
    {
      aux = basic_ipower (x, y / 2);
      return (aux * aux);
    }
}

/*--------------------------------------------------------------------------*/

#if defined (__sgi) || defined (__convex__) || defined (_IBMR2)

double log2 (x)
     double x;
{
  return (log (x) / log (2.0));
}

double exp2 (x)
     double x;
{
  return (exp (x * log (2.0)));
}

double exp10 (x)
     double x;
{
  return (exp (x * log (10.0)));
}

#endif  /* #ifdef */
