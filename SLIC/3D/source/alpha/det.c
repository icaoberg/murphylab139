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
/* lia/det.c  ---  Lia determinants. */

/* All procedures lia_detI (for I=2,3,4) calculate small determinants of Lia
   arguments employing expansion by minors (here: by the last column).  The
   expansion is implemented in a bottom-up (dynamic programming) fashion.
   According to the tables Harald Rosenberger's PhD thesis, this IS the best
   method to compute exact d-by-d determinants for d <= 5.  It clearly defeats
   Gaussian elimination for d <= 10. */

/* NOTE: In contrast to the basic routines of the Lia Library, lia_detI()'s
         output parameters are at the end of the parameter list! */

/*--------------------------------------------------------------------------*/

#include "basic.h"
#include "lia.h"
#include "base.h"

/*--------------------------------------------------------------------------*/

/* Workspace. */

#define W_SIZE  15

static Lia_ptr *w = NULL;
static int li = 0;

#define D12   w[5]
#define D13   w[6]
#define D14   w[7]
#define D23   w[8]
#define D24   w[9]
#define D34   w[10]
#define D123  w[11]
#define D124  w[12]
#define D134  w[13]
#define D234  w[14]

/*--------------------------------------------------------------------------*/

/* Macro for 2-by-2 determinant.
   RESULT = X11 * X22 - X12 * X21. */

#define DET2(X11,X12,X21,X22,RESULT) \
do { \
     lia_mul (w[0], X11, X22); \
     lia_mul (w[1], X12, X21); \
     lia_sub (RESULT, w[0], w[1]); \
   } once

/*--------------------------------------------------------------------------*/

/* Macro for expansion of 3-by-3 determinant by last row.
   RESULT = X13 * SUB23 - X23 * SUB13 + X33 * SUB12. */

#define EX3(X13,SUB23,X23,SUB13,X33,SUB12,RESULT) \
do { \
     lia_mul (w[0], X33, SUB12); \
     lia_mul (w[1], X23, SUB13); \
     lia_sub (w[3], w[0], w[1]); \
     lia_mul (w[0], X13, SUB23); \
     lia_add (RESULT, w[3], w[0]); \
   } once

/*--------------------------------------------------------------------------*/

void lia_det ()
     /* Module initialization: allocate the w[W_SIZE] working space. */
{
  if (not w)
    {
      int i;
      li = lia_get_maximum ();
      w = MALLOC (Lia_ptr, W_SIZE);
      MARK (w, -LIA_MAGIC);
      upfor (i, 0, W_SIZE)
        {
          w[i] = MALLOC (Lia, li);
          BZERO (w[i],   Lia, li);
          MARK  (w[i], -LIA_MAGIC);
        }
    }
}

/*--------------------------------------------------------------------------*/

void lia_det2 (a, b, 
               c, d,  result)
     Lia_obj a, b;
     Lia_obj c, d;
     Lia_obj result;  /* output */
{
  DET2 (a, b,  c, d,  result);
}

/*--------------------------------------------------------------------------*/

void lia_det3 (x11, x12, x13,
               x21, x22, x23,
               x31, x32, x33, result)
     Lia_obj x11, x12, x13;
     Lia_obj x21, x22, x23;
     Lia_obj x31, x32, x33;
     Lia_obj result;  /* output */
{
  DET2 (x11, x12,  x21, x22,  D12);
  DET2 (x11, x12,  x31, x32,  D13);
  DET2 (x21, x22,  x31, x32,  D23);
  ;
  EX3 (x13, D23,  x23, D13,  x33, D12,  result);
}

/*--------------------------------------------------------------------------*/

void lia_det4 (x11, x12, x13, x14,
               x21, x22, x23, x24,
               x31, x32, x33, x34,
               x41, x42, x43, x44, result)
     Lia_obj x11, x12, x13, x14;
     Lia_obj x21, x22, x23, x24;
     Lia_obj x31, x32, x33, x34;
     Lia_obj x41, x42, x43, x44;
     Lia_obj result;  /* output */
{
  DET2 (x11, x12,  x21, x22,  D12);
  DET2 (x11, x12,  x31, x32,  D13);
  DET2 (x11, x12,  x41, x42,  D14);
  DET2 (x21, x22,  x31, x32,  D23);
  DET2 (x21, x22,  x41, x42,  D24);
  DET2 (x31, x32,  x41, x42,  D34);
  ;
  EX3 (x13, D23,  x23, D13,  x33, D12,  D123);
  EX3 (x13, D24,  x23, D14,  x43, D12,  D124);
  EX3 (x13, D34,  x33, D14,  x43, D13,  D134);
  EX3 (x23, D34,  x33, D24,  x43, D23,  D234);
  ;
  lia_mul (w[0], x44, D123);
  lia_mul (w[1], x34, D124);
  lia_sub (w[3], w[0], w[1]);
  lia_mul (w[0], x24, D134);
  lia_mul (w[1], x14, D234);
  lia_sub (w[4], w[0], w[1]);
  lia_add (result, w[3], w[4]);  
}
