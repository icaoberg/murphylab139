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
/*
 * weighted.c -- code for computing sizes of different simplices and
 *               attachment tests
 */

#ifndef lint
static char author[] = "Mike Facello";
#endif

/* #define __MDEBUG__  /* Uncomment this, except when debugging */

/*
 * date --  8/2/93
 *        11/30/93 -- size.c and attached.c combined
 *         2/10/94 -- merged into univeral (weighed & unweighted) mkalf
 *
 * NOTE: alf_w_begin() MUST be called before any functions in this module
 *       will work; and alf_w_close() should be called when finished.
 *
 */

/*----------------------------------------------------------------------*/

#include "mkalf.h"
#include "basic.h"
#include "sos.h"
#include "internal.h"

#define lia_square()       lia_ipower (2)

#ifdef __MDEBUG__
 static Alf_float value ();
#endif

/*----------------------------------------------------------------------*/

/* Global variables for partial result storage */

static Lia *results[5][5];
static Lia *d0, *d1, *d2, *d3, *d4, *m123;

/*----------------------------------------------------------------------*/

/* Global Variables for checking primitives */

/* maybe not needed anymore! */

static int math_test = FALSE;
static FILE *check_file;

/*--------------------------------------------------------------------------*/

/* counters */

static int size[4] = { 0 };
static int hidden[4] = { 0 };

/*----------------------------------------------------------------------*/

/*
 * alf_w_begin () -- Initializes the module by allocating partial result
 *                   storage.
 */
void alf_w_begin ()
{
  int a, b;

  upfor (a, 0, 3)
    size[a] = hidden[a] = 0;

  upfor (a, 0, 4)
    upfor (b, 0, 4)
      results[a][b] = MALLOC (Lia, lia_get_maximum());

  d0 = MALLOC (Lia, lia_get_maximum());
  d1 = MALLOC (Lia, lia_get_maximum());
  d2 = MALLOC (Lia, lia_get_maximum());
  d3 = MALLOC (Lia, lia_get_maximum());
  d4 = MALLOC (Lia, lia_get_maximum());
  m123 = MALLOC (Lia, lia_get_maximum());

  if (math_test)
    check_file = basic_fopen (basic_cb_frmt ("mkalf.math"), "w");
}

/*----------------------------------------------------------------------*/

/*
 * alf_w_end () -- Closes the module by freeing partial result
 *                 storage.
 */
void alf_w_end ()
{
  int a, b;

  /* Free up memory */

  upfor (a, 0, 4)
    upfor (b, 0, 4)
      FREE(results[a][b]);

  FREE(d0);
  FREE(d1);
  FREE(d2);
  FREE(d3);
  FREE(d4);
  FREE(m123);

  if (math_test)
    basic_fclose (check_file);
}

/*--------------------------------------------------------------------------*/

void alf_w_calls (size_cnt, hidden_cnt, dim)
     int **size_cnt;
     int **hidden_cnt;
     int *dim;
     /* Returns counters of corresponding function calls.
        See mkalf/mkalf.c for sample usage. */
        
{
  *size_cnt = size;
  *hidden_cnt = hidden;
  *dim = 3;
}

/*----------------------------------------------------------------------*
 *                                                                      *
 *                           Size Computations                          *
 *                                                                      *
 *----------------------------------------------------------------------*/

/*
 * alf_w_size0 (i, num, den) -- Computes the size of the vertex i.  The result
 *                       is stored in size.  The size is just the negative
 *                       weight.  To compute this, though, we have to
 *                       determine the weight from the fourth coordinate.
 *
 *                       A possible change:  store weight from input.
 */
void alf_w_size0 (i, num, den)
     int i;
     Lia_obj num, den;  /* Output */
{
  Assert (lia_stack_empty ());
  size[0]++;

  lia_push (sos_lia (i, 4));
  lia_push (sos_lia (i, 1));
  lia_square ();
  lia_minus ();
  lia_push (sos_lia (i, 2));
  lia_square ();
  lia_minus ();
  lia_push (sos_lia (i, 3));
  lia_square ();
  lia_minus ();

  lia_pop (num);
  lia_load(den, 1);

#ifdef __MDEBUG__
  {
    Alf_float v = (float) value(num, den);
    print ("Size0(%d) = ", i);
    print ("%f \n", v);
  }
#endif
}

/*----------------------------------------------------------------------*/

/*
 * alf_w_size1 (i, j, num, den) -- Computes the size of the edge with endpoints
 *                              i and j, storing the output num/den 
 *                              in den and num.
 *
 *       NOTE: Due to the way that the two additional hyperplanes are
 *             determined (by expanding parallel to the x2- and x3-axes),
 *             it is necessary to rotate the points by rotating the
 *             coordinates.  The array c[] contains the rotated
 *             coordinates.
 */
void alf_w_size1 (i, j, num, den)
     int i, j;
     Lia_obj num, den;  /* Output */
{
  int a, b;
  int coord[5];       /* Stores the coordinates after rotation */
  size[1]++;

#ifdef __MDEBUG2__
  if (sos_proto_flag)
    print ("size1 (%d,%d)\n", i, j);
#endif

  (void) basic_isort2 (&i, &j);

  /* First, determine proper ordering for the coordinates. */

  coord[4] = 4;
  if (not lia_eq (sos_lia (i, 1), sos_lia (j, 1)))
    {
      coord[1] = 1; coord[2] = 2; coord[3] = 3;
    }
  else if (not lia_eq (sos_lia (i, 2), sos_lia (j, 2)))
    {
      coord[1] = 2; coord[2] = 3; coord[3] = 1;
    }
  else if (not lia_eq (sos_lia (i, 3), sos_lia (j, 3)))
    {
      coord[1] = 3; coord[2] = 1; coord[3] = 2;
    }
  else
    basic_error ("size1: identical coordinates %d and %d.\n", i, j);


  /* Compute and store all the 2 x 2 subdeterminants needed. */

  /* when 0 is a column, it must come last.  results[0][b] stores
     M_{b,0} */

  upfor (b, 1, 4)
    lia_assign (results[0][b], sos_minor2 (i, j, coord[b], 0));

  upfor (a, 1, 3)
    upfor (b, a+1, 4)
      lia_assign (results[a][b], sos_minor2 (i, j, coord[a], coord[b]));

  /* Compute d0 */
  
  Assert (lia_stack_empty ());
  ;
  lia_push (results[0][1]);
  upfor (a, 1, 3) {
    lia_push (results[0][a]);
    lia_square ();
  }
  ;
  lia_plus ();
  lia_plus ();
  lia_times ();
  lia_push (lia_const (-2));
  lia_times ();
  ;
  lia_pop (d0);

  /* Compute d1 */
  
  Assert (lia_stack_empty ());

  lia_push (results[0][1]);

  lia_push (results[0][3]);
  lia_push (results[1][3]);
  lia_times ();
  lia_push (results[0][2]);
  lia_push (results[1][2]);
  lia_times();
  lia_plus ();
  lia_push (lia_const(2));
  lia_times ();

  lia_push (results[0][1]);
  lia_push (results[0][4]);
  lia_times ();
  lia_minus ();

  lia_times ();

  lia_pop (d1);

  /* Compute d2 */

  Assert (lia_stack_empty ());

  lia_push (results[1][2]);
  lia_push (results[0][1]);
  lia_square ();
  lia_push (results[0][3]);
  lia_square ();
  lia_plus ();
  lia_times ();
  lia_push (lia_const (-2));
  lia_times ();

  lia_push (results[0][2]);
  lia_push (results[0][1]);
  lia_push (results[0][4]);
  lia_times ();
  lia_push (results[0][3]);
  lia_push (results[1][3]);
  lia_push (lia_const (-2));
  lia_times ();
  lia_times ();
  lia_plus ();
  lia_times ();

  lia_minus ();

  lia_pop (d2);

  /* Compute d3 */

  Assert (lia_stack_empty ());

  lia_push (results[1][3]);
  lia_push (results[0][1]);
  lia_square ();
  lia_push (results[0][2]);
  lia_square ();
  lia_plus ();
  lia_times ();
  lia_push (lia_const (-2));
  lia_times ();

  lia_push (results[0][3]);
  lia_push (results[0][1]);
  lia_push (results[0][4]);
  lia_times ();
  lia_push (results[0][2]);
  lia_push (results[1][2]);
  lia_push (lia_const (-2));
  lia_times ();
  lia_times ();
  lia_plus ();
  lia_times ();

  lia_minus ();

  lia_pop (d3);

  /* Compute d4 */

  Assert (lia_stack_empty ());

  upfor (a, 1, 3) {
    lia_push (results[0][a]);
    lia_push (results[a][4]);
    lia_times ();
  }
  lia_plus ();
  lia_plus ();
  lia_push (results[0][1]);
  lia_push (lia_const (2));
  lia_times ();
  lia_times ();

  lia_push (results[0][3]);
  lia_push (results[1][2]);
  lia_times ();
  lia_push (results[0][2]);
  lia_push (results[1][3]);
  lia_times ();
  lia_minus ();
  lia_push (results[2][3]);
  lia_times ();

  lia_push (results[1][2]);
  lia_square ();
  lia_push (results[1][3]);
  lia_square ();
  lia_plus ();
  lia_push (results[0][1]);
  lia_times ();
  lia_minus ();
  lia_push (lia_const (4));
  lia_times ();

  lia_plus ();

  lia_pop (d4);

  /* Compute top of fraction */

  Assert (lia_stack_empty ());

  lia_push (d1);
  lia_square ();
  lia_push (d2);
  lia_square ();
  lia_push (d3);
  lia_square ();
  lia_plus ();
  lia_plus ();
  lia_push (d0);
  lia_push (d4);
  lia_times ();
  lia_minus ();

  lia_pop (num);

  /* Compute bottom of fraction */
  
  Assert (lia_stack_empty ());

  lia_push (d0);
  lia_square ();
  lia_pop (den);

  if (math_test)
      fprint (check_file, "testeq [size1[{%d, %d}], %.0f / %.0f]\n",
              i, j, lia_real(num), lia_real(den));

#ifdef __MDEBUG__
  {
    Alf_float v = (float) value(num, den);
    print ("Size1(%d, %d) = ", i, j);
    print ("%f = ", v);
    print ("%s / %s\n", lia_deci(num), lia_deci(den));
    lia_clear();
  }
#endif
}


/*----------------------------------------------------------------------*/

/*
 * alf_w_size2 (i, j, k, num, den) -- Computes the size of the edge with 
 *                              endpoints i and j, storing the size
 *                              num/den in den and num.
 */
void alf_w_size2 (i, j, k, num, den)
     int i, j;
     Lia_obj num, den;  /* Output */
{
  int a, b;
  size[2]++;
  
#ifdef __MDEBUG2__
  if (sos_proto_flag)
    print ("size2 (%d,%d,%d)\n", i, j, k);
#endif

  (void) basic_isort3 (&i, &j, &k);

  /* Determinant usage

     Used once:
     124,134,234

     Used more than once:
     120,130,140,230,240,340,123
   */

  /* Allocate and pre-compute determinants used more than once */

  upfor (a, 0, 2)
    upfor (b, a, 2)
      lia_assign (results[a][b], sos_minor3 (i, j, k, a+1, b+2, 0));

  lia_assign (m123, sos_minor3 (i, j, k, 1, 2, 3));

  /* Compute d0 */

  Assert (lia_stack_empty ());

  lia_push (results[0][0]);     /* M_{1,2,0} */
  lia_square ();
  lia_push (results[0][1]);     /* M_{1,3,0} */
  lia_square ();
  lia_push (results[1][1]);     /* M_{2,3,0} */
  lia_square ();
  lia_plus ();
  lia_plus ();
  lia_push (lia_const (4));
  lia_times ();

  lia_pop (d0);

  /* Compute d1 */

  Assert (lia_stack_empty ());

  lia_push (results[0][1]);     /* M_{1,3,0} */
  lia_push (results[2][2]);     /* M_{3,4,0} */
  lia_times ();
  lia_push (results[0][0]);     /* M_{1,2,0} */
  lia_push (results[1][2]);     /* M_{2,4,0} */
  lia_times ();
  lia_push (m123);
  lia_push (results[1][1]);     /* M_{2,3,0} */
  lia_push (lia_const (-2));
  lia_times ();
  lia_times ();
  lia_plus ();
  lia_plus ();

  lia_push (lia_const (-2));
  lia_times ();

  lia_pop (d1);

  /* Compute d2 */

  Assert (lia_stack_empty ());

  lia_push (results[0][0]);     /* M_{1,2,0} */
  lia_push (results[0][2]);     /* M_{1,4,0} */
  lia_times ();
  lia_push (results[1][1]);     /* M_{2,3,0} */
  lia_push (results[2][2]);     /* M_{3,4,0} */
  lia_times ();
  lia_minus ();
  lia_push (m123);
  lia_push (results[0][1]);     /* M_{1,3,0} */
  lia_push (lia_const (-2));
  lia_times ();
  lia_times ();
  lia_plus ();

  lia_push (lia_const (2));
  lia_times ();

  lia_pop (d2);

  /* Compute d3 */

  Assert (lia_stack_empty ());

  lia_push (results[1][1]);     /* M_{2,3,0} */
  lia_push (results[1][2]);     /* M_{2,4,0} */
  lia_times ();
  lia_push (results[0][1]);     /* M_{1,3,0} */
  lia_push (results[0][2]);     /* M_{1,4,0} */
  lia_times ();
  lia_push (m123);
  lia_push (results[0][0]);     /* M_{1,2,0} */
  lia_push (lia_const (2));
  lia_times ();
  lia_times ();
  lia_plus ();
  lia_plus ();

  lia_push (lia_const (2));
  lia_times ();

  lia_pop (d3);

  /* Compute d4 */

  Assert (lia_stack_empty ());

  lia_push (results[0][0]);             /* M_{1,2,0} */
  lia_push (sos_minor3 (i,j,k,1,2,4));  /* M_{1,2,4} */
  lia_times ();
  lia_push (results[0][1]);             /* M_{1,3,0} */
  lia_push (sos_minor3 (i,j,k,1,3,4));  /* M_{1,3,4} */
  lia_times ();
  lia_push (results[1][1]);             /* M_{2,3,0} */
  lia_push (sos_minor3 (i,j,k,2,3,4));  /* M_{2,3,4} */
  lia_times ();
  lia_push (m123);
  lia_square ();
  lia_push (lia_const (-2));
  lia_times ();
  lia_plus ();
  lia_plus ();
  lia_plus ();

  lia_push (lia_const (-4));
  lia_times ();

  lia_pop (d4);

  /* Compute top of fraction */

  Assert (lia_stack_empty ());

  lia_push (d1);
  lia_square ();
  lia_push (d2);
  lia_square ();
  lia_push (d3);
  lia_square ();
  lia_plus ();
  lia_plus ();
  lia_push (d0);
  lia_push (d4);
  lia_times ();
  lia_minus ();

  lia_pop (num);

  /* Compute bottom of fraction */
  
  Assert (lia_stack_empty ());

  lia_push (d0);
  lia_square ();
  lia_pop (den);

  if (math_test) {
      fprint (check_file, "testeq [size2[{%d, %d, %d}], %.0f / %.0f]\n",
              i, j, k, lia_real(num), lia_real(den));
/*      fprint (check_file, "d1 = %.0f,", lia_real(d1));
      fprint (check_file, "d2 = %.0f,", lia_real(d2));
      fprint (check_file, "d3 = %.0f,", lia_real(d3));
      fprint (check_file, "d4 = %.0f\n", lia_real(d4));
*/
    }

#ifdef __MDEBUG__
  {
    Alf_float v = (float) value(num, den);
    print ("Size2(%d, %d, %d) = ", i, j, k);
    print ("%f \n", v);
  }
#endif
}

/*----------------------------------------------------------------------*/

/*
 * alf_w_size3 (i, j, k, l, num, den) -- Computes the size of the tetrahedron 
 *                       ijkl.  The result is stored in num and den.
 */
void alf_w_size3 (i, j, k, l, num, den)
     int i, j, k, l;
     Lia_obj num, den;  /* Output */
{
  Assert (lia_stack_empty ());
  size[3]++;

  lia_push (sos_minor4(i,j,k,l,4,2,3,0));
  lia_square ();
  lia_push (sos_minor4(i,j,k,l,1,4,3,0));
  lia_square ();
  lia_push (sos_minor4(i,j,k,l,1,2,4,0));
  lia_square ();
  lia_plus ();
  lia_plus ();

  lia_push (sos_minor4(i,j,k,l,1,2,3,0));

  /* In the meantime, compute the bottom of the fraction */
  lia_pushtop ();
  lia_square ();
  lia_push(lia_const(4));
  lia_times ();
  lia_pop (den);

  lia_push (sos_minor4(i,j,k,l,1,2,3,4));
  lia_push (lia_const (4));
  lia_times ();
  lia_times ();
  lia_plus ();

  lia_pop (num);

  if (math_test)
      fprint (check_file, "testeq [size3[{%d, %d, %d, %d}], %.0f / %.0f]\n",
              i, j, k, l, lia_real(num), lia_real(den));

#ifdef __MDEBUG__
  {
    Alf_float v = (float) value(num, den);
    print ("Size3(%d, %d, %d, %d) = ", i, j, k, l);
    print ("%f = %s / %s\n", v, lia_deci(num), lia_deci(den));
  }
#endif
}


/*----------------------------------------------------------------------*
 *                                                                      *
 *                         Attachment Tests                             *
 *                                                                      *
 *----------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

/*
 * alf_w_hidden2 (i, j, k, p) -- Returns TRUE if triangle ijk is 
 *      attached with respect to p, i.e. if the power distance from p to 
 *      y_T is less than 0, where T = {i,j,k}.
 *
 * (Actually, this is exactly the same procedure as in_triangle_sphere.)
 */
int alf_w_hidden2 (i, j, k, p)
     /*i*/ int i, j, k, p;
{
  Lia results[3][4];
  int iota;
  Assert (lia_stack_empty ());
  hidden[2]++;

  /* Compute multiply-used determinants */

  lia_assign (results[0], sos_minor3 (i, j, k, 2, 3, 0));
  lia_assign (results[1], sos_minor3 (i, j, k, 1, 3, 0));
  lia_assign (results[2], sos_minor3 (i, j, k, 1, 2, 0));

  /* Compute Det(Gamma) */

  upfor (iota, 0, 2)
    {
      lia_push (results[iota]);
      lia_square();
    }
  lia_plus();
  lia_plus();

  /* Compute Det(Lambda) */

  lia_push (sos_minor4 (i, j, k, p, 2, 3, 4, 0));
  lia_push (results[0]);
  /* lia_push (sos_minor3 (i, j, k, 2, 3, 0)); */
  lia_times ();
  ;
  lia_push (sos_minor4 (i, j, k, p, 1, 3, 4, 0));
  lia_push (results[1]);
  /* lia_push (sos_minor3 (i, j, k, 1, 3, 0)); */
  lia_times ();
  ;
  lia_push (sos_minor4 (i, j, k, p, 1, 2, 4, 0));
  lia_push (results[2]);
  /* lia_push (sos_minor3 (i, j, k, 1, 2, 0)); */
  lia_times ();
  ;
  lia_plus ();
  lia_plus ();
  ;
  lia_push (lia_const (2));
  lia_push (sos_minor4 (i, j, k, p, 1, 2, 3, 0));
  lia_push (sos_minor3 (i, j, k, 1, 2, 3));
  lia_times ();
  lia_times ();
  ;
  lia_minus ();

  /* Multiply Det(Gamma) * Det(Lambda) */

  lia_times ();

#ifdef __MDEBUG__
  print ("Triangle (%d, %d, %d) is ", i, j, k);
#endif

  if (math_test)
      fprint (check_file, "testeq [attached2[{%d, %d, %d, %d}], ", i, j, k, p);


  switch (lia_sign (lia_popf ()))       /* The sign is reversed because */
    {                                   /* a negative sign is omitted   */
      case  1: 
        { 
#ifdef __MDEBUG__
          print ("attached to %d.\n", p);
#endif

          if (math_test)
            fprint (check_file, "True]\n");

          return (TRUE);
        }
      case  0: 
        if (math_test)
          fprint (check_file, "True]\n");

#ifdef __DEBUG__
        print ("Degeneracy in attachment test!\n");
#endif
        return (TRUE);
      case -1: 
        { 
#ifdef __MDEBUG__
          print ("not attached to %d.\n", p);
#endif

          if (math_test)
            fprint (check_file, "False]\n");
          return (FALSE);
        }
      default:
       basic_error ("wrong sign");
       return (TRUE);  /* just to silent up lint */
    }
}

/*--------------------------------------------------------------------------*/

/*
 * alf_w_hidden1 (i, j, p) -- Returns TRUE if edge ij is 
 *      attached with respect to p, i.e. if the power distance from p to 
 *      y_T is less than 0, where T = {i,j}.
 */
int alf_w_hidden1 (i, j, p)
     /*i*/ int i, j, p;
{
  int iota;
  int coord[5];
  hidden[1]++;

  /* First, determine proper ordering for the coordinates. */

  coord[4] = 4;
  if (not lia_eq (sos_lia (i, 1), sos_lia (j, 1)))
    {
      coord[1] = 1; coord[2] = 2; coord[3] = 3;
    }
  else if (not lia_eq (sos_lia (i, 2), sos_lia (j, 2)))
    {
      coord[1] = 2; coord[2] = 3; coord[3] = 1;
/*      print ("size1: using x1- and x3- axes\n"); */
    }
  else if (not lia_eq (sos_lia (i, 3), sos_lia (j, 3)))
    {
      coord[1] = 3; coord[2] = 1; coord[3] = 2;
/*      print ("size1: using x1- and x2- axes\n"); */
   }
  else
    basic_error ("is_attached_edge: identical coordinates %d and %d.\n", 
                 i, j);

  /* Pre-compute the determinants that occur more than once. */

  upfor (iota, 0, 2)
    lia_assign (results[iota][0], sos_minor2(i, j, coord[iota+1], 0));

  lia_assign (results[3][0], sos_minor2(i, j, coord[1], coord[2]));
  lia_assign (results[4][0], sos_minor2(i, j, coord[1], coord[3]));

  /* Compute Det(Gamma) */

  Assert (lia_stack_empty ());

  lia_push (results[0][0]);                     /* minor2(i,j,1,0) */
  lia_pushtop ();
  lia_square ();
  lia_push (results[1][0]);                     /* minor2(i,j,2,0) */
  lia_square ();
  lia_push (results[2][0]);                     /* minor2(i,j,3,0) */
  lia_square ();
  lia_plus ();
  lia_plus ();
  /* lia_push (lia_const (-2));
     lia_times ();              */
  lia_times ();

  /* Compute Det(Lambda) */

  upfor (iota, 1, 3)
    {
      lia_push (results[iota-1][0]);                /* minor2(i,j,iota,0) */
      lia_push (sos_minor3 (i, j, p, coord[iota], coord[4], 0));
      lia_times();
    }
  lia_plus ();
  lia_plus ();
  lia_push (results[4][0]);                    /* sos_minor2 (i, j, 1, 3) */
  lia_push (sos_minor3 (i, j, p, coord[1], coord[3], 0));
  lia_times ();
  lia_push (results[3][0]);                    /* sos_minor2 (i, j, 1, 2) */
  lia_push (sos_minor3 (i, j, p, coord[1], coord[2], 0));
  lia_times ();
  lia_plus();
  lia_push (lia_const (2));
  lia_times ();
  lia_minus ();

  lia_push (results[0][0]);                    /* sos_minor2 (i, j, 1, 0) */
  lia_times ();

  lia_push (results[3][0]);                    /* sos_minor2 (i, j, 1, 2) */
  lia_push (results[2][0]);                    /* sos_minor2 (i, j, 3, 0) */
  lia_times ();
  lia_push (results[4][0]);                    /* sos_minor2 (i, j, 1, 3) */
  lia_push (results[1][0]);                    /* sos_minor2 (i, j, 2, 0) */
  lia_times ();
  lia_minus ();
  lia_push (sos_minor3 (i, j, p, coord[2], coord[3], 0));
  lia_push (lia_const (2));
  lia_times ();
  lia_times ();

  lia_plus ();

  /* Multiply Det(Gamma) * Det(Lambda) */

  lia_times ();

#ifdef __MDEBUG__
  print ("Edge (%d, %d) is ", i, j);
#endif

  if (math_test)
      fprint (check_file, "testeq [attached1[{%d, %d, %d}], ", i, j, p);

  switch (lia_sign (lia_popf ()))       /* The sign is reversed because */
    {                                   /* a negative sign is omitted   */
      case  1: 
        { 
#ifdef __MDEBUG__
          print ("attached to %d.\n", p);
#endif
          if (math_test)
            fprint (check_file, "True]\n");
          return (TRUE);
        }
      case  0: 
        {
          if (math_test)
            fprint (check_file, "True]\n");
          
#ifdef __DEBUG__
          print ("Degeneracy in attachment test!\n");
#endif
          return (TRUE);
        }
      case -1: 
        { 
#ifdef __MDEBUG__
          print ("not attached to %d.\n", p);
#endif
          if (math_test)
            fprint (check_file, "False]\n");
          return (FALSE);
        }
      default:
       basic_error ("wrong sign");
       return (TRUE);  /* just to silent up lint */
    }
}

/*--------------------------------------------------------------------------*/

/*
 * alf_w_hidden0 (i, p) -- Returns TRUE if vertex i is 
 *      attached with respect to p, i.e. if the power distance from p to 
 *      y_T is less than 0, where T = {i}.
 */
int alf_w_hidden0 (i, p)
     /*i*/ int i, p;
{
  int a;
  Assert (lia_stack_empty ());
  hidden[0]++;

  upfor (a, 1, 3) {
    lia_push (sos_lia (i, a));
    lia_pushtop ();
    lia_push (sos_lia (p, a));
    lia_minus ();
    lia_times ();
  }
  lia_plus ();
  lia_plus ();

  lia_push (lia_const (2));
  lia_times ();

  lia_push (sos_lia (p, 4));
  lia_plus ();
  lia_push (sos_lia (i, 4));
  lia_minus ();

#ifdef __MDEBUG__
  print ("Vertex %d is ", i);
#endif

  if (math_test)
      fprint (check_file, "testeq [attached0[{%d, %d}], ", i, p);

  switch (lia_sign (lia_popf ()))
    {                                   
      case  1: 
        { 
#ifdef __MDEBUG__
          print ("not attached to %d.\n", p);
#endif
          if (math_test)
            fprint (check_file, "False]\n");
          return (FALSE);
        }
      case  0: 
        {
          if (math_test)
            fprint (check_file, "True]\n");
          
#ifdef __DEBUG__
          print ("Degeneracy in attachment test!\n");
#endif
          return (TRUE);
        }
      case -1: 
        {
#ifdef __MDEBUG__
          print ("attached to %d.\n", p);
#endif
          if (math_test)
            fprint (check_file, "True]\n");
          
          return (TRUE);
        }
      default:
       basic_error ("wrong sign");
       return (TRUE);  /* just to silent up lint */
    }
}

/*--------------------------------------------------------------------------*/

#ifdef __MDEBUG__

static Alf_float value (a, b)
     Lia_ptr a, b;
     /* Returns real value of a/b. */
{
  double return_val;

  if (lia_sign (b) == 0)
    { /* special ratio a/0 */
      return_val = If ((lia_sign (a) == 0), 0.0, 
                       If ((lia_sign(a) < 0), -ALF_INFINITY, ALF_INFINITY));
    }
  else
    { 
      return_val = lia_real (a) / lia_real (b);

      if (Abs (return_val) > ALF_MAXIMUM)
        return_val = Sign (return_val) * ALF_INFINITY;
      if (Abs (return_val) < ALF_MINIMUM)
        return_val = ALF_ZERO;
    }
#ifdef __DEBUG__
  if (return_val == -ALF_INFINITY)
    print ("*** value: returning -ALF_INFINITY = %lf", -ALF_INFINITY);
#endif
  return ((Alf_float) return_val);
}

#endif
