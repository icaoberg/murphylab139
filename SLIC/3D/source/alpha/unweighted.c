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
/* mkalf/unweighted.c --- Non-SoS primitives for (unweighted) alpha shapes. */

/* NOTE: - Computations are done with Lia but without SoS.
           There are no degenerate tetrahedra, triangles, edges!
         - Some dummy procedures/functions were added to keep things
           compatible to mkalf/weighted.c for weighted alpha shapes. */

#ifndef lint
 static char author[] = "Ernst Mucke";
#endif

/*--------------------------------------------------------------------------*/

#include "mkalf.h"

/*--------------------------------------------------------------------------*/

/* counters */

static int rho[4] = { 0 };
static int hidden[4] = { 0 };

/*--------------------------------------------------------------------------*/

void alf_uw_begin ()
     /* Initializes this module.  Dummy procedure. */
{
  int i;
  upfor (i, 0, 3)
    rho[i] = hidden[i] = 0;
      
}

/*--------------------------------------------------------------------------*/

void alf_uw_end ()
     /* Closes this module.  Dummy procedure. */
{
}

/*--------------------------------------------------------------------------*/

void alf_uw_calls (rho_cnt, hidden_cnt, dim)
     int **rho_cnt;
     int **hidden_cnt;
     int *dim;
     /* Returns counters of corresponding function calls.
        See mkalf/mkalf.c for sample usage. */
        
{
  *rho_cnt = rho;
  *hidden_cnt = hidden;
  *dim = 3;
}

/*--------------------------------------------------------------------------*/

/*ARGSUSED*/

void alf_uw_rho0 (i, a, b)
     int i;
     Lia_obj a, b;  /* output */
     /* Computes (unweighted) rho value (a/b) of vertex.  Dummy procedure.
        Will always return zero; ie, a/b == 0/1. */
     /* GENERAL ASSUMPTIONS FOR EACH OF THE FOLLOWING rhoI (<INDICES>, a, b)
        FUNCTIONS:
        - <INDICES> in proper range and pairwise different
        - a and b point to Lia objects of sufficient size
        - a != 0 and b > 0 always holds after the call
        - actually: in the unweighted case, a > 0 */
     
{
  rho[0] ++;
  lia_assign (a, lia_const (0));
  lia_assign (b, lia_const (1));
}

/*--------------------------------------------------------------------------*/

void alf_uw_rho1 (i, j, a, b)
     int i, j;
     Lia_obj a, b;  /* output */
     /* Computes (unweighted) rho value (a/b) of edge in R^3;
        ie, radius of circumsphere of edge given by indices i, j. */
{
  rho[1] ++;
  ;
  lia_push (sos_minor2 (i, j, 1, 0));
  lia_ipower (2);
  lia_push (sos_minor2 (i, j, 2, 0));
  lia_ipower (2);
  lia_push (sos_minor2 (i, j, 3, 0));
  lia_ipower (2);
  lia_plus ();
  lia_plus ();
  lia_assign (a, lia_popf ());
  Assert_always (lia_sign (a) > 0);
  ;
  lia_assign (b, lia_const (4));
  Assert_always (lia_sign (b) > 0);
}

/*--------------------------------------------------------------------------*/

void alf_uw_rho2 (i, j, k, a, b)
     int i, j, k;
     Lia_obj a, b;  /* output */
     /* Computes (unweighted) rho value (a/b) of triangle in R^3;
        ie, radius of circumsphere of triangle given by indices i, j, k. */
{
  rho[2] ++;
  ;
  lia_push (sos_minor2 (i, j, 1, 0));
  lia_ipower (2);
  lia_push (sos_minor2 (i, j, 2, 0));
  lia_ipower (2);
  lia_push (sos_minor2 (i, j, 3, 0));
  lia_ipower (2);
  lia_plus ();
  lia_plus ();
  lia_push (sos_minor2 (i, k, 1, 0));
  lia_ipower (2);
  lia_push (sos_minor2 (i, k, 2, 0));
  lia_ipower (2);
  lia_push (sos_minor2 (i, k, 3, 0));
  lia_ipower (2);
  lia_plus ();
  lia_plus ();
  lia_push (sos_minor2 (j, k, 1, 0));
  lia_ipower (2);
  lia_push (sos_minor2 (j, k, 2, 0));
  lia_ipower (2);
  lia_push (sos_minor2 (j, k, 3, 0));
  lia_ipower (2);
  lia_plus ();
  lia_plus ();
  lia_times ();
  lia_times ();
  lia_assign (a, lia_popf ());
  Assert_always (lia_sign (a) > 0);
  ;
  lia_push (lia_const (4));
  lia_push (sos_minor3 (i, j, k, 1, 2, 0));
  lia_ipower (2);
  lia_push (sos_minor3 (i, j, k, 1, 3, 0));
  lia_ipower (2);
  lia_push (sos_minor3 (i, j, k, 2, 3, 0));
  lia_ipower (2);
  lia_plus ();
  lia_plus ();
  lia_times ();
  lia_assign (b, lia_popf ());
  Assert_always (lia_sign (b) > 0);  
}

/*--------------------------------------------------------------------------*/

void alf_uw_rho3 (i, j, k, l, a, b)
     int i, j, k, l;
     Lia_obj a, b;  /* output */
     /* Computes (unweighted) rho value (a/b) of tetrahedron in R^3;
        ie, radius of circumsphere of tetrahedron given by i, j, k, l. */
{
  rho[3] ++;
  ;
  lia_push (lia_const (4));
  lia_push (sos_minor4 (i, j, k, l, 1, 2, 3, 0));
  lia_push (sos_minor4 (i, j, k, l, 1, 2, 3, 4));
  lia_times ();
  lia_times ();
  lia_push (sos_minor4 (i, j, k, l, 1, 2, 4, 0));
  lia_ipower (2);
  lia_push (sos_minor4 (i, j, k, l, 1, 3, 4, 0));
  lia_ipower (2);
  lia_push (sos_minor4 (i, j, k, l, 2, 3, 4, 0));
  lia_ipower (2);
  lia_plus ();
  lia_plus ();
  lia_plus ();
  lia_assign (a, lia_popf ());
  Assert_always (lia_sign (a) > 0);
  ;
  lia_push (lia_const (4));
  lia_push (sos_minor4 (i, j, k, l, 1, 2, 3, 0));
  lia_ipower (2);
  lia_times ();
  lia_assign (b, lia_popf ());
  Assert_always (lia_sign (b) > 0);
}

/*--------------------------------------------------------------------------*/

/*ARGSUSED*/

int alf_uw_hidden0 (i, p)
     int i, p;
     /* Dummy function.  Returns always FALSE. */
{
  hidden[0] ++;
  return (FALSE);
}

/*--------------------------------------------------------------------------*/

int alf_uw_hidden1 (i, j, p)
     int i, j, p;
     /* Returns TRUE (ie, nonzero) iff point p is "hidden" by (or: lies inside)
        the (unweighted) smallest sphere spanned by the edge with endpoints
        i, j.
        In the nondegenerate case, edge ij is "attached" iff there exists
        a adjacent point p that is hidden. Since the triangulation is Delaunay,
        there is always some adjacent point p which is not hidden; but this is
        only true in the nondegenerate case!
        In the degenerate case (ie, if p lies on the sphere), the function 
        returns also TRUE (in fact, it really returns TRUE + TRUE == 2 :). */
{
  int iota;
  hidden[1] ++;
  ;
  upfor (iota, 1, 3)     /* push d^2 (i, j) */
    {
      lia_push (sos_minor2 (i, j, iota, 0));
      lia_ipower (2);
    }
  lia_plus ();
  lia_plus ();
  ;
  upfor (iota, 1, 3)    /* subtract d^2 (midpoint (i, j), p) */
    {
      lia_push (sos_minor2 (i, p, iota, 0));
      lia_push (sos_minor2 (j, p, iota, 0));
      lia_plus ();
      lia_ipower (2);
      lia_minus ();
    }    
  ;
  switch (lia_sign (lia_popf ()))
    {
      case  1: return (TRUE);
      case  0: return (TRUE + TRUE);
      case -1: return (FALSE);
      default:
       Assert_always (FALSE);
       return (FALSE);  /* just to silent up lint */
    }
}

/*--------------------------------------------------------------------------*/

int alf_uw_hidden2 (i, j, k, p)
     int i, j, k, p;
     /* Returns TRUE (ie, nonzero) iff point p is "hidden" by (or: lies inside)
        the (unweighted) smallest sphere spanned by the triangle with endpoints
        i, j, k.
        In the nondegenerate case, triangle ijk is "attached" iff there exists
        a adjacent point p that is hidden. Since the triangulation is Delaunay,
        there is always some adjacent point p which is not hidden; but this is
        only true in the nondegenerate case!
        In the degenerate case (ie, if p lies on the sphere), the function 
        returns also TRUE (in fact, it really returns TRUE + TRUE == 2 :). */
{
  hidden[2] ++;
  ;
  lia_push (sos_minor4 (i, j, k, p, 2, 3, 4, 0));
  lia_push (sos_minor3 (i, j, k, 2, 3, 0));
  lia_times ();
  ;
  lia_push (sos_minor4 (i, j, k, p, 1, 3, 4, 0));
  lia_push (sos_minor3 (i, j, k, 1, 3, 0));
  lia_times ();
  ;
  lia_push (sos_minor4 (i, j, k, p, 1, 2, 4, 0));
  lia_push (sos_minor3 (i, j, k, 1, 2, 0));
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
  ;
  switch (lia_sign (lia_popf ()))
    {
      case  1: return (TRUE);
      case  0: return (TRUE + TRUE);
      case -1: return (FALSE);
      default:
       Assert_always (FALSE);
       return (FALSE);  /* just to silent up lint */
    }
}
