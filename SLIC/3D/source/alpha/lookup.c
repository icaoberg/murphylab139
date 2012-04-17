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
/* mkalf/lookup.c --- Alf Library: interface to spectrum & rank tables. */

/*--------------------------------------------------------------------------*/

#include "mkalf.h"

#define SQRT  alf_sqrt

/*--------------------------------------------------------------------------*/

Alf_float alf_sqrt (v)
     Alf_float v;
     /* The "alpha shape" squareroot. */
{
  int negative = (v < 0);
  Alf_float abs_v = If (negative, -v, v);
  if (abs_v >= ALF_INFINITY)
    return (v);
  else if (negative)
    return (-sqrt (abs_v));
  else
    return (sqrt (abs_v));
}

/*--------------------------------------------------------------------------*/

Alf_float alf_threshold (r)
     int r;
     /* Returns alpha = SQRT(spectrum[r]) (= threshold), for 1 <= r <= ranks.
        Returns -ALF_INFINITY for r < 1, and ALF_INFINITY for r >= ranks.
        NOTE: alpha == alf_threshold (r) iff r == alf_rank (alpha) */
{
  return (SQRT (alf_threshold2 (r)));
}

/*--------------------------------------------------------------------------*/

Alf_float alf_threshold2 (r)
     int r;
     /* Returns v = spectrum[r] (ie, square of threshold), for 1 <= r <= ranks.
        Returns -ALF_INFINITY for r < 1, and ALF_INFINITY for r >= ranks.
        (In the unweighted case: v == alpha^2 for alpha shape of rank r.)
        NOTE: v == alf_threshold2 (r) iff r == alf_rank2 (v) */
{
  r = Max (0, r);
  r = Min (r, alf_context->ranks);
  return (alf_context->spectrum[r]);
}

/*--------------------------------------------------------------------------*/

void alf_interval (r, a, b)
     int r;
     Alf_float *a, *b;  /* output */
     /* Returns interval [&a, &b] = [SQRT(spectrum[r]), SQRT(spectrum[r+1])]
        of threshold values corresponding to rank r.
        Returns [SQRT(spectrum[ranks-1]), SQRT(spectrum[ranks])], if r >= ranks
        Returns [SQRT(spectrum[0]), SQRT(spectrum[1])], if r < 1. */
{
  r = Max (0, r);
  r = Min (r, alf_context->ranks - 1);
  alf_interval2 (r, a, b);
  *a = SQRT (*a);
  *b = SQRT (*b);
}

/*--------------------------------------------------------------------------*/

void alf_interval2 (r, a, b)
     int r;
     Alf_float *a, *b;  /* output */
     /* Returns interval [&a, &b] = [spectrum[r], spectrum[r+1]] of squares
        of threshold values corresponding to rank r.
        Returns [spectrum[ranks-1], spectrum[ranks]], if r >= ranks.
        Returns [spectrum[0], spectrum[1]], if r < 1. */
{
  r = Max (0, r);
  r = Min (r, alf_context->ranks - 1);
  *a = alf_context->spectrum[r];
  *b = alf_context->spectrum[r + 1];
}

/*--------------------------------------------------------------------------*/

int alf_rank (alpha)
     Alf_float alpha;
     /* Cf, alf_rank2 (). */
{
  return (alf_rank2 (Sign (alpha) * alpha * alpha));
}

/*--------------------------------------------------------------------------*/

int alf_rank2 (v)
     Alf_float v;
     /* Returns rank r of value v wrt intervals in spectrum[1..ranks].
        The following holds for result r:
        - spectrum[r] <= v < spectrum[r+1], for v in range of spectrum;
        - r == 1, for v < spectrum[1];
        - r == ranks, for v >= spectrum[ranks].
        NOTE: This means that rank 0 (denoting empty set) is unaccessible. */
{
  int left = 1, right = alf_context->ranks, mid;
  if (v <= alf_context->spectrum[left])
    return (left);
  else if (v >= alf_context->spectrum[right])
    return (right);
  else
    {
      loop
        {
          if (right - left <= 1)
            break;
          else
            {
              mid = (left + right) / 2;
              if (v < alf_context->spectrum[mid])
                right = mid;
              else
                left = mid;
            }
        }
      return (left);
    }
}

/*--------------------------------------------------------------------------*/

int alf_is_redundant (i)
     int i;
     /* Returns TRUE (nonzero) iff vertex with index i is redundant/dumped,
        based on rank-table lookup: mu1 == mu2 <== 0 ==> vertex redundant. */
{
  Assert (alf_context and (1 <= i) and (i <= alf_context->dt->n));
  return (    (alf_context->v_rank[i].mu1 == 0)
          and (alf_context->v_rank[i].mu2 == 0));
}

/*--------------------------------------------------------------------------*/

int alf_is_valid (f_type, i)
     int f_type, i;
     /* Returns TRUE iff index i is valid index of face f_type;
        f_type is one of {ALF_VERTEX, ALF_EDGE, ALF_TRIANGLE, ALF_TETRA}
        (which is also assumed in all alf_is_*() functions below!).
        NOTE: redundant/dumped vertices are not cosidered valid! */
{
  Assert (alf_proper (alf_context));
  if (i < 1)
    return (FALSE);
  else
    {
      switch (f_type)
        {
         case ALF_VERTEX:
          return (    (i <= alf_context->dt->n)
                  and (  alf_context->v_rank[i].mu1
                       + alf_context->v_rank[i].mu2 != 0));
         case ALF_EDGE:
          return (    (i <= alf_context->e_hash_m)
                  and (  alf_context->e_rank[i].mu1
                       + alf_context->e_rank[i].mu2 != 0));
         case ALF_TRIANGLE:
          return (i <= alf_context->dt->num.f);
         case ALF_TETRA:
          return (    (i <= alf_context->t_hash_m)
                  and (alf_context->t_rank[i].rho != 0));
        }
      return (FALSE);
    }
}

/*--------------------------------------------------------------------------*/

int alf_is_singular (f_type, rank, i)
          int f_type, i, rank;
          /* Looks up rank tables, and returns TRUE iff face with index i
             is singular wrt given rank; dimension of face is given by f_type.
             (Obviously, a tetrahedron can never be singular.) */
{
  Assert (alf_is_valid (f_type, i));
  switch (f_type)
    {
     case ALF_VERTEX:   return (    (alf_context->v_rank[i].rho >  0)
                                and (alf_context->v_rank[i].rho <= rank)
                                and (alf_context->v_rank[i].mu1 >  rank));
     case ALF_EDGE:     return (    (alf_context->e_rank[i].rho >  0)
                                and (alf_context->e_rank[i].rho <= rank)
                                and (alf_context->e_rank[i].mu1 >  rank));
     case ALF_TRIANGLE: return (    (alf_context->f_rank[i].rho >  0)
                                and (alf_context->f_rank[i].rho <= rank)
                                and (alf_context->f_rank[i].mu1 >  rank));
    }
  return (FALSE);
}

/*--------------------------------------------------------------------------*/

int alf_is_interior (f_type, rank, i)
     int f_type, i, rank;
     /* Looks up rank tables, and Returns TRUE iff face with index i
        is interior wrt given rank; dimension of face is given by f_type. */
{
  Assert (alf_is_valid (f_type, i));
  if (alf_is_on_hull (f_type, i))
    return (FALSE);
  else
    {
      switch (f_type)
        {
         case ALF_VERTEX:   return (alf_context->v_rank[i].mu2 <= rank);
         case ALF_EDGE:     return (alf_context->e_rank[i].mu2 <= rank);
         case ALF_TRIANGLE: return (alf_context->f_rank[i].mu2 <= rank);
         case ALF_TETRA:    return (alf_context->t_rank[i].rho <= rank);
        }
    }
  return (FALSE);
}

/*--------------------------------------------------------------------------*/

int alf_is_in_complex (f_type, rank, i)
  int f_type, i, rank;
  /* Looks up rank tables, and Returns TRUE iff face with index i
     is in alpha complex (ie, interior, regular, or singular) wrt given rank;
     dimension of face is given by its f_type. */
{
  Assert (alf_is_valid (f_type, i));
  switch (f_type)
    {
     case ALF_VERTEX:
      {
        if (alf_context->v_rank[i].rho)
          return (alf_context->v_rank[i].rho <= rank);
        else
          return (alf_context->v_rank[i].mu1 <= rank);
      }
     case ALF_EDGE:
      {
        if (alf_context->e_rank[i].rho)
          return (alf_context->e_rank[i].rho <= rank);
        else
          return (alf_context->e_rank[i].mu1 <= rank);
      }
     case ALF_TRIANGLE:
      {
        if (alf_context->f_rank[i].rho)
          return (alf_context->f_rank[i].rho <= rank);
        else
          return (alf_context->f_rank[i].mu1 <= rank);
      }
     case ALF_TETRA:
      {
        return (alf_context->t_rank[i].rho <= rank);
      }
    }
  return (FALSE);
}

/*--------------------------------------------------------------------------*/

int alf_is_on_hull (f_type, i)
     int f_type, i;
     /* Returns TRUE iff face with index i and dimension f_type lies on the
        (convex) hull; decison is based on rank tables.
        (In this context: tetrahedra never lie on the hull!) */
{
  Assert (alf_is_valid (f_type, i));
  switch (f_type)
    {
     case ALF_VERTEX:   return (alf_context->v_rank[i].mu2 == 0);
     case ALF_EDGE:     return (alf_context->e_rank[i].mu2 == 0);
     case ALF_TRIANGLE: return (alf_context->f_rank[i].mu2 == 0);
     case ALF_TETRA:    return (FALSE);
    }
  return (FALSE);
}
