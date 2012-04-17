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
/* sig/combinatorics.c --- Standard combinatoric signatures. */

/*--------------------------------------------------------------------------*/

#ifndef lint
 static char author[] = "Ernst Mucke";
#define _no_status_
#endif

/*--------------------------------------------------------------------------*/

#include "alf.h"

/*--------------------------------------------------------------------------*/

static int is_okay();

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_t_info ()
{
  static Sig_info sm = {0};
  sm.id   = 6;
  sm.name = "#tetrahedra";
  sm.text = "Number of tetrahedra in alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_t;
  return (&sm);
}

Sig_int* sig_num_t ()
{
  Sig_int *s = sig_int_open (sig_num_t_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_TETRA)
              if (alf_ml_r_type () == ALF_RHO)
                { /* tetra, rho: none --> interior */
                  num ++;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  s->min_value = s->value[1];
  s->max_value = s->value[ranks - 1];
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_f_boundary_info ()
{
  static Sig_info sm = {0};
  sm.id   = 21;
  sm.name = "#tri";
  sm.text = "Number of boundary triangles of alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_f_boundary;
  return (&sm);
}

Sig_int* sig_num_f_boundary ()
{
  Sig_int *s = sig_int_open (sig_num_f_boundary_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_TRIANGLE)
              switch (alf_ml_r_type ())
                {
                 case ALF_RHO:
                  /* triangle, rho: none --> singular */
                  num ++;
                  break;
                 case ALF_MU1:  
                  if (not alf_ml_is_attached ())
                    { /* unattached triangle, mu1: singular --> regular */
                      ;
                    }
                  else
                    { /* attached triangle mu1: none --> regular */
                      num ++;
                    }
                  break;
                 case ALF_MU2:
                  /* triangle, mu2: regular --> interior */
                  num --;
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_f_regular_info ()
{
  static Sig_info sm = {0};
  sm.id   = 8;
  sm.name = "#tri_rg";
  sm.text = "Number of regular triangles of alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_f_regular;
  return (&sm);
}

Sig_int* sig_num_f_regular ()
{
  Sig_int *s = sig_int_open (sig_num_f_regular_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_TRIANGLE)
              switch (alf_ml_r_type ())
                {
                 case ALF_MU1:
                  /* triangle, mu1: none/singular --> regular */
                  num ++;
                  break;
                 case ALF_MU2:
                  /* triangle, mu2: regular --> interior */
                  num --;
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_f_singular_info ()
{
  static Sig_info sm = {0};
  sm.id   = 7;
  sm.name = "#tri_sg";
  sm.text = "Number of singular triangles of alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_f_singular;
  return (&sm);
}

Sig_int* sig_num_f_singular ()
{
  Sig_int *s = sig_int_open (sig_num_f_singular_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_TRIANGLE)
              switch (alf_ml_r_type ())
                {
                 case ALF_RHO:
                  /* triangle, rho: none --> singular */
                  num ++;
                  break;
                 case ALF_MU1:
                  if (not alf_ml_is_attached ())
                    { /* unattached triangle, mu1: singular --> regular */
                      num --;
                    }
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_f_interior_info ()
{
  static Sig_info sm = {0};
  sm.id   = 9;
  sm.name = "#tri_in";
  sm.text = "Number of interior triangles in alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_f_interior;
  return (&sm);
}

Sig_int* sig_num_f_interior ()
{
  Sig_int *s = sig_int_open (sig_num_f_interior_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_TRIANGLE)
              switch (alf_ml_r_type ())
                {
                 case ALF_MU2:
                  /* triangle, mu2: regular --> interior */
                  num ++;
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_e_boundary_info ()
{
  static Sig_info sm = {0};
  sm.id   = 22;
  sm.name = "#edges";
  sm.text = "Number of boundary edges of alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_e_boundary;
  return (&sm);
}

Sig_int* sig_num_e_boundary ()
{
  Sig_int *s = sig_int_open (sig_num_e_boundary_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_EDGE)
              switch (alf_ml_r_type ())
                {
                 case ALF_RHO:
                  /* edge, rho: none -> singular */
                  num ++;
                  break;
                 case ALF_MU1:
                  if (not alf_ml_is_attached ())
                    { /* unattached edge, mu1: singular --> regular */
                      ;
                    }
                  else
                    { /* attached edge, mu1: none --> regular */
                      num ++;
                    }
                  break;
                 case ALF_MU2:
                  /* edge, mu2: regular --> interior */
                  num --;
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_e_regular_info ()
{
  static Sig_info sm = {0};
  sm.id   = 11;
  sm.name = "#edges_rg";
  sm.text = "Number of regular edges of alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_e_regular;
  return (&sm);
}

Sig_int* sig_num_e_regular ()
{
  Sig_int *s = sig_int_open (sig_num_e_regular_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_EDGE)
              switch (alf_ml_r_type ())
                {
                 case ALF_MU1:
                  /* edge, mu1: none/singular --> regular */
                  num ++;
                  break;
                 case ALF_MU2:
                  /* edge, mu2: regular --> interior */
                  num --;
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_e_singular_info ()
{
  static Sig_info sm = {0};
  sm.id   = 10;
  sm.name = "#edges_sg";
  sm.text = "Number of singular edges of alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_e_singular;
  return (&sm);
}

Sig_int* sig_num_e_singular ()
{
  Sig_int *s = sig_int_open (sig_num_e_singular_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_EDGE)
              switch (alf_ml_r_type ())
                {
                 case ALF_RHO:
                  /* edge, rho: none -> singular */
                  num ++;
                  break;
                 case ALF_MU1:
                  if (not alf_ml_is_attached ())
                    { /* unattached edge, mu1: singular --> regular */
                      num --;
                    }
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_e_interior_info ()
{
  static Sig_info sm = {0};
  sm.id   = 12;
  sm.name = "#edges_in";
  sm.text = "Number of interior edges in alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_e_interior;
  return (&sm);
}

Sig_int* sig_num_e_interior ()
{
  Sig_int *s = sig_int_open (sig_num_e_interior_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_EDGE)
              switch (alf_ml_r_type ())
                {
                 case ALF_MU2:
                  /* edge, mu2: regular --> interior */
                  num ++;
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_v_boundary_info ()
{
  static Sig_info sm = {0};
  sm.id   = 23;
  sm.name = "#vert";
  sm.text = "Number of boundary vertices of alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_v_boundary;
  return (&sm);
}

Sig_int* sig_num_v_boundary ()
{
  Sig_int *s = sig_int_open (sig_num_v_boundary_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_VERTEX)
              switch (alf_ml_r_type ())
                {
                 case ALF_RHO:
                  /* vertex, rho1: none --> singular */
                  num ++;
                  break;
                 case ALF_MU1:
                  if (not alf_ml_is_attached ())
                    { /* unattached vertex, mu1 : singular --> regular */
                      ;
                    }
                  else
                    { /* attached vertex: none --> regular */
                      num ++;
                    }
                  break;
                 case ALF_MU2:
                  /* vertex, mu2: regular --> interior */
                  num --;
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_v_regular_info ()
{
  static Sig_info sm = {0};
  sm.id   = 14;
  sm.name = "#vert_rg";
  sm.text = "Number of regular vertices of alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_v_regular;
  return (&sm);
}

Sig_int* sig_num_v_regular ()
{
  Sig_int *s = sig_int_open (sig_num_v_regular_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_VERTEX)
              switch (alf_ml_r_type ())
                {
                 case ALF_MU1:
                  /* vertex, mu1: none/singular --> regular */
                  num ++;
                  break;
                 case ALF_MU2:
                  /* vertex, mu2: regular --> interior */
                  num --;
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_v_singular_info ()
{
  static Sig_info sm = {0};
  sm.id   = 13;
  sm.name = "#vert_sg";
  sm.text = "Number of singular vertices of alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_v_singular;
  return (&sm);
}

Sig_int* sig_num_v_singular ()
{
  Sig_int *s = sig_int_open (sig_num_v_singular_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_VERTEX)
              switch (alf_ml_r_type ())
                {
                 case ALF_RHO:
                  /* vertex, rho1: none --> singular */
                  num ++;
                  break;
                 case ALF_MU1:
                  if (not alf_ml_is_attached ())
                    { /* unattached vertex, mu1: singular --> regular */
                      num --;
                    }
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_num_v_interior_info ()
{
  static Sig_info sm = {0};
  sm.id   = 15;
  sm.name = "#vert_in";
  sm.text = "Number of interior vertices in alpha complex.";
  sm.type = SIG_INT;
  sm.int_func = sig_num_v_interior;
  return (&sm);
}

Sig_int* sig_num_v_interior ()
{
  Sig_int *s = sig_int_open (sig_num_v_interior_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_VERTEX)
              switch (alf_ml_r_type ())
                {
                 case ALF_MU2:
                  /* vertex, mu2: --> interior */
                  num ++;
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = num;
    }
  sig_int_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

static int is_okay (s)
     Sig_int *s;
     /* Returns TRUE iff combinatoric signature s is (well... looks) okay. */
{
  int i;
  int num, min = MAXINT, max = -MAXINT;
  upfor (i, 1, s->high - 1)
    {
      num = s->value[i];
      min = Min (min, num);
      max = Max (max, num);
      if (num < 0.0)
        {
          print ("   negative signature value [%d] == %d\n", i, num);
          return (FALSE);
        }
      else if (num < s->min_value)
        {
          print ("   value [%d] == %d < min == %d\n", i, num, s->min_value);
          return (FALSE);
        }
      else if (num > s->max_value)
        {
          print ("   value [%d] == %d > max == %d\n", i, num, s->max_value);
          return (FALSE);
        }
    }
  if ((min != s->min_value) or (max != s->max_value))
    {
      print ("   sig's (min) max == (%d) %d\n", s->min_value, s->max_value);
      print ("                   != (%d) %d\n", min, max);
      return (FALSE);
    }
  if (   (s->value[0] > s->value[1])
      or (s->value[s->high - 1] != s->value[s->high]))
    return (FALSE);
  return (TRUE);
}
