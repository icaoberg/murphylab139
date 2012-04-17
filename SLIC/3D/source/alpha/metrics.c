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
/* sig/metrics.c --- Standard metric signatures. */

/*--------------------------------------------------------------------------*/

#ifndef lint
 static char author[] = "Ernst Mucke";
#define _no_status_
#endif

/*--------------------------------------------------------------------------*/

#include "alf.h"

/*--------------------------------------------------------------------------*/

static Alf_float volume(), area(), length(), sub();
static int is_okay();
static int corrections_0, corrections_1;

extern char *sig_prefix_string;  /* silently imported from sig.c */

/*--------------------------------------------------------------------------*/

Sig_info* sig_volume_info ()
{
  static Sig_info sm;
  sm.id   = 1;
  sm.name = "volume";
  sm.text = "Volume of alpha shape.";
  sm.type = SIG_FLOAT;
  sm.float_func = sig_volume;
  return (&sm);
}

Sig_float* sig_volume ()
{
  Sig_float *s = sig_float_open (sig_volume_info ());
  int rank, ranks = alf_ml_ranks ();
  Alf_float value = 0.0;
  upfor (rank, 1, ranks)
    {
      if (alf_ml_sublist (rank))
        do
          {
            if (alf_ml_f_type () == ALF_TETRA)
              if (alf_ml_r_type () == ALF_RHO)
                { /* tetra, rho: {} --> interior */
                  value += volume (alf_ml_face ());
                }
          } while (alf_ml_next ());
      s->value[rank] = value;
    }
  sig_float_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_area_f_boundary_info ()
{
  static Sig_info sm;
  sm.id   = 2;
  sm.name = "area";
  sm.text = "Boundary area of alpha shape.";
  sm.type = SIG_FLOAT;
  sm.float_func = sig_area_f_boundary;
  return (&sm);
}

Sig_float* sig_area_f_boundary ()
{
  Sig_float *s = sig_float_open (sig_area_f_boundary_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  Alf_float value_m = 0.0;
  Alf_float value_p = 0.0;
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
                  value_p += 2.0 * area (alf_ml_face ());
                  break;
                 case ALF_MU1:  
                  if (not alf_ml_is_attached ())
                    { /* unattached triangle, mu1: singular --> regular */
                      value_m += area (alf_ml_face ());
                      /* NOTE: boundary area = regular + 2 * singular */
                    }
                  else
                    { /* attached triangle mu1:  {} --> regular */
                      num ++;
                      value_p += area (alf_ml_face ());
                    }
                  break;
                 case ALF_MU2:
                  /* triangle, mu2: regular --> interior */
                  num --;
                  value_m += area (alf_ml_face ());
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = sub (value_p, value_m, num);
    }
  sig_float_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_area_f_regular_info ()
{
  static Sig_info sm;
  sm.id   = 18;
  sm.name = "area_rg";
  sm.text = "Regular boundary area of alpha shape.";
  sm.type = SIG_FLOAT;
  sm.float_func = sig_area_f_regular;
  return (&sm);
}

Sig_float* sig_area_f_regular ()
{
  Sig_float *s = sig_float_open (sig_area_f_regular_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  Alf_float value_m = 0.0;
  Alf_float value_p = 0.0;
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
                  value_p += area (alf_ml_face ());
                  break;
                 case ALF_MU2:
                  /* triangle, mu2: regular --> interior */
                  num --;
                  value_m += area (alf_ml_face ());
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = sub (value_p, value_m, num);
    }
  sig_float_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_area_f_singular_info ()
{
  static Sig_info sm;
  sm.id   = 19;
  sm.name = "area_sg";
  sm.text = "Singular boundary area of alpha shape.";
  sm.type = SIG_FLOAT;
  sm.float_func = sig_area_f_singular;
  return (&sm);
}

Sig_float* sig_area_f_singular ()
{
  Sig_float *s = sig_float_open (sig_area_f_singular_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  Alf_float value_m = 0.0;
  Alf_float value_p = 0.0;
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
                  value_p += area (alf_ml_face ());
                  break;
                 case ALF_MU1:
                  if (not alf_ml_is_attached ())
                    { /* unattached, triangle, mu1: singular --> regular */
                      num --;
                      value_m += area (alf_ml_face ());
                    }
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = sub (value_p, value_m, num);
    }
  sig_float_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_length_e_singular_info ()
{
  static Sig_info sm;
  sm.id   = 20;
  sm.name = "length_sg";
  sm.text = "Length of singular edges of alpha shape.";
  sm.type = SIG_FLOAT;
  sm.float_func = sig_length_e_singular;
  return (&sm);
}

Sig_float* sig_length_e_singular ()
{
  Sig_float *s = sig_float_open (sig_length_e_singular_info ());
  int rank, ranks = alf_ml_ranks ();
  int num = 0;
  Alf_float value_m = 0.0;
  Alf_float value_p = 0.0;
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
                  value_p += length (alf_ml_face ());
                  break;
                 case ALF_MU1:
                  if (not alf_ml_is_attached ())
                    { /* unattached edge, mu1: singular --> regular */
                      num --;
                      value_m += length (alf_ml_face ());
                    }
                  break;
                }
          } while (alf_ml_next ());
      s->value[rank] = sub (value_p, value_m, num);
    }
  sig_float_range (s);
  Assert (is_okay (s));
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_w_info ()
{
  static Sig_info sm;
  sm.id   = 24;
  sm.name = "w";
  sm.text = STRDUP (basic_cb_frmt ("A weighted combination of %s, %s, and %s.",
                                   sig_volume_info ()->name,
                                   sig_area_f_boundary_info ()->name,
                                   sig_length_e_singular_info ()->name));
  sm.type = SIG_FLOAT;
  sm.float_func = sig_w;
  return (&sm);
}

Sig_float* sig_w ()
{
  Sig_float *s = sig_float_open (sig_w_info ());
  int ranks = alf_ml_ranks ();
  char *save = sig_prefix_string;
  sig_prefix_string = STRDUP (basic_cb_frmt ("  %s", sig_prefix_string));
  { /* w = 4 area^{1/2} - ... */
    Sig_float *x = sig_area_f_boundary ();
    Alf_float *v = x->value;
    Alf_float *w = s->value;
    int r = 0;
    while (r < ranks)
      {
        v ++;
        w ++;
        r ++;
        *w = (Alf_float) 4.0 * sqrt ((double) *v);
      }
    sig_FREE (x);
  }
  { /* w = 4 area^{1/2} - 8 volume^{1/3} - ... */
    Sig_float *x = sig_volume ();
    Alf_float *v = x->value;
    Alf_float *w = s->value;
    int r = 0;
    while (r < ranks)
      {
        v ++;
        w ++;
        r ++;
        *w = *w - (Alf_float) 8.0 * cbrt ((double) *v);
      }
    sig_FREE (x);
  }
  { /* w = Max {0, 4 area^{1/2} - 8 volume^{1/3} - 2.0 length_sg}*/
    Sig_float *x = sig_length_e_singular ();
    Alf_float *v = x->value;
    Alf_float *w = s->value;
    int r = 0;
    while (r < ranks)
      {
        v ++;
        w ++;
        r ++;
        *w = *w - (Alf_float) 2.0 * (*v);
        if (*w < 0.0)
          *w = 0.0;
      }
    sig_FREE (x);
  }
  sig_float_range (s);
  Assert (is_okay (s));
  FREE (sig_prefix_string);
  sig_prefix_string = save;
  return (s);
}

/*--------------------------------------------------------------------------*/

static Alf_float volume (ef)
     int ef;
{
  int a, b, c;
  trist_triangle (ef, &a, &b, &c);
  return (alf_volume (a, b, c, Dest (Enext (Fnext (ef)))));
}

/*--------------------------------------------------------------------------*/

static Alf_float area (ef)
     int ef;
{
  int a, b, c;
  trist_triangle (ef, &a, &b, &c);
  return (alf_triangle_area (a, b, c));
}

/*--------------------------------------------------------------------------*/

static Alf_float length (ef)
     int ef;
{
  return (alf_edge_length (Org (ef), Dest (ef)));
}

/*--------------------------------------------------------------------------*/

static Alf_float sub (metric_p, metric_m, num)
     Alf_float metric_p, metric_m;
     int num;
     /* Returns: metric_p - metric_m; but possibly with "zero correction",
        in case the combinatorical number num is zero. */
{
  Alf_float m = metric_p - metric_m;
  Assert ((metric_p >= 0.0) and (metric_m >= 0.0) and (num >= 0));
  if (num == 0)
    {
      if (m < 0.0)
        corrections_0 ++;  /* count only the negative ones! */
      m = 0.0;
    }
  if (m < 0.0)
    {
      corrections_1 ++;
      m = 0.0;  /* even this still happens because of numerical error! */
    }
  Assert (m >= 0.0);
  return (m);
}

/*--------------------------------------------------------------------------*/

void sig_get_metrics_corrections (c0, c1)
     int *c0, *c1;
{
  *c0 = corrections_0;
  *c1 = corrections_1;
}

/*--------------------------------------------------------------------------*/

static int is_okay (s)
     Sig_float *s;
     /* Returns TRUE iff metric signature s looks okay. */
{
  int i;
  Alf_float val, min = MAXFLOAT, max = -MAXFLOAT;
  upfor (i, 1, s->high - 1)
    {
      val = s->value[i];
      min = Min (min, val);
      max = Max (max, val);
      if (val < 0.0)
        {
          print ("   negative signature value [%d] == %f\n", i, val);
          return (FALSE);
        };
      if (val < s->min_value)
        {
          print ("   value [%d] == %f < min == %f", i, val, s->min_value);
          return (FALSE);
        }
      else if (val > s->max_value)
        {
          print ("   value [%d] == %f > max == %f", i, val, s->max_value);
          return (FALSE);
        }
    }
  if ((min != s->min_value) or (max != s->max_value))
    {
      print ("   sig's (min) max == (%f) %f\n", s->min_value, s->max_value);
      print ("                   != (%f) %f\n", min, max);
      return (FALSE);
    }
  if (   (s->value[0] > s->value[1])
      or (s->value[s->high - 1] != s->value[s->high]))
    return (FALSE);
  return (TRUE);
}
