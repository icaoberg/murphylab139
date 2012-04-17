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
/* sig/sig.c --- Core module for alpha-shape signatures. */

/*--------------------------------------------------------------------------*/

#ifndef lint
 static char author[] = "Ernst Mucke";
#define _no_status_
#endif

/*--------------------------------------------------------------------------*/

#include "mkalf.h"

/*--------------------------------------------------------------------------*/

static Alf_vect *coord;
static Alf_float fminor3one(), fminor4one();

char *sig_prefix_string = "Signature";  /* Internal use only! */

/*--------------------------------------------------------------------------*/

void sig_core ()
     /* Initializes the sig.c module. */
{
  coord = alf_get_coords ();
}

/*--------------------------------------------------------------------------*/

Sig_float* sig_float_open (sm)
     Sig_info *sm;
     /* Allocates space for signature and initializes it. */
{
  Sig_float *s;
  print ("%s: \"%s\" ...\n", sig_prefix_string, sm->name);
  s = MALLOC (Sig_float, 1);
  s->high = alf_context->ranks;
  s->min_value =  MAXFLOAT;
  s->max_value = -MAXFLOAT;
  s->value = MALLOC (Alf_float, s->high + 1);
  s->value[0] = 0;
  sig_core ();
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_int* sig_int_open (sm)
     Sig_info *sm;
     /* Allocates space for signature and initializes it. */
{
  Sig_int *s;
  print ("%s: \"%s\" ...\n", sig_prefix_string, sm->name);
  s = MALLOC (Sig_int, 1);
  s->high = alf_context->ranks;
  s->min_value =  MAXINT;
  s->max_value = -MAXINT;
  s->value = MALLOC (int, s->high + 1);
  s->value[0] = 0;
  sig_core ();
  return (s);
}

/*--------------------------------------------------------------------------*/

void sig_int_range (s)
     Sig_int *s;  /* input/output */
     /* computes s->min_value and s->max_value of s->value[1..s->high] */
{
  int i = s->high;
  int *v = &(s->value[i]);
  int min = *v, max = *v;
  while (i > 1)
    {
      i --;
      v --;
      max = Max (max, *v);
      min = Min (min, *v);
    }
  s->min_value = min;
  s->max_value = max;
}

/*--------------------------------------------------------------------------*/

void sig_float_range (s)
     Sig_float *s;  /* input/output */
     /* computes s->min_value and s->max_value of s->value[1..s->high] */
{
  int i = s->high;
  Alf_float *v = &(s->value[i]);
  Alf_float min = *v, max = *v;
  while (i > 1)
    {
      i --;
      v --;
      max = Max (max, *v);
      min = Min (min, *v);
    }
  s->min_value = min;
  s->max_value = max;
}

/*--------------------------------------------------------------------------*/

#if defined (sgi) || defined (__GNUC__)

static Alf_float square (Alf_float x)
{
  return (x * x);
}

#else

static Alf_float square (x)
     Alf_float x;
{
  return (x * x);
}

#endif

/*--------------------------------------------------------------------------*/

#if defined (sgi)
# define sqroot  fsqrt
#else
# define sqroot   sqrt
#endif

/*--------------------------------------------------------------------------*/

Alf_float alf_volume (i, j, k, l)
     int i, j, k, l;
{
  return (fabs (fminor4one (i, j, k, l, ALF_X, ALF_Y, ALF_Z)
                / (Alf_float) 6.0));
}

/*--------------------------------------------------------------------------*/

Alf_float alf_triangle_area (i, j, k)
     int i, j, k;
{
  return (sqroot ((  square (fminor3one (i, j, k, ALF_Y, ALF_Z))
                   + square (fminor3one (i, j, k, ALF_Z, ALF_X))
                   + square (fminor3one (i, j, k, ALF_X, ALF_Y)))
                  / (Alf_float) 4.0));
}

/*--------------------------------------------------------------------------*/

Alf_float alf_edge_length (i, j)
     int i, j;
{
  return (sqroot (  square (coord[i][ALF_X] - coord[j][ALF_X])
                  + square (coord[i][ALF_Y] - coord[j][ALF_Y])
                  + square (coord[i][ALF_Z] - coord[j][ALF_Z])));
}

/*--------------------------------------------------------------------------*/

static Alf_float fminor3one (i, j, k, a, b)
     int i, j, k, a, b;
{
  Alf_float o, p, q, r;
  o = coord[j][a] - coord[i][a];
  p = coord[j][b] - coord[i][b];  /* 2nd line - 1st line */
  q = coord[k][a] - coord[i][a];
  r = coord[k][b] - coord[i][b];  /* 3rd line - 1st line */
  return (o * r - q * p);  
}

/*--------------------------------------------------------------------------*/

static Alf_float fminor4one (i, j, k, l, a, b, c)
     int i, j, k, l, a, b, c;
{
  Alf_float o, p, q, r, s, t, u, v, w;
  o = coord[j][a] - coord[i][a];
  p = coord[j][b] - coord[i][b];
  q = coord[j][c] - coord[i][c];  /* 2nd line - 1st line */
  r = coord[k][a] - coord[i][a];
  s = coord[k][b] - coord[i][b];
  t = coord[k][c] - coord[i][c];  /* 3rd line - 1st line */
  u = coord[l][a] - coord[i][a];
  v = coord[l][b] - coord[i][b];
  w = coord[l][c] - coord[i][c];  /* 4th line - 1st line */
  return (  q * s * u - o * s * w
          + o * t * v - p * t * u
          + p * r * w - q * r * v);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_spectrum_info ()
{
  static Sig_info sm;
  sm.id   = 17;
  sm.name = "spectrum";
  sm.text = "Threshold values of alpha spectrum (radii squared).";
  sm.type = SIG_FLOAT;
  sm.float_func = sig_spectrum;
  return (&sm);
}

Sig_float* sig_spectrum ()
     /* This is the basic signature, the alpha-spectrum of the current Alf.
        NOTE: Obviously, alf_set() needs to be called beforehand. */
{
  int i = 0;
  Sig_float *s = sig_float_open (sig_spectrum_info ());
  upfor (i, 0, s->high)
    s->value[i] = alf_context->spectrum[i];
  s->min_value = s->value[1];
  s->max_value = s->value[s->high - 1];
  return (s);
}

/*--------------------------------------------------------------------------*/

Sig_info* sig_alpha_info ()
{
  static Sig_info sm;
  sm.id   = 16;
  sm.name = "alpha";
  sm.text = "Radius used to define the current shape or complex.";
  sm.type = SIG_FLOAT;
  sm.float_func = sig_alpha;
  return (&sm);
}

Sig_float* sig_alpha ()
     /* This is the basic signature, the alpha-spectrum of the current Alf.
        NOTE: Obviously, alf_set() needs to be called beforehand. */
{
  int i = 0;
  Sig_float *s = sig_float_open (sig_alpha_info ());
  upfor (i, 0, s->high)
    s->value[i] = alf_sqrt (alf_context->spectrum[i]);
  s->min_value = s->value[1];
  s->max_value = s->value[s->high - 1];
  return (s);
}
