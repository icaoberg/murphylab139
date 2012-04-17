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
/* mkalf/alf.c  --- Alf Library: primary interface to Alf. */

/*--------------------------------------------------------------------------*/

#ifndef lint
 static char version[] = "@(#) Alf Library 2.2";
 static char author[] = "Ernst Mucke";
#include "copyright.h"
#endif

/*--------------------------------------------------------------------------*/

#include "mkalf.h"

static int user_decimals = 0;

#define Free_iit(IIT)  basic_iit_kill (IIT);  IIT = NULL;

/*--------------------------------------------------------------------------*/

/* function variables
   for (generic == weighted or unweighted) geometric primitive operations */

void (* alf_calls)() = NULL;
void (* alf_rho0)() = NULL;
void (* alf_rho1)() = NULL;
void (* alf_rho2)() = NULL;
void (* alf_rho3)() = NULL;
int  (* alf_hidden0)() = NULL;
int  (* alf_hidden1)() = NULL;
int  (* alf_hidden2)() = NULL;

/*--------------------------------------------------------------------------*/

Alf_info* alf_info ()
     /* Returns pointer to Alf_info structure wrt current Alf. */
     /* NOTE: The returned address, which points to the info structure,
              is a constant.  DO NOT FREE() IT and consider the fields
              of the structure as read-only. */
{
  static Alf_info ai = { 0 };
  ai = alf_inforec;
  if (alf_context)
    {
      int bytes = sizeof (Alf);
      bytes += (alf_context->dt->n + 1)     * sizeof (Alf_vect); 
      bytes += (alf_context->e_hash_m)      * sizeof (Alf_edge_key);
      bytes += (alf_context->t_hash_m)      * sizeof (Alf_tetra_key);
      bytes += (alf_context->t_hash_m + 1)  * sizeof (Alf_tetra_rank);
      bytes += (alf_context->dt->num.f + 1) * sizeof (Alf_triangle_rank);
      bytes += (alf_context->e_hash_m + 1)  * sizeof (Alf_edge_rank);
      bytes += (alf_context->dt->n + 1)     * sizeof (Alf_vertex_rank);
      bytes += (alf_context->ranks + 1)     * sizeof (Alf_float);
      bytes += (alf_context->ranks + 1)     * sizeof (int);
      bytes += (alf_context->entries + 1)   * sizeof (Alf_master_node);
      bytes += (alf_context->entries + 1)   * sizeof (Alf_master_type);
      ai.bytes = bytes;
      ;
      ai.n =               alf_context->dt->n;
      ai.dt_num =          alf_context->dt->num;
      ai.ranks =           alf_context->ranks;
      ai.entries =         alf_context->entries;
      ai.edge_index_max =  alf_context->e_hash_m;
      ai.tetra_index_max = alf_context->t_hash_m;
      ai.t_hash_m =        alf_context->t_hash_m;
      ai.e_hash_m =        alf_context->e_hash_m;
    }
  ai.e_hash_fnexts =   trist_info ()->min_ef_fnexts;
  return (&ai);
}

/*--------------------------------------------------------------------------*/

void alf_set_decimals (decimals)
     int decimals;
     /* Call this before alf_load_all() to overide the data's maximum
        number of decimals per coordinate. */
{
  if (decimals > 0)
    user_decimals = decimals;
}

/*--------------------------------------------------------------------------*/

Alf_adt alf_load_all (data_path, dt_path, alf_path)
     char data_path[], dt_path[], alf_path[];
     /* Loads coordinates from file data_path, starts SoS (parameter matrix),
        and loads Alf structure from files alf_path and dt_path and makes it
        the current Alf structure (aka alf_context). */
{
  Dt_input_scan *data;
  Alf_adt alp;
  time_t t1 = basic_modtime (data_path);
  time_t t2 = basic_modtime (dt_path);
  time_t t3 = basic_modtime (alf_path);
  if (not ((t1 < t2) and (t2 < t3)))
    print ("WARNING: %s\n %17u %s\n %17u %s\n %17u %s\n\n",
           "Time stamps on files are out of order!",
           t1, data_path, t2, dt_path, t3, alf_path);
  data = dt_input_scan (data_path);
  if (user_decimals > data->decimals)
    { /* cf: alf_set_decimals() */
      data->decimals = user_decimals;
    }
  sos_matrix (data->n, 4,
              data->scale, alf_sos_len (data), alf_sos_lenp (data));
  dt_input_load (data_path, data);
  alp = alf_load (alf_path, dt_path);
  alf_begin ();
  Assert_always ((alf_context == (Alf *) alp) and alf_proper (alf_context));
  return (alp);
}

/*--------------------------------------------------------------------------*/

void alf_kill (alp)
     Alf_adt alp;  /* side effect! */
     /* Deallocate Alf structure... and "shut down" SoS. */
{
  if (alp)
    {
      Alf *a = (Alf *) alp;
      Assert_always (alf_proper (a));
      alf_end ();
      if (alf_context == a)
        alf_context = NULL;
      FREE (a->t_hash_table);
      FREE (a->t_hash_a);
      FREE (a->e_hash_table);
      FREE (a->e_hash_a);
      FREE (a->v_rank);
      FREE (a->e_rank);
      FREE (a->f_rank);
      FREE (a->t_rank);
      dt_kill (a->dt);
      FREE (a->spectrum);
      FREE (a->master_p);
      FREE (a->master_node);
      FREE (a->master_type);
      FREE (a->coord);
      Free_iit (a->t0_iit);
      Free_iit (a->f0_iit);  Free_iit (a->f1_iit);  Free_iit (a->f2_iit);
      Free_iit (a->e0_iit);  Free_iit (a->e1_iit);  Free_iit (a->e2_iit);
      Free_iit (a->v0_iit);  Free_iit (a->v1_iit);  Free_iit (a->v2_iit);
      FREE (a);
    }
  sos_shutdown ();
}

/*--------------------------------------------------------------------------*/

Dt* alf_dt ()
     /* Returns pointer to Dt structure of current Alf structure. */
{
  Assert_always (alf_proper (alf_context));
  return ((Dt *) alf_context->dt);
}

/*--------------------------------------------------------------------------*/

Alf_vect* alf_get_coords ()
     /* Returns pointer to matrix of floating-point coordinates. */
{
  Assert_always (alf_proper (alf_context));
  if (alf_context->coord == NULL)
    { /* The floating-point matrix does not yet exist; load it from SoS. */
      int i, n = alf_context->dt->n;
      Alf_coord scale = (Alf_coord) sos_scale ();
      alf_context->coord = MALLOC (Alf_vect, n + 1);
      upfor (i, 1, n)
        {
          alf_context->coord[i][ALF_X] = scale * (Alf_coord) sos_fp (i, 1);
          alf_context->coord[i][ALF_Y] = scale * (Alf_coord) sos_fp (i, 2);
          alf_context->coord[i][ALF_Z] = scale * (Alf_coord) sos_fp (i, 3);
          if (alf_context->is_weighted)
            { /* re-compute input weight from 4-th coordinate,
                 ie: x^2 + y^2 + z^2 - sign (w) * w^2 */
              lia_push (sos_lia (i, 1));
              lia_ipower (2);
              lia_push (sos_lia (i, 2));
              lia_ipower (2);
              lia_push (sos_lia (i, 3));
              lia_ipower (2);
              lia_plus ();
              lia_plus ();
              lia_push (sos_lia (i, 4));
              lia_minus ();
              alf_context->coord[i][ALF_W] =
                scale * alf_sqrt ((Alf_coord) lia_real (lia_popf ()));
            }
          else
            alf_context->coord[i][ALF_W] = 0.0;
        }
    }
  return (alf_context->coord);
}

/*--------------------------------------------------------------------------*/

int alf_edge_index (ef)
     int ef;
     /* Returns edge index [1..e_hash_m] to edge given by edfacet ef.
        Based on edge hash table [0..e_hash_m-1]. */
{
  int found_flag, h = alf_edge_hash (alf_edge_key (ef), &found_flag);
  Assert_always (found_flag);
  Assert (h + 1 <= alf_context->e_hash_m);
  return (h + 1);
}

/*--------------------------------------------------------------------------*/

int alf_edge_ef (ix)
     int ix;
     /* Returns edfacet representing edge with given index ix. */
{
  Assert ((1 <= ix) and (ix <= alf_context->e_hash_m));
  return (alf_context->e_hash_table[ix-1].min_ef);
}

/*--------------------------------------------------------------------------*/

int alf_tetra_index (ef)
     int ef;
     /* Returns tetra index [1..t_hash_m] to tetrahedra ef^+.
        Based on tetra hash table [0..t_hash_m-1].
        Returns 0 if e^+ outside hull. */
{
  if (trist_hull_facet (ef))
    return (0);
  else
    {
      int found_flag, h = alf_tetra_hash (alf_tetra_key (ef), &found_flag);
      Assert_always (found_flag);
      Assert (h + 1 <= alf_context->t_hash_m);
      return (h + 1);
    }
}

/*--------------------------------------------------------------------------*/

int alf_tetra_ef (ix)
     int ix;
     /* Returns edfacet representing tetrahedron with given index ix. */
{
  Assert ((1 <= ix) and (ix <= alf_context->t_hash_m));
  return (alf_context->t_hash_table[ix-1].min_ef);
}

/*--------------------------------------------------------------------------*/

char *alf_value2str (value)
     Alf_float value;
     /* Takes alpha value and returns *temporary* pointer to string "%e",
        or, if applicable, "%e (-infinity)", "0.0", or "%e (+infinity)".
        NOTE: Use STRDUP() to save temporary strings! */
{
  if (value == -ALF_INFINITY)
    return (basic_cb_frmt ("%e (-infinity)", value));
  else if (value == 0.0)
    return (basic_cb_frmt ("0.0"));
  else if (value == ALF_INFINITY)
    return (basic_cb_frmt ("%e (+infinity)", value));
  else
    return (basic_cb_frmt ("%e", value));
}
