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
/* mkalf/ml.c --- Alf Library: interface to master list. */

/*--------------------------------------------------------------------------*/

#include "mkalf.h"

/*--------------------------------------------------------------------------*/

static Alf_master_node *ml_node;
static Alf_master_type *ml_type;
static int f_type;
static int j;

/*--------------------------------------------------------------------------*/

int alf_ml_ranks ()
     /* Returns number of different ranks in current alpha shape.
        NOTE: This function must be called at least once, before you can
        access the master list.  (It initializes module as a side effect!) */
{
  (void) alf_ml_sublist (0);  /* This is tricky!?! */
  return (alf_context->ranks);
}

/*--------------------------------------------------------------------------*/

int alf_ml_sublist (r)
     /* Sets mlp ("master-list pointer") to beginning of sublist for rank r.
        Returns NULL iff sublist is empty.  Valid ranks: 1..alf_ml_ranks(). */
{
  Assert (alf_proper (alf_context) and (0 <= r) and (r <= alf_context->ranks));
  j = alf_context->master_p[r];
  if (j)
    {
      Assert ((0 <= j) and (j <= alf_context->entries));
      ml_node = &(alf_context->master_node[j]);
      ml_type = &(alf_context->master_type[j]);
      f_type = ml_type->f_type;
      return (j);
    }         
  else
    {
      j = If ((r >= alf_context->ranks), alf_context->entries, 0);
      ml_node = &(alf_context->master_node[j]);
      ml_type = &(alf_context->master_type[j]);
      f_type = ALF_BLANK;
      return (NULL);
    }
}

/*--------------------------------------------------------------------------*/

int alf_ml_next ()
     /* Increments mlp. Returns NULL iff end of sublist was reached. */
{
  Assert (ml_node and ml_type);
  if ((f_type == ALF_BLANK) or ml_type->last)
    return (NULL);
  else
    {
      j ++;
      Assert ((0 <= j) and (j <= alf_context->entries));
      ml_node ++;
      ml_type ++;
      f_type = ml_type->f_type;
      return (j);
    }
}

/*--------------------------------------------------------------------------*/

int alf_ml_eofsublist (r)
     /* Sets mlp ("master-list pointer") to end of sublist for rank r.
        Returns NULL iff sublist is empty.  Valid ranks: 1..alf_ml_ranks(). */
{
  Assert ((0 <= r) and (r <= alf_context->ranks));
  j = alf_context->master_p[r];
  if (j)
    {
      Assert ((0 <= j) and (j <= alf_context->entries));
      { /* compute "next" j */
        int j_prime;
        r ++;
        while ((r <= alf_context->ranks) and (not alf_context->master_p[r]))
          r ++;
        j_prime = If ((r > alf_context->ranks),
                      alf_context->entries, alf_context->master_p[r] - 1);
        Assert (j_prime >= j);
        j = j_prime;
      }
      ml_node = &(alf_context->master_node[j]);
      ml_type = &(alf_context->master_type[j]);
      f_type = ml_type->f_type;
      Assert ((0 <= j) and (j <= alf_context->entries));
      return (j);
    }         
  else
    {
      j = If ((r >= alf_context->ranks), alf_context->entries, 0);
      ml_node = &(alf_context->master_node[j]);
      ml_type = &(alf_context->master_type[j]);
      f_type = ALF_BLANK;
      return (NULL);
    }
}

/*--------------------------------------------------------------------------*/

int alf_ml_prev ()
     /* Decrements mlp. Returns NULL iff beginning of sublist was reached. */
{
  Alf_master_type *ml_type_prime = ml_type - 1;
  Assert (ml_node and ml_type);
  if ((f_type == ALF_BLANK) or (ml_type_prime->last))
    return (NULL);
  else 
    {
      j --;
      Assert ((0 <= j) and (j <= alf_context->entries));
      ml_node --;
      ml_type --;
      f_type = ml_type->f_type;
      return (j);
    }
}

/*--------------------------------------------------------------------------*/

int alf_ml_f_type ()
     /* Returns f_type (in {ALF_TETRA, ALF_TRIANGLE, ALF_EDGE, ALF_VERTEX})
        of the face corresponding to current mlp.
        Returns ALF_BLANK, if undefined. */
{
  return (f_type);
}

/*--------------------------------------------------------------------------*/

int alf_ml_r_type ()
     /* Returns r_type (in {ALF_RHO, ALF_MU1, ALF_MU2}) of the face
        corresponding to current mlp. Returns ALF_BLANK, if undefined. */
{
  return (If (f_type == ALF_BLANK, ALF_BLANK, ml_type->r_type));
}

/*--------------------------------------------------------------------------*/

int alf_ml_face ()
     /* Returns the edfacet (or vertex index) of corresponding face
        (or vertex) of current mlp.
        Assumes there is such a face (ie, alf_ml_f_type () != ALF_BLANK);
        otherwise, returns -1; */
{
  switch (f_type)
    {
     case ALF_TETRA:
      return (alf_context->t_hash_table[ml_node->ix - 1].min_ef);
     case ALF_TRIANGLE:
      return (EdFacet (ml_node->ix, 0));
     case ALF_EDGE:
      return (alf_context->e_hash_table[ml_node->ix - 1].min_ef);
     case ALF_VERTEX:
      return (ml_node->ix);
     default:
      return (-1);
    }
}

/*--------------------------------------------------------------------------*/

int alf_ml_index ()
     /* Returns the "index" of face (or vertex) corresponding to current mlp.
        Note:
                tetrahdron ==> alf_tetra_index (alf_ml_face ())
                triangle   ==>         TrIndex (alf_ml_face ())
                edge       ==>  alf_edge_index (alf_ml_face ())
                vertex     ==>                  alf_ml_face ()) */
{
#if __DEBUG__
  switch (f_type)
    {
     case ALF_TETRA:
      Assert (alf_tetra_index (alf_ml_face ()) == ml_node->ix);
      break;
     case ALF_TRIANGLE:
      Assert (TrIndex (alf_ml_face ()) == ml_node->ix);
      break;
     case ALF_EDGE:
      Assert (alf_edge_index (alf_ml_face ()) == ml_node->ix);
      break;
    }
#endif
  return (ml_node->ix);
}

/*--------------------------------------------------------------------------*/

int alf_ml_is_attached ()
     /* Returns TRUE iff face corresponding to current mlp is attached. */
{
  switch (f_type)
    {
     case ALF_TRIANGLE:
      return (alf_context->f_rank[ml_node->ix].rho == 0);
     case ALF_EDGE:
      return (alf_context->e_rank[ml_node->ix].rho == 0);
     case ALF_VERTEX:
      return (alf_context->v_rank[ml_node->ix].rho == 0);
     default:
      return (FALSE);
    }
}

/*--------------------------------------------------------------------------*/

int alf_ml_is_first ()
     /* Returns TRUE iff entry corresponding to current mlp is first
        (leftmost) occurance of corresponding face. */
{
  return (    (alf_ml_r_type () == ALF_RHO)
          or ((alf_ml_r_type () == ALF_MU1) and (alf_ml_is_attached ())));
}
