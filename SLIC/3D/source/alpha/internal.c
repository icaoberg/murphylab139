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
/* mkalf/internal.c --- Alf Library: internal core of Mkalf & Alf. */

/*--------------------------------------------------------------------------*/

#include "mkalf.h"

/*--------------------------------------------------------------------------*/

/* Global variables.  (For internal use only!) */

Alf     *alf_context = NULL;  /* Pointer to current Alf structure. */
Alf_info alf_inforec = { 0 };


/*--------------------------------------------------------------------------*/

/* Local variables and types. */

#define ALF_MAGIC  130862016

typedef union edge_union_type
{
  Alf_edge_key key;
  Basic_byte bytes [sizeof (Alf_edge_key)];
} Edge_union;

typedef union tetra_union_type
{
  Alf_tetra_key key;
  Basic_byte bytes [sizeof (Alf_tetra_key)];
} Tetra_union;

/*--------------------------------------------------------------------------*/

Alf* alf_malloc ()
     /* Allocates Alf record, zeros it, and sets magic number; nothing else! */
{
  Alf *a = MALLOC (Alf, 1);
  BZERO (a,        Alf, 1);
  a->magic = ALF_MAGIC;
  a->float_size = sizeof (Alf_float);
  return (a);
}

/*--------------------------------------------------------------------------*/

int alf_proper (a)
     Alf *a;
     /* Returns TRUE iff Alf *a has right magic number, etc. */
{
  return (     a
          and (a->magic == ALF_MAGIC)
          and (a->float_size = sizeof (Alf_float))
          and (a->dt)
          and (a->is_weighted == (abs (a->dt->type) == DT_WEIGHTED)));
}

/*--------------------------------------------------------------------------*/

void alf_begin ()
{
  if (alf_context->is_weighted)
    {
      alf_w_begin ();
      alf_calls =   alf_w_calls;
      alf_rho0 =    alf_w_size0;
      alf_rho1 =    alf_w_size1;
      alf_rho2 =    alf_w_size2;
      alf_rho3 =    alf_w_size3;
      alf_hidden0 = alf_w_hidden0;
      alf_hidden1 = alf_w_hidden1;
      alf_hidden2 = alf_w_hidden2;
    }
  else
    {
      alf_uw_begin ();
      alf_calls =   alf_uw_calls;
      alf_rho0 =    alf_uw_rho0;
      alf_rho1 =    alf_uw_rho1;
      alf_rho2 =    alf_uw_rho2;
      alf_rho3 =    alf_uw_rho3;
      alf_hidden0 = alf_uw_hidden0;
      alf_hidden1 = alf_uw_hidden1;
      alf_hidden2 = alf_uw_hidden2;
    }
}

/*--------------------------------------------------------------------------*/

void alf_end ()
{
  if (alf_context->is_weighted)
    {
      alf_w_end ();
    }
  else
    {
      alf_uw_end ();
    }
}

/*--------------------------------------------------------------------------*/

void alf_save (alf_path, alp)
     char alf_path[];
     Alf_adt alp;
     /* Writes (Alf *) alp to file alf_path.
        NOTE: ((Alf *) alp)->dt is to be saved seperately! */
{
  Alf *a = (Alf *) alp;
  FILE *f = basic_fopen (alf_path, "w");
  int m, r;
  Assert (a and a->dt);
  print ("Saving %salpha shapes to binary file \"%s\" ...\n",
         If (a->is_weighted, "weighted ", ""), alf_path);
  ;
  Assert_always (alf_proper (a));
  binfwrite (&(a->magic), 1, f);
  binfwrite (&(a->float_size), 1, f);
  ;
  binfwrite (a->v_rank, a->dt->n + 1, f);
  ;
  r = sizeof (Alf_edge_key);
  m = a->e_hash_m;
  binfwrite (a->e_hash_a, r, f);
  binfwrite (&m, 1, f);
  binfwrite (a->e_rank, m + 1, f);
  ;
  binfwrite (a->f_rank, a->dt->num.f + 1, f);
  ;
  r = sizeof (Alf_tetra_key);
  m = a->t_hash_m;
  binfwrite (a->t_hash_a, r, f);
  binfwrite (&m, 1, f);
  binfwrite (a->t_rank, m + 1, f);
  ;
  binfwrite (&(a->entries), 1, f);
  binfwrite (a->master_node, a->entries + 1, f);
  binfwrite (a->master_type, a->entries + 1, f);
  ;
  binfwrite (&(a->ranks), 1, f);
  binfwrite (a->spectrum, a->ranks + 1, f);
  /* master_p[] can be re-computed if all i with master_p[i] == 0 are known;
     NOTE: master_p[i] == 0 iff (i == 0) or (i == ranks);
     Cf: alf_load() */
  ;
  basic_fclose (f);
#if 1
  { /* This test won't be needed in future! */
    int i;
    upfor (i, 1, a->ranks - 1)
      Assert_always ((a->master_p[i] >= 1) and (a->master_p[i] <= a->entries));
    Assert_always ((a->master_p[0] == 0) and (a->master_p[a->ranks] == 0));
  }
#endif
}

/*--------------------------------------------------------------------------*/

Alf_adt alf_load (alf_path, dt_path)
     char alf_path[], dt_path[];
     /* Loads Alf structure (including triangulation!) from files
        alf_path and dt_path and makes it the current Alf structure
        (aka alf_context).
        NOTE: user should call alf_load_all() instead. */
{
  FILE *f;
  int m, r, correct = FALSE;
  int v;
  alf_context = alf_malloc ();
  alf_context->dt = dt_load (dt_path);
  trist_set (alf_context->dt->trist);
  switch (abs (alf_context->dt->type))
    {
     case DT:
     case DT_WEIGHTED:
      correct = TRUE;
    }
  v = trist_vertices ();
  if ((not correct)
      or (alf_context->dt->num.t_flat != 0)
      or (not trist_hull_facet (alf_context->dt->hull_ef))
      or (alf_context->dt->num.v != v))
    basic_error ("load_alf: corrupt triangulation in \"%s\" (%d,%d; %d,%d)",
                 dt_path,
                 alf_context->dt->type, alf_context->dt->num.t_flat,
                 alf_context->dt->num.v, v);
  alf_context->is_weighted = (abs (alf_context->dt->type) == DT_WEIGHTED);
  ;
  f = basic_fopen (alf_path, "r");
  print ("Reading %salpha shapes from binary file \"%s\" ...\n",
         If (alf_context->is_weighted, "weighted ", ""),
         If (basic_fopen_zpath, basic_fopen_zpath, alf_path));
  binfread (&(alf_context->magic), 1, f);
  binfread (&(alf_context->float_size), 1, f);
  if (not alf_proper (alf_context))
    basic_error ("alf_load: wrong magic number (%d,%d)\n",
                 alf_context->magic, alf_context->float_size);
  ;
  alf_context->v_rank = MALLOC (Alf_vertex_rank, alf_context->dt->n + 1);
  binfread (alf_context->v_rank, alf_context->dt->n + 1, f);
  ;
  r = sizeof (Alf_edge_key);
  alf_context->e_hash_a = MALLOC (int, r);
  binfread (alf_context->e_hash_a, r, f);
  binfread (&m, 1, f);
  alf_context->e_hash_m = m;
  alf_context->e_hash_table = MALLOC (Alf_edge_key, m);
  BZERO (alf_context->e_hash_table,   Alf_edge_key, m);
  alf_edge_hash_build (NULL_HOOK);
  alf_context->e_rank = MALLOC (Alf_edge_rank, m + 1);
  binfread (alf_context->e_rank, m + 1, f);
  ;
  alf_context->f_rank = MALLOC (Alf_triangle_rank, alf_context->dt->num.f + 1);
  binfread (alf_context->f_rank, alf_context->dt->num.f + 1, f);
  ;
  r = sizeof (Alf_tetra_key);
  alf_context->t_hash_a = MALLOC (int, r);
  binfread (alf_context->t_hash_a, r, f);
  binfread (&m, 1, f);
  alf_context->t_hash_m = m;
  alf_context->t_hash_table = MALLOC (Alf_tetra_key, m);
  BZERO (alf_context->t_hash_table,   Alf_tetra_key, m);
  alf_tetra_hash_build (NULL_HOOK);
  alf_context->t_rank = MALLOC (Alf_tetra_rank, m + 1);
  binfread (alf_context->t_rank, m + 1, f);
  ;
  binfread (&(alf_context->entries), 1, f);
  alf_context->master_node = MALLOC (Alf_master_node,
                                     alf_context->entries + 1);
  alf_context->master_type = MALLOC (Alf_master_type,
                                     alf_context->entries + 1);
  binfread (alf_context->master_node, alf_context->entries + 1, f);
  binfread (alf_context->master_type, alf_context->entries + 1, f);
  ;
  binfread (&(alf_context->ranks), 1, f);
  alf_context->spectrum = MALLOC (Alf_float, alf_context->ranks + 1);
  alf_context->master_p = MALLOC (int,       alf_context->ranks + 1);
  binfread (alf_context->spectrum, alf_context->ranks + 1, f);
  ;
  { /* re-compute master_p[]; Cf: alf_save() */
    int i, j = 0;
    Alf_master_type *ml_type;
    print ("  Computing master_p[] ...\n");
    alf_context->master_p[0] = 0;
    ml_type = &(alf_context->master_type[j]);
    upfor (i, 1, alf_context->ranks - 1)
      {
        j ++;
        ml_type ++;
        alf_context->master_p[i] = j;
        while (not ml_type->last)
          {
            j ++;
            ml_type ++;
          }
      }
    alf_context->master_p[alf_context->ranks] = 0;
  }
  ;
  /* The time stamp is a heuristic to distinguish different Alf instances.
     The memory address is not enough, since different data structures might
     end up at the same memory allocation...
     Cf. Alvis' "reload" feature, alf_kill(), mkalf/scan.c:initialize(), ... */
  alf_context->stamp = basic_seconds ();
  ;
  /* The following fields stay undefined at this point:
     alf_context->coord;
     alf_context->f1_iit;  alf_context->f2_iit;
     alf_context->e1_iit;  alf_context->e2_iit;
     alf_context->v1_iit;  alf_context->v2_iit; */
  ;
  basic_fclose (f);
  return ((Alf_adt) alf_context);
}

/*--------------------------------------------------------------------------*/

void alf_edge_hash_build (e_hook)
     void (* e_hook) ();
     /* Builds alf_context->e_hash_table[].  If (e_hook != NULL_HOOK),
        "e_hook (ef, ix)" is called for each new entry.
        Assumptions: - alf_context->e_hash_m,e_hash_a[] is already set;
                     - alf_context->e_hash_table[] allocated & zeroed!
        Cf: edges() [mkalf.c], load_alf [alf.c]. */
{
  int found_flag;
  int t, ef, ef1, hx, ix;
  Alf_edge_key k;
  if (not e_hook)
    print ("  Computing e_hash_table[] ...\n");
  Assert_always (alf_proper (alf_context));
  trist_forall (t)
    {
      ef1 = ef = EdFacet (t, 0);
      do
        {
          k = alf_edge_key (ef);
          hx = alf_edge_hash (k, &found_flag);
          ix = hx + 1;
          if (not found_flag)
            {
              Assert_always (hx != -MAXINT);  /* otherwise: hash overflow! */
              alf_context->e_hash_table[hx] = k;
              if (e_hook)
                e_hook (ef, ix);
            }
          ef = Enext (ef);
        } until (ef == ef1);
    }
}

/*--------------------------------------------------------------------------*/

int alf_edge_hash (k, found_flag)
     Alf_edge_key k;
     int *found_flag; /* output */
     /* Returns edge hash index with given key k wrt current Alf.  If output
        parameter *found_flag is set to TRUE, k was found in the hash table;
        OTHERWISE, result is index where k is to be inserted in hash table...
        *UNLESS* -MAXINT is returned: then the table is full! */
{ /* METHOD: Open addressing with double hashing.  
     Cf: Thomas H Cormen, Charles E Leiserson, and Ronald L Rivest.
         "Introduction to Algorithms."  The MIT Press, 1990, p232ff. */
  Edge_union u;
  Alf_edge_key *t;
  static Alf_edge_key nil = { 0 };
  int j, i, ha, hb, m;
  Assert (sizeof (Edge_union) == sizeof (Alf_edge_key));
  u.key = k;
  t = alf_context->e_hash_table;
  m = alf_context->e_hash_m;
  ha = basic_uhash (alf_context->e_hash_a,
                    (int) sizeof (Alf_edge_key), m, u.bytes);
  hb = 1 + ((Abs (k.min_ef)) mod (m - 2));
  alf_inforec.e_hash_queries ++;
  *found_flag = FALSE;
  i = 0;
  j = ha - hb;
  do
    {
      i ++;
      j = (j + hb) mod m;  /* the i-th probe probes t[j] */
      if (alf_edge_key_eq (t[j], k))
        *found_flag = TRUE;
    } until (*found_flag or alf_edge_key_eq (t[j], nil) or (i > m));
  alf_inforec.e_hash_probes += i;
  if (*found_flag)
    return (j);
  else
    return (If ((i == m), -MAXINT, j));
}

/*--------------------------------------------------------------------------*/

Alf_edge_key alf_edge_key (ef)
     int ef;
     /* Returns key <min_ef> for edge given by edge facet ef.
        Cf: trist_edge_min_ef(). */
{
  Alf_edge_key k;
  k.min_ef = trist_edge_min_ef (ef);
  return (k);
}

/*--------------------------------------------------------------------------*/

int alf_edge_key_eq (k1, k2)
     Alf_edge_key k1, k2;
     /* Compares two edge keys. */
{
  return (k1.min_ef == k2.min_ef);
}

/*--------------------------------------------------------------------------*/

void alf_tetra_hash_build (t_hook)
     void (* t_hook) ();
     /* Builds alf_context->t_hash_table[].  If (t_hook != NULL_HOOK),
        "t_hook (ef, ix)" is called for each new entry.
        Assumptions: - alf_context->t_hash_m,t_hash_a[] set;
                     - alf_context->t_hash_table[] allocated & zeroed!
        Cf: tetrahedra() [mkalf.c], load_alf [alf.c]. */
{
  int found_flag;
  int t, ef, v, hx, ix;
  Alf_tetra_key k;
  if (not t_hook)
    print ("  Computing t_hash_table[] ...\n");
  Assert_always (alf_proper (alf_context));
  trist_forall (t)
    upfor (v, 0, 1)
      if (not trist_hull_facet (ef = EdFacet (t, v)))
        {
          k = alf_tetra_key (ef);
          hx = alf_tetra_hash (k, &found_flag);
          ix = hx + 1;
          if (not found_flag)
            {
              Assert_always (hx != -MAXINT);  /* otherwise: hash overflow! */
              alf_context->t_hash_table[hx] = k;
              if (t_hook)
                t_hook (ef, ix);
            }
        }
}

/*--------------------------------------------------------------------------*/

int alf_tetra_hash (k, found_flag)
     Alf_tetra_key k;
     int *found_flag; /* output */
     /* Returns tetra hash index with given key k wrt current Alf.  If output
        parameter *found_flag is set to TRUE, k was found in the hash table;
        OTHERWISE, result is index where k is to be inserted in hash table...
        *UNLESS* -MAXINT is returned: then the table is full! */
{ /* METHOD: open addressing with double hashing.
     Cf: Thomas H Cormen, Charles E Leiserson, and Ronald L Rivest.
         "Introduction to Algorithms."  The MIT Press, 1990, p232ff. */
  Tetra_union u;
  Alf_tetra_key *t;
  static Alf_tetra_key nil = {0};
  int h, i, ha, hb, m;
  Assert (sizeof (Tetra_union) == sizeof (Alf_tetra_key));
  u.key = k;
  t = alf_context->t_hash_table;
  m = alf_context->t_hash_m;
  ha = basic_uhash (alf_context->t_hash_a,
                    (int) sizeof (Alf_tetra_key), m, u.bytes);
  hb = 1 + ((Abs (k.min_ef)) mod (m - 2));
  alf_inforec.t_hash_queries ++;
  *found_flag = FALSE;
  i = 0;
  h = ha - hb;
  do
    {
      i ++;
      h = (h + hb) mod m;  /* the i-th probe probes t[h] */
      if (alf_tetra_key_eq (t[h], k))
        *found_flag = TRUE;
    } until (*found_flag or alf_tetra_key_eq (t[h], nil) or (i > m));
  alf_inforec.t_hash_probes += i;
  if (*found_flag)
    return (h);
  else
    return (If ((i > m), -MAXINT, h));
}

/*--------------------------------------------------------------------------*/

Alf_tetra_key alf_tetra_key (ef)
     int ef;
     /* Returns key <min_ef> for tetra on top of edge facet ef.
        Cf: trist_tetra_min_ef(). */
{
  Alf_tetra_key k;
  k.min_ef = trist_tetra_min_ef (ef);
  return (k);
}

/*--------------------------------------------------------------------------*/

int alf_tetra_key_eq (k1, k2)
     Alf_tetra_key k1, k2;
     /* Compares two tetra keys. */
{
  return (k1.min_ef == k2.min_ef);
}

/*--------------------------------------------------------------------------*/

int alf_ratio_compare (a, b, c, d)
     Lia_obj a, b, c, d;
     /* Compares a/b with c/d returning -1, 0, or 1.
        NOTE: Only special case: 1/0 = infinity ... or ???? */
{
  alf_inforec.ratio_comparisons ++;
  if (lia_eq (a, c) and lia_eq (b, d))  
    return (0);  /* trivially equal */
  else if (lia_sign (b) == 0)
    { /* a/b == infinity */
      Assert_always (FALSE);  /*    <--- does no longer occur, or??? */
      if (lia_sign (d) == 0)
        { /* both are special ratios */
          lia_sub (lia_pushf (LIA_NULL), a, c);
          return (lia_sign (lia_popf ()));
        }
      else
        switch (lia_sign (a))   /* only a/b is special */
          {
           case -1: return (-1);  /* -infinity < anything */
           case  1: return ( 1);  /*  infinity > anything */
           case  0: return (-lia_sign (c));  /* compare 0 to something !=0 */
           default: return (-99);  /* for lint only */
          }
    }
  else if (lia_sign (d) == 0)  /* only second ratio is special */
    {
      Assert_always (FALSE);  /*    <--- does no longer occur, or??? */
      return (-alf_ratio_compare (c, d, a, b));  /* swap ratios! */
    }
  else
    { /* normal case: a/b versus c/d */
      lia_mul (lia_pushf (LIA_NULL), a, d);
      lia_mul (lia_pushf (LIA_NULL), b, c);
      lia_minus ();  /* ad - bc */
      return (lia_sign (lia_popf ()));
    }
}
