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
/***************************************************************************/
/**************************    betti.c    **********************************/
/***************************************************************************/
/*  This version was completed 2/25/94 at 6:00pm                           */
/*  It is edited for alf version 2.0.                                      */
/*  It returns only the betti number of the desired dimension.             */
/*  betti is allocated to hold only one signature.                         */
/*  The numbers beside some comments correspond to the numbered            */
/*  statements in the general algorithm.                                 */
/***************************************************************************/

#ifndef lint
 static char author[] = "Cecil Delfinado";
#define _no_status_
#endif

/* See the end of this file for the exported functions. */

#include "alf.h"


/***************************************************************************/
/***************************************************************************/
/****  Definition of External Static Variables, Types and Functions  *******/
/***************************************************************************/
/***************************************************************************/

#define NIL -1  /* used as null pointers in linked-lists                   */

#define Switch(i, j) { int tmp; tmp = i; i = j; j = tmp; }
                /* used to implement weight balancing in the               */
                /* union find data structure                               */

typedef struct uf_record_type  /* union-find data structure data type      */
{
  int p;        /* if this node is not the root                            */
                /*   then this is the index of the parent of this node     */
                /*   else this is the negative of                          */
                /*        the number of descendants + 1 of this root node  */
} Uf_record;

static Alf_info* info;               /* pointer to structure Alf_info      */

static int* ml_marks;                /* pointer to an array of integers    */

static Uf_record* ufds;              /* union-find data structure          */
                                     /* used first with points as elements */
                                     /* then with tetrahedra as elements   */

static Sig_int* betti;               /* will point to a 1 element array of */
                                     /* Sig_int corresponding either to    */
                                     /* betti0, betti1 or betti2           */


/***************************************************************************/
/***************************************************************************/
/****  Useful procedures for the Union-Find data structure         *********/
/***************************************************************************/
/***************************************************************************/


/***************************************************************************/
/*  (1.1.1) and (1.2.1)                                                    */
/*  This procedure creates a yet empty union-find structure of size n      */
/***************************************************************************/

static void uf_create (n)
     int n;
{
  if (ufds == NULL)  
      ufds = MALLOC(Uf_record, n);
  else            /* if called a second time */
    REALLOC(ufds, Uf_record, n); 

  BZERO(ufds, Uf_record, n); 
}


/***************************************************************************/
/*  (1.1.8) and (1.2.11)                                                   */
/*  This procedure returns (the name of) the set that contains the         */
/*  simplex indexed by i and does path compression when possible           */
/***************************************************************************/

static int uf_find (i)
     int i;
{
  int j,t;
  j = i;
  while (ufds[j].p >= 0) {     /* find the root of this set                */
    j = ufds[j].p;
  }                            /* j is the index of the root of set of i   */
  while (ufds[i].p >= 0) {     /* PATH COMPRESSION :  make nodes on path   */
    t = i;                     /*   from i to root have root as parent     */
    i = ufds[i].p;
    ufds[t].p = j;
  }
  return (j);                  /* j is the index of the root of set of i   */
}


/***************************************************************************/
/*  (1.1.2), (1.2.2) and (1.2.6)                                           */
/*  This procedure adds the simplex indexed by i as the only element       */
/*  of a new set to the system                                             */
/***************************************************************************/

static void uf_add (i)
     int i;
{
  ufds[i].p = -1;               /* element i defines a set by itself       */
                                /* and has size 1                          */
}


/***************************************************************************/
/*  (1.1.8) and (1.2.11)                                                   */
/*  This procedure tries to form the union of two sets containing          */
/*  elements indexed by i and j.  If the union is actually performed       */
/*  1 is returned, 0 otherwise                                             */
/***************************************************************************/

static int uf_union (i, j)
     int i;
     int j;
{
  i = uf_find(i);                /* find the root of set containing i      */
  j = uf_find(j);                /* find the root of set containing j      */
  if (i == j) { 
    return (0);                  /* no union effected                      */
  }
  else                           /* if different sets                      */
    {
                                 /* WEIGHT BALANCING                       */
      if (ufds[i].p              /* ufds[i].p = -1 *                       */
                                 /*            number of descendants of i  */
          > ufds[j].p             
         )
         Switch(i, j);           /* i must be the larger set               */
      ufds[i].p = ufds[i].p + ufds[j].p;
                                 /* no. of descendants of i increases      */
      ufds[j].p = i;             /* j now points to i                      */

      return (1);                /* union effected                         */
    }
}


/***************************************************************************/
/***************************************************************************/
/****  Useful procedures for the Alpha shape data structure        *********/
/***************************************************************************/
/***************************************************************************/


/***************************************************************************/
/*  (1.1.7)                                                                */
/*  given a pointer to an edge in the master_list                          */
/*  this procedure returns the indices of the endpoints of this edge       */
/*  by making i and j point to these vertices                              */
/***************************************************************************/

static void endpoints (i, j)
     int * i;
     int * j;
{
  *i = Org (alf_ml_face());   
               /* index of the origin vertex in the vertex rank table      */
  *j = Dest(alf_ml_face());
               /* index of the destination vertex in the vertex rank table */
}


/***************************************************************************/
/***************************************************************************/
/****  Main procedures for betti computation                       *********/
/***************************************************************************/
/***************************************************************************/


/***************************************************************************/
/*  This procedure marks the array pointed to by ml_marks to compute       */
/*    betti_1.  The marks are                                              */
/*           0 if ml entry is new but unmarked, and                        */
/*           1 if ml entry is new but marked.                              */
/*          -1 if ml entry is otherwise.                                   */
/***************************************************************************/
/***************************************************************************/

static void mark_master_list_1 ()
{
  int rank,                    /* the current rank number of alpha shape   */
      ranks,                   /* the number of ranks of the alpha shape   */
      mlix,                    /* linearized master list index             */
      t1,t2,                   /* pointers to tetrahedra                   */
      v1,v2,                   /* pointers to vertices                     */
      ef,                      /* edge-facet of a triangle                 */
      did_unite;               /* if union effected                        */

  /******* (1.1)  forward direction  ***************************************/

  uf_create(info->n+1);  
               /* (1.1.1) create union find data structure with n sets     */
               /* the indices used are indices of the undumped vertices    */
 
  ranks = alf_ml_ranks (); /* set ranks and initialize master list pointer */
  mlix = 0;

  /* (1.1.3) Traverse the master list ranks in increasing rank number      */
  upfor (rank, 1, ranks) {
    if (alf_ml_sublist (rank))
      do {    /* (1.1.4) Traverse the elements in a master-list rank       */
              /* from first element in rank to last element in rank        */
              
        mlix++;    /* counter to make ml_marks parallel with master-list   */

              /* Process the new vertices in this interval.                */

        if (alf_ml_f_type() == ALF_VERTEX) 
          if (alf_ml_is_first()) 
            {
              ml_marks[mlix] = 1;                               /* mark it */
              uf_add(alf_ml_index());   
                         /* add the vertex in union find data structure    */
            }
          else 
            {
              ml_marks[mlix] = -1;                       /* not first time */
            }

        else
        if (alf_ml_f_type() == ALF_EDGE)                       /* (1.1.5)  */
          if (alf_ml_is_first())
                         /* (1.1.6) Check if first occurrence              */
            {
              endpoints(&v1,&v2);
                         /* (1.1.7) Find the endpoints                     */
              did_unite = uf_union(v1,v2);
                         /* (1.1.8) Unite the set containing the endpoints */
              if (did_unite)                                   /* (1.1.9)  */  
                ml_marks[mlix] = 0;          /* (1.1.10) union effected    */
              else 
                ml_marks[mlix] = 1;          /* (1.1.11) no union effected */
                                                               /* (1.1.12) */
            }                                                
          else {
            ml_marks[mlix] = -1;             /* (1.1.13) not first time    */
          }                                                    /* (1.1.14) */
                                                               /* (1.1.16) */
      }
      while (alf_ml_next ());                                  /* (1.1.15) */
  }                                                            /* (1.1.17) */
                                                               /* (1.1.18) */

  /*******  (1.2)  backward direction  *************************************/

  uf_create(info->t_hash_m+1);  
             /* (1.2.1) create union find data structure with t_hash_m sets*/
             /* the indices used are indices of the tetrahedra rank table  */

  uf_add(0); /* (1.2.2) add the solid outside the convex hull              */

  ranks = alf_ml_ranks (); /* set ranks and initialize master list pointer */
  mlix = info->entries + 1;

  /* (1.2.3) Traverse the master list ranks in decreasing rank number      */
  downfor (rank, ranks, 1) {
    if (alf_ml_eofsublist (rank))
      do {    /* (1.2.4) Traverse the elements in a master-list rank       */
              /* from last element in rank to first element in rank        */

        mlix--;    /* counter to make ml_marks parallel with master-list   */

        if (alf_ml_f_type() == ALF_TETRA)                      /* (1.2.5)  */
              /* all tetrahedra appear in the master list exactly once     */
          {

            ml_marks[mlix] = 0;     
              /* (1.2.7) all bounded tetrahedra are unmarked               */
            uf_add(alf_ml_index());                            /* (1.2.6)  */
          }                                                   

        else
        if (alf_ml_f_type() == ALF_TRIANGLE)                   /* (1.2.8)  */
          if (alf_ml_is_first())
                         /* (1.2.9)  Check if first occurrence             */
            { ef = alf_ml_face();
              t1 = alf_tetra_index(ef);
              t2 = alf_tetra_index(Sym(ef));
                         /* (1.2.10) Find the incident tetrahedra          */
              did_unite = uf_union(t1,t2);
                         /* (1.2.11) Unite the set containing              */
                         /* the incident tetrahedra                        */
              if (did_unite)                                   /* (1.2.12) */
                ml_marks[mlix] = 1;         /* (1.2.13) union effected     */
              else 
                ml_marks[mlix] = 0;         /* (1.2.14) no union effected  */
            }                                                  /* (1.2.15) */
          else 
            {
              ml_marks[mlix] = -1;           /* (1.2.16) not first time    */
            }                                                  /* (1.2.17) */
                                                               /* (1.2.19) */
                                                               /* (1.2.20) */
      }
      while (alf_ml_prev ());                                  /* (1.2.18) */
  }                                                            /* (1.2.21) */
                                                               /* (1.2.22) */
  FREE (ufds);                                                    /* (3.2) */
  ufds = NULL;
}

/***************************************************************************/
/* The following procedures compute the betti signatures.                  */
/* For dimension 0 and 2, the signatures are computed directly.            */
/* For dimension 1, the simplices are first marked before sig is computed. */
/***************************************************************************/

static void build_betti_0 (sm)
     Sig_info *sm;
{
  int rank,                    /* the current rank number of alpha shape   */
      ranks,                   /* the number of ranks of the alpha shape   */
      v1,v2,                   /* pointers to vertices                     */
      did_unite;               /* if union effected                        */

  /******* (1.1)  forward direction  ***************************************/

  ranks = alf_ml_ranks (); /* set ranks and initialize master list pointer */

  betti = sig_int_open (sm);
    /* betti signature for betti 0 */
  betti->min_value = 0;
  betti->max_value = 0;
  betti->value[1]  = 0;

  uf_create(info->n+1);  
               /* (1.1.1) create union find data structure with n sets     */
               /* the indices used are indices of the undumped vertices    */
 
  /* (1.1.3) Traverse the master list ranks in increasing rank number      */
  upfor (rank, 1, ranks) {
    if (alf_ml_sublist (rank)) {
      do {    /* (1.1.4) Traverse the elements in a master-list rank       */
              /* from first element in rank to last element in rank        */
              
        if (alf_ml_f_type() == ALF_VERTEX) {
          if (alf_ml_is_first()) 
            {
              uf_add(alf_ml_index());   
                         /* add the vertex in union find data structure    */
              betti->value[rank]++;
            }
        }
        else
          if (alf_ml_f_type() == ALF_EDGE)                     /* (1.1.5)  */
            if (alf_ml_is_first())
                         /* (1.1.6) Check if first occurrence              */
              {
                endpoints(&v1,&v2);
                         /* (1.1.7) Find the endpoints                     */
                did_unite = uf_union(v1,v2);
                         /* (1.1.8) Unite the set containing the endpoints */
                if (did_unite)                                 /* (1.1.9)  */
                  betti->value[rank]--;      /* (1.1.10) union effected    */
              }                                                
      }                                                        /* (1.1.16) */
      while (alf_ml_next ());                                  /* (1.1.15) */
      if (betti->value[rank] > betti->max_value)
        betti->max_value = betti->value[rank];
    }
    if (rank < ranks)
      betti->value[rank+1] = betti->value[rank];
  }                                                            /* (1.1.17) */
                                                               /* (1.1.18) */
  FREE (ufds);                                                    /* (3.2) */
  ufds = NULL;
}


/***************************************************************************/

static void build_betti_1 (sm)
     Sig_info *sm;
{
  int dim,                              /* the dimension of a simplex      */
      mlix,                             /* linearized master list index    */
      rank,                             /* the current alpha shape rank    */
      ranks;                            /* the number of alpha shape ranks */
      
  betti = sig_int_open (sm);
      /*  betti signature for betti1 */

  /* initialize some signature fields */
  betti->min_value = 0;
  betti->max_value = 0;
  betti->value[1]  = 0;

  /* (2.8) traverse array of marks parallel to master list                 */

  ranks = alf_ml_ranks ();
  mlix = 0;
 
  upfor (rank, 1, ranks) {
    if (alf_ml_sublist(rank)) {
      do {
        mlix++;    /* counter to make ml_marks parallel with master-list   */

        if (ml_marks[mlix]>=0)                                   /* (2.9)  */
          {
            switch (alf_ml_f_type()) {                           /* (2.10) */
              case ALF_VERTEX:   dim = 0; break;
              case ALF_EDGE:     dim = 1; break;
              case ALF_TRIANGLE: dim = 2; break;
              case ALF_TETRA:    dim = 3; 
            }
            if (dim == 1) {                                      /* (2.11) */
              if (ml_marks[mlix] == 1)
                betti->value[rank]++;                            /* (2.12) */
            }
            else if (dim == 2)
              if (ml_marks[mlix] == 0)
                betti->value[rank]--;                            /* (2.13) */
                                                                 /* (2.14) */
          }                                                      /* (2.15) */
      }
      while (alf_ml_next ());                                    /* (2.16) */
      if (betti->value[rank] > betti->max_value)                 /* (2.18) */
        betti->max_value = betti->value[rank];
    }
    if (rank < ranks)
      betti->value[rank+1] = betti->value[rank];                 /* (2.17) */
  }
}                                                                /* (2.21) */


/***************************************************************************/

static void build_betti_2 (sm)
     Sig_info *sm;
{
  int rank,                    /* the current rank number of alpha shape   */
      ranks,                   /* the number of ranks of the alpha shape   */
      t1,t2,                   /* pointers to tetrahedra                   */
      ef,                      /* edge-facet of a triangle                 */
      did_unite;               /* if union effected                        */

  /*******  (1.2)  backward direction  *************************************/

  ranks = alf_ml_ranks (); /* set ranks and initialize master list pointer */

  betti = sig_int_open (sm);
    /* betti signature for betti 2 */
  betti->min_value = 0;
  betti->max_value = 0;
  betti->value[ranks]  = 0;

  uf_create(info->t_hash_m+1);  
             /* (1.2.1) create union find data structure with t_hash_m sets*/
             /* the indices used are indices of the tetrahedra rank table  */

  uf_add(0); /* (1.2.2) add the solid outside the convex hull              */


  /* (1.2.3) Traverse the master list ranks in decreasing rank number      */
  downfor (rank, ranks, 1) {
    betti->value[rank-1] = betti->value[rank];
    if (alf_ml_eofsublist (rank)) {
      do {    /* (1.2.4) Traverse the elements in a master-list rank       */
              /* from last element in rank to first element in rank        */

        if (alf_ml_f_type() == ALF_TETRA)                      /* (1.2.5)  */
              /* all tetrahedra appear in the master list exactly once     */
          {
            uf_add(alf_ml_index());                            /* (1.2.6)  */
            betti->value[rank-1]++;
          }                                                   

        else
        if (alf_ml_f_type() == ALF_TRIANGLE)                   /* (1.2.8)  */
          if (alf_ml_is_first())
                         /* (1.2.9)  Check if first occurrence             */
            { ef = alf_ml_face();
              t1 = alf_tetra_index(ef);
              t2 = alf_tetra_index(Sym(ef));
                         /* (1.2.10) Find the incident tetrahedra          */
              did_unite = uf_union(t1,t2);
                         /* (1.2.11) Unite the set containing              */
                         /* the incident tetrahedra                        */
              if (did_unite)                                   /* (1.2.12) */
                betti->value[rank-1]--;
            }                                                  /* (1.2.15) */
      }
      while (alf_ml_prev ());                                  /* (1.2.18) */
      if (betti->value[rank-1] > betti->max_value)
        betti->max_value = betti->value[rank-1];
    }
  }                                                            /* (1.2.21) */
                                                               /* (1.2.22) */
  FREE (ufds);                                                    /* (3.2) */
  ufds = NULL;
}


/***************************************************************************/
/*  These procedures are called from the outside to compute the            */
/*  betti number signature of given dimension.                             */
/***************************************************************************/

Sig_info* sig_betti_0_info ()
{
  static Sig_info sm = {0};
  sm.id   = 3;
  sm.name = "#components";
  sm.text = "Betti[0], number of connected components.";
  sm.type = SIG_INT;
  sm.int_func = sig_betti_0;
  return (&sm);
}

Sig_int* sig_betti_0 ()
{
  info = alf_info ();
  build_betti_0 (sig_betti_0_info ());
  return (betti);
}

Sig_info* sig_betti_1_info ()
{
  static Sig_info sm = {0};
  sm.id   = 4;
  sm.name = "#tunnels";
  sm.text = "Betti[1], number of independent non-bounding 1-cycles.";
  sm.type = SIG_INT;
  sm.int_func = sig_betti_1;
  return (&sm);
}

Sig_int* sig_betti_1 ()
{
  info = alf_info ();
  ml_marks = MALLOC (int, info->entries + 1);
  mark_master_list_1 ();
  build_betti_1 (sig_betti_1_info ());
  FREE (ml_marks); 
  return (betti);
}

Sig_info* sig_betti_2_info ()
{
  static Sig_info sm = {0};
  sm.id   = 5;
  sm.name = "#voids";
  sm.text = "Betti[2], number of independent non-bounding 2-cycles.";
  sm.type = SIG_INT;
  sm.int_func = sig_betti_2;
  return (&sm);
}

Sig_int* sig_betti_2 ()
{
  info = alf_info ();
  build_betti_2(sig_betti_2_info ());
  return (betti);
}


/***************************************************************************/
/******************    END OF betti.c   ************************************/
/***************************************************************************/

