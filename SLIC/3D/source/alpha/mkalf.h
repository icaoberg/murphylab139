/* mkalf/mkalf.h --- Common header file for Mkalf and internal Alf core. */

#ifndef __MKALF_H__  /* Include this file only once! */
#define __MKALF_H__

/*--------------------------------------------------------------------------*/

#include "alf.h"

/*--------------------------------------------------------------------------*/

/* internal type definitions: alf.c & mkalf.c */

typedef struct alf_edge_key_type
{
  int min_ef;
} Alf_edge_key;

typedef struct alf_tetra_key_type
{
  int min_ef;
} Alf_tetra_key;

typedef struct alf_vertex_rank_type
{
  int rho, mu1, mu2;
} Alf_vertex_rank;

typedef struct alf_edge_rank_type
{
  int rho, mu1, mu2;
} Alf_edge_rank;

typedef struct alf_triangle_rank_type
{
  int rho, mu1, mu2;
} Alf_triangle_rank;

typedef struct alf_tetra_rank_type
{
  int rho;
} Alf_tetra_rank;

/* NOTE: * generally, 0 <= rho <= mu1 <= mu2 <= ranks;
         * hull vertex, edge, triangle iff mu2 == 0;
         * attached vertex, edge, triangle iff rho == 0;
         * redundant/dumped (or valid) vertex (or face) iff mu1 == mu2 == 0. */

typedef struct alf_master_node_struct
{
  int ix;
} Alf_master_node;

typedef struct alf_master_type_struct
{
  Basic_byte f_type : 3;
  Basic_byte r_type : 2;
  Basic_byte last   : 1;
} Alf_master_type;

/* NOTE!
   These typedefs will be used for the dynamic ARRAYs; see also: spectrum.c.
   To save memory (roughly 18%) it's necessary to keep two parallel arrays
   ARRAY_node[] *and* ARRAY_type[] instead of just one ARRAY[].  Otherwise,
   if the structure is maintained as one record, the f_type and r_type entries
   waste bytes...) */

typedef struct alf_type
{ 
  int magic;                        /*                 filed */
  short float_size;                 /*                 filed */
  short is_weighted;
  Alf_vertex_rank   *v_rank;        /* [1..dt->n],     filed */
  int                e_hash_m;      /* filed */
  int               *e_hash_a;      /* filed */
  Alf_edge_key      *e_hash_table;  /* [0..e_hash_m-1] */
  Alf_edge_rank     *e_rank;        /* [1...e_hash_m], filed */
  Alf_triangle_rank *f_rank;        /* [1..dt->num.f], filed */
  int                t_hash_m;      /*                 filed */
  int               *t_hash_a;      /*                 filed */
  Alf_tetra_key     *t_hash_table;  /* [0..t_hash_m-1] */
  Alf_tetra_rank    *t_rank;        /* [1...t_hash_m]  filed */
  int                entries;       /* filed */
  Alf_master_node   *master_node;   /* [1..entries] \  filed */
  Alf_master_type   *master_type;   /* [1..entries] /  filed */
  int                ranks;         /*                 filed */
  Alf_float         *spectrum;      /* [1..ranks] \    filed */
  int               *master_p;      /* [1..ranks] / */
  Dt *dt;
  double stamp;
  Alf_vect *coord;
  Basic_iit t0_iit;
  Basic_iit f0_iit, f1_iit, f2_iit;
  Basic_iit e0_iit, e1_iit, e2_iit;
  Basic_iit v0_iit, v1_iit, v2_iit;
} Alf;

/*--------------------------------------------------------------------------*/

/* macros for lia_length parameters to sos_matrix(); cf: mkalf.c, alf.c */

#define alf_sos_len(DATA) \
  If ((DATA)->has_weights,  \
      Lia_DIGITS (18 * (DATA)->decimals + 12), \
      Lia_DIGITS (14 * (DATA)->decimals +  8))   /* cf: alf_ratio_compare() */

#define alf_sos_lenp(DATA) \
  Lia_DIGITS (2 * (DATA)->decimals + 1)       /* cf: sum of square + weight */

/*--------------------------------------------------------------------------*/

/* mkalf/internal.c */

extern Alf*     alf_context;
extern Alf_info alf_inforec;

Alf*          alf_malloc();
int           alf_proper();
void          alf_begin();
void          alf_end();
void          alf_save();
Alf_adt       alf_load();

void          alf_edge_hash_build();
int           alf_edge_hash();
Alf_edge_key  alf_edge_key();
int           alf_edge_key_eq();

void          alf_tetra_hash_build();
int           alf_tetra_hash();
Alf_tetra_key alf_tetra_key();
int           alf_tetra_key_eq();

int           alf_ratio_compare();

/*--------------------------------------------------------------------------*/

/* mkalf/weighted.c */

void alf_w_begin();
void alf_w_end();
void alf_w_calls();
void alf_w_size0();
void alf_w_size1();
void alf_w_size2();
void alf_w_size3();
int  alf_w_hidden0();
int  alf_w_hidden1();
int  alf_w_hidden2();

/*--------------------------------------------------------------------------*/

/* mkalf/unweighted.c */

void alf_uw_begin();
void alf_uw_end();
void alf_uw_calls();
void alf_uw_rho0();
void alf_uw_rho1();
void alf_uw_rho2();
void alf_uw_rho3();
int  alf_uw_hidden0();
int  alf_uw_hidden1();
int  alf_uw_hidden2();

/*--------------------------------------------------------------------------*/

/* mkalf/spectrum.c */

void spectrum_open();
void spectrum_tetra();
void spectrum_triangle();
void spectrum_edge();
void spectrum_vertex();
void spectrum_close();
void spectrum_memory_print();

/*--------------------------------------------------------------------------*/

#endif  /* #ifndef __MKALF_H__ */
