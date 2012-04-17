/* detri/trist.h --- Trist header file. */

#ifndef __TRIST_H__  /* Include this file only once! */
#define __TRIST_H__ 

/*--------------------------------------------------------------------------*/

#include "basic.h"
#include "sos.h"

/*--------------------------------------------------------------------------*/

/* trist.c typedefs */

typedef  unsigned int  Trist_org;

/* NOTE: Trist_org is the type of the origin[0..2] entries in the Trist_record.
         It is assumed to be either "unsigned short" or "unsigned int".
         In the 1st case, there are 30 bytes per triangles;
         in the 2nd case, 36 bytes, and MAX_ORG will be identical to MAXINT. */

typedef struct trist_record_type
{
  Trist_org origin[3];
  int        fnext[6];
} Trist_record;

typedef struct trist_type
{
  int magic;
  int max_org;
  int max_triangle;
  Trist_record *triangle;  /* [1..last_triangle]; entry [0] is unused! */
  int last_triangle;
  int used_triangles;
  int next_reusable_triangle;
  int data_size;
  char *data;  /* pointer to optional data array, parallel to triangle[] */
} Trist;

typedef struct trist_info_type
{
  Basic_counter orgs, syms, enexts, fnexts, fsplices;
  int min_ef_fnexts;
  int bpt, bytes, maxbytes;
} Trist_info;

typedef struct trist_num_type
{
  int t, f, e, v;
  int fh, eh, vh;
  int t_proper;
  int t_flat;
} Trist_num;

/*--------------------------------------------------------------------------*/

/* some macros for convenience */

#define trist_Div8  >> 3
#define trist_Mul8  << 3
#define trist_Mod8  & 07

/* NOTE: These bit operators should only be used in trist.c & macros below. */

#define trist_trindex(E)   ( (E) trist_Div8)
#define trist_trversion(E) ( (E) trist_Mod8)
#define trist_edfacet(T,V) (((T) trist_Mul8) + (V))

/* NOTE: Above 3 macros act like functions! */

/*--------------------------------------------------------------------------*/

/* loop macros, scanning through the Trist structure */

#define trist_for(T)    upfor (T, 1, trist_last ()) if (not trist_deleted (T))
#define trist_forall(T) upfor (T, 1, trist_last ())

/* NOTE: - trist_for (t) { ... } loops through triangle records t,
           but doesn't execute block if trist_deleted (t) == TRUE;
           eg, Detri needs to use this!
         - trist_forall (t) { ... } loops through all triangle records t
           without implicitly calling trist_deleted();
           this works if Trist is packed; eg, if it's loaded from a file like
           in Mkalf (in this case, trist_deleted() is always FALSE! */

/*--------------------------------------------------------------------------*/

/* trist.c --- core routines */

int    trist_upper_bound();
Trist* trist_alloc();
void   trist_kill();
void   trist_set();
void   trist_copy();
void   trist_clear();
Trist* trist_current();

int  trist_max();
int  trist_last();
int  trist_vertices();
int  trist_max_org();

/* NOTE: Above routines are implemented as functions rather than macros
   which will cost a little bit (eg, 3 secs out of 500) ... but it's a
   lot cleaner. */

void trist_max_org_set();

int  trist_make();
void trist_delete();
int  trist_deleted();

int trist_org();
int trist_sym();
int trist_enext();
int trist_fnext();

int trist_enext2();
int trist_turn();
int trist_dest();

void trist_fsplice();

void trist_triangle();
void trist_tetra();

void trist_hull_facet_set();
int  trist_hull_triangle();
int  trist_hull_facet();

int trist_tetra_min_ef();
int trist_edge_min_ef();

int         trist_vertices();
Trist_num   trist_num();
int         trist_num_eq();

Trist_info* trist_info();

/*--------------------------------------------------------------------------*/

/* trist.h --- auxiliary functions */

void trist_print();

int  trist_pack();
int  trist_pack_n_keep();
int  trist_pack_n_sort();
void trist_permute();
void trist_permute_hook();
void trist_hull_copy();

void trist_color_set();
char trist_color();

void  trist_data_size();
char *trist_data_addr();
void  trist_data_zero();

/*--------------------------------------------------------------------------*/

/* trist.c --- normal vectors */

void trist_nvx_push();
void trist_nvx();
void trist_nvx_pop();

/*--------------------------------------------------------------------------*/

/* peel.c */

int trist_peel();

/*--------------------------------------------------------------------------*/

/* trist.c --- debugging */

#ifdef __DEBUG__
 extern int trist_test_flag;
 extern int trist_io_check_flag;
       void trist_io_check();
#endif

/*--------------------------------------------------------------------------*/

/* abbreviations for the most used Trist functions */

#define Org trist_org
#define Sym trist_sym
#define Enext trist_enext
#define Fnext trist_fnext
#define Enext2 trist_enext2
#define Turn trist_turn
#define Dest trist_dest
#define Fsplice trist_fsplice

#define TrIndex trist_trindex
#define TrVersion trist_trversion
#define EdFacet trist_edfacet

/*--------------------------------------------------------------------------*/

/*     ........................................................................
       trist, sb.[1] Obs. Also 3-5 triste, 4-5 tryst(e, (5 thrist t)). [App.
       etymologically related to TRAIST, TRUST; but the nature of the relation
       is not clear; see further under TRUST sb.] Confidence, faith; confident
       expectation, hope: = TRUST sb. 1, 2.
       --pat/OED2 server on oed2.cso.uiuc.edu (Oxford English Dictionary)
       ......................................................................*/

#endif  /* #ifndef __TRIST_H__ */
