/* mkalf/alf.h  --- Alf header file. */

#ifndef __ALF_H__  /* Include this file only once! */
#define __ALF_H__

/*--------------------------------------------------------------------------*/

#include "basic.h"
#include "lia.h"
#include "sos.h"
#include "trist.h"
#include "dt.h"

/*--------------------------------------------------------------------------*/

/* Mkalf file paths: <data_PATH>.<EXTENSION> */

#define  alf_PATH(data_PATH)  basic_cb_frmt ("%s.alf", data_PATH)
#define sigs_PATH(data_PATH)  basic_cb_frmt ("%s.sig", data_PATH)

/* NOTE that above *_PATH macros return *temporary* strings only! */

/*--------------------------------------------------------------------------*/

/* data types */

typedef  char*  Alf_adt;  /* abstract data type! */

typedef struct alf_info_type
{
  int bytes;
  int n;
  Trist_num dt_num;
  int ranks, entries;
  int edge_index_max, tetra_index_max;
  int ratio_comparisons;
  int e_hash_fnexts;
  int e_hash_queries, e_hash_probes;
  int t_hash_queries, t_hash_probes;
  int e_hash_m, t_hash_m;
} Alf_info;

typedef  float      Alf_float;
typedef  float      Alf_coord;    /* same as "Coord" type in SGI's <gl.h> */
typedef  Alf_coord  Alf_vect[4];  /* ALF_X,Y,Z(,W) coordinates */

/*--------------------------------------------------------------------------*/

/* constants, symbols, enumerators */

/* indexing ALF_X,Y,Z,W coordinates */
#define ALF_X 0
#define ALF_Y 1
#define ALF_Z 2
#define ALF_W 3

/* some floating-point constants */
#define ALF_MAXIMUM  ((Alf_float) MAXFLOAT)
#define ALF_MINIMUM  ((Alf_float) MINFLOAT)
#define ALF_ZERO     ((Alf_float) 0.0)
#define ALF_INFINITY ((Alf_float) MAXFLOAT)

/* f_type/r_type enumerators; f_type == dimension + 1 */
#define ALF_TETRA    4
#define ALF_TRIANGLE 3
#define ALF_EDGE     2
#define ALF_VERTEX   1    /* ^                                       */
#define ALF_BLANK    0    /* |- 3 bits f_type  and  2 bits |- r_type */
#define ALF_RHO      1    /*                               V         */
#define ALF_MU1      2
#define ALF_MU2      3

/*--------------------------------------------------------------------------*/

/* mkalf/alf.c */

Alf_info*  alf_info();
void       alf_set_decimals();
Alf_adt    alf_load_all();
void       alf_kill();

Dt*        alf_dt();
Alf_vect*  alf_get_coords();

int        alf_edge_index();
int        alf_edge_ef();

int        alf_tetra_index();
int        alf_tetra_ef();

char*      alf_value2str();

/* generic primitives */
extern void (* alf_calls)();
extern void (* alf_rho0)();
extern void (* alf_rho1)();
extern void (* alf_rho2)();
extern void (* alf_rho3)();
extern int  (* alf_hidden0)();
extern int  (* alf_hidden1)();
extern int  (* alf_hidden2)();

/*--------------------------------------------------------------------------*/

/* alf/ml.c */

int  alf_ml_ranks();
int  alf_ml_sublist();
int  alf_ml_next();
int  alf_ml_eofsublist();
int  alf_ml_prev();
int  alf_ml_f_type();
int  alf_ml_r_type();
int  alf_ml_face();
int  alf_ml_index();
int  alf_ml_is_attached();
int  alf_ml_is_first();

/*--------------------------------------------------------------------------*/

/* mkalf/lookup.c */

Alf_float  alf_sqrt();
Alf_float  alf_threshold();
Alf_float  alf_threshold2();
void       alf_interval();
void       alf_interval2();
int        alf_rank();
int        alf_rank2();

int alf_is_redundant();
int alf_is_valid();
int alf_is_singular();
int alf_is_interior();
int alf_is_in_complex();
int alf_is_on_hull();

/*--------------------------------------------------------------------------*/

/* mkalf/print_alf.c */

void alf_print_cmd();
int  alf_is_print_cmd();

/*--------------------------------------------------------------------------*/

/* mkalf/check.c */

void check_dt();
void check_alf_hash();
void check_ranks(); 

/*--------------------------------------------------------------------------*/

/* mkalf/scan.c */

void alf_scan_alpha_f1(), alf_scan_alpha_f2();
void alf_scan_alpha_e1(), alf_scan_alpha_e2();
void alf_scan_alpha_v1(), alf_scan_alpha_v2();

extern char* alf_scan_iit_info_path;
extern int   alf_scan_iit_proto_flag;

/*--------------------------------------------------------------------------*/

/* mkalf/old_scan.c */

void alf_scan_alpha();  /* obsolete! */
void alf_scan_diff();

/*--------------------------------------------------------------------------*/

#include "sig.h"   /* Include the signature stuff. */

/*--------------------------------------------------------------------------*/

#endif  /* #ifndef __ALF_H__ */
