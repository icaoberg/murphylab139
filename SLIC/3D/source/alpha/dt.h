/* detri/dt.h --- Dt (and Trist) header file. (This is part of the Alf lib!) */

#ifndef __DT_H__  /* Include this file only once! */
#define __DT_H__ 

/*--------------------------------------------------------------------------*/

#include "trist.h"  /* Include the header file for the Trist module. */

/*--------------------------------------------------------------------------*/

/* file paths for *.dt files, etc: <data_PATH>.<EXTENSION> */

#define       dt_PATH(data_PATH)  basic_cb_frmt ("%s.dt",     data_PATH)
#define     dt_f_PATH(data_PATH)  basic_cb_frmt ("%s.dt_f",   data_PATH)
#define     info_PATH(data_PATH)  basic_cb_frmt ("%s.info",   data_PATH)
#define   status_PATH(data_PATH)  basic_cb_frmt ("%s.status", data_PATH)

#define   tetra_dt_PATH(data_PATH)  basic_cb_frmt ("%s.tl",   data_PATH)
#define tetra_dt_f_PATH(data_PATH)  basic_cb_frmt ("%s.tl_f", data_PATH)

#define fl_default_PATH(data_PATH,TYPE) \
  basic_cb_frmt ("%s.%s.fl", data_PATH,TYPE)

/*** NOTE that above *_PATH macros return pointers to *temporary* strings! ***/
/***      Most of the time you'll have to use it with STRDUP() !!!         ***/

/*--------------------------------------------------------------------------*/

/* detri/dt.c, typedefs */

typedef struct dt_record
{
  short type; /* triangulation type; negative if SoS artefacts were removed! */
  short bpt;
  int n;
  int redundant;  /* n == trist->max_org == num.v + redundant */
  int hull_ef;
  int last_hull;
  Trist_num num;
  Trist *trist;
} Dt;

/* triangulation-type enumerators */
#define DT          101
#define DT_CLOSEST  101
#define DT_FURTHEST 102
#define DT_WEIGHTED 103
#define DT_REGULAR  103

typedef struct dt_input_scan_struct
{
  char *name, *title;
  int lines, n;
  int decimals;
  int fix_w, fix_a;
  double scale;
  int has_weights;
} Dt_input_scan;

/*--------------------------------------------------------------------------*/

/* detri/dt.c */

void dt_save();
Dt*  dt_load();
void dt_kill();

Dt_input_scan* dt_input_scan();
void           dt_input_load();

void dt_test_hooks();
int  dt_test();
int  dt_test_triangle();
int  dt_test_open_tetra();
int  dt_test_pedantic();

#ifdef __DEBUG__
 extern int dt_test_proto_flag;
#endif

void dt_print_info_sec();

/*--------------------------------------------------------------------------*/

/* detri/search.c */

int  dt_search();
int  dt_search_get_tests();

/*--------------------------------------------------------------------------*/

/* other */

void dt_print_cmd();

/*--------------------------------------------------------------------------*/

#endif  /* #ifndef __DT_H__ */
