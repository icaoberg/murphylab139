/* sos/internal.h --- Headerfile for internal SoS routines. */

/* For internal use only! */

#ifndef __SOS_INTERNAL_H__ /* Include this file only once! */
#define __SOS_INTERNAL_H__

/*--------------------------------------------------------------------------*/

#include "basic.h"
#include "sos.h"

/*--------------------------------------------------------------------------*/

#define SOS_MAGIC  130862004

typedef struct sos_common_type
{
  int d, len, lenp;
  int n, n_max;
  double scale;
  int high_minor;
  int *minor;
} SoS_common;

/*--------------------------------------------------------------------------*/

/* sos.c */

extern  SoS_common sos_common;
Lia_ptr sos_lia_0();
Lia_ptr sos_lia_0_0();
void    sos_set_last_star();

/*--------------------------------------------------------------------------*/

/* primitive.c */

void    sos_new_depth_counters();
int     sos_epsilon_compare();

/*--------------------------------------------------------------------------*/

#endif  /* #ifndef __SOS_INTERNAL_H__ */
