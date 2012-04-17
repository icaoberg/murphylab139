/* sig/sig.h  --- Signatures header file. */

#ifndef __SIG_H__  /* Include this file only once! */
#define __SIG_H__

/*--------------------------------------------------------------------------*/

/*        NOTE: This header files is included inside the "alf.h" file.      */

/*--------------------------------------------------------------------------*/

/* sig.c */

#define SIG_FLOAT 1
#define SIG_INT   2

typedef struct sig_float_type
{
  int high;
  Alf_float max_value;
  Alf_float min_value;
  Alf_float *value; /* [0..high] (**) */
} Sig_float;

typedef struct sig_int_type
{
  int high;
  int max_value;
  int min_value;
  int *value; /* [0..high] (**) */
} Sig_int;

/* (**) A signature s->value[] array goes from index 0 to s->high,
   and:
   - s->high == alf_info ()->ranks, for any signature;
   - s->min_value = Min { s->value[i], 0  < i < s->high };
   - s->max_value = Max { s->value[i], 0  < i < s->high };
   - alpha[0] == -infinity;
   - alpha[s->high] == +infinity */

typedef struct sig_info_type
{
  int  id;
  char *name;
  char *text;
  int  type;
  Sig_int*   (*   int_func) ();
  Sig_float* (* float_func) ();
} Sig_info;

Sig_int*   sig_int_open();
Sig_float* sig_float_open();
void       sig_core();

#define    sig_FREE(S)  do { FREE ((S)->value); FREE (S); } once

void       sig_int_range();
void       sig_float_range();

/* Some metric funtions. */
Alf_float  alf_volume();
Alf_float  alf_triangle_area();
Alf_float  alf_edge_length();
/* ... Actually, they would be quite useful in a more general setting ...
   ... Alf geometric primitives.                                      ... */ 

/* The "core" signatures. */
Sig_float* sig_spectrum();     Sig_info* sig_spectrum_info();
Sig_float* sig_alpha();        Sig_info* sig_alpha_info();

/*--------------------------------------------------------------------------*/

/* metrics.c --- Geometric signatures. */

Sig_float* sig_volume();             Sig_info* sig_volume_info();

Sig_float* sig_area_f_boundary();    Sig_info* sig_area_f_boundary_info();
Sig_float* sig_area_f_regular();     Sig_info* sig_area_f_regular_info();
Sig_float* sig_area_f_singular();    Sig_info* sig_area_f_singular_info();

Sig_float* sig_length_e_singular();  Sig_info* sig_length_e_singular_info();

Sig_float* sig_w();                  Sig_info* sig_w_info();

void       sig_get_metrics_corrections();

/*--------------------------------------------------------------------------*/

/* combinatorics.c --- Combinatorical signatures. */

Sig_int* sig_num_t();           Sig_info* sig_num_t_info();

Sig_int* sig_num_f_boundary();  Sig_info* sig_num_f_boundary_info();
Sig_int* sig_num_f_regular();   Sig_info* sig_num_f_regular_info();
Sig_int* sig_num_f_singular();  Sig_info* sig_num_f_singular_info();
Sig_int* sig_num_f_interior();  Sig_info* sig_num_f_interior_info();

Sig_int* sig_num_e_boundary();  Sig_info* sig_num_e_boundary_info();
Sig_int* sig_num_e_regular();   Sig_info* sig_num_e_regular_info();
Sig_int* sig_num_e_singular();  Sig_info* sig_num_e_singular_info();
Sig_int* sig_num_e_interior();  Sig_info* sig_num_e_interior_info();

Sig_int* sig_num_v_boundary();  Sig_info* sig_num_v_boundary_info();
Sig_int* sig_num_v_regular();   Sig_info* sig_num_v_regular_info();
Sig_int* sig_num_v_singular();  Sig_info* sig_num_v_singular_info();
Sig_int* sig_num_v_interior();  Sig_info* sig_num_v_interior_info();

/*--------------------------------------------------------------------------*/

/* betti.c --- Betti numbers. */

Sig_int* sig_betti_0();  Sig_info* sig_betti_0_info();
Sig_int* sig_betti_1();  Sig_info* sig_betti_1_info();
Sig_int* sig_betti_2();  Sig_info* sig_betti_2_info();

/*--------------------------------------------------------------------------*/

#endif  /* #ifndef __SIG_H__ */
