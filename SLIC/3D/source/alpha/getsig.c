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
#include <stdio.h>
#include "alf.h"

#define NUM_FLOAT_SIGS 7
#define NUM_INT_SIGS 16
#define N_ARGS 2
void print_usage( int argc, char* argv[])
{
	printf( "\n  Usage: %s <inputfile> <outputfile> \n\n",argv[0]);
}

int main( int argc, char *argv[])
{
  int i;
  long num_written;
  long no_of_alphas, num_float_sigs, num_int_sigs;
  FILE *f;
  char *data_name;
  char *output_name;
  char *dt_path;
  char *alf_path;
  Alf_adt alp;
  Sig_float *floatsig, *alpha;
  Sig_int *intsig;
  Sig_float* (*float_sigfunc)( void);
  Sig_int* (*int_sigfunc)( void);
  Sig_float* (*float_sigfunc_array[NUM_FLOAT_SIGS])( void) = 
    {sig_alpha, sig_volume, sig_area_f_boundary, sig_area_f_regular, 
     sig_area_f_singular, sig_length_e_singular, sig_w};
  // Do not use sig_spectrum, it's just alpha^2, not interesting
  Sig_int* (*int_sigfunc_array[NUM_INT_SIGS])( void) = 
    {sig_betti_0, sig_betti_1, sig_betti_2, sig_num_t, sig_num_f_boundary, 
     sig_num_f_singular, sig_num_f_regular, sig_num_f_interior,
     sig_num_e_boundary, sig_num_e_singular, sig_num_e_regular, 
     sig_num_e_interior, sig_num_v_boundary, sig_num_v_singular, 
     sig_num_v_regular, sig_num_v_interior};

  if( argc != N_ARGS+1){ print_usage(argc,argv); return 4;}

  data_name = argv[1];
  output_name = argv[2];
  dt_path = STRDUP( dt_PATH( data_name));
  alf_path = STRDUP( alf_PATH( data_name));
  alp = alf_load_all( data_name, dt_path, alf_path);

  //volume = sig_volume();
  //sigfunc = sigfunc_array[1];
  //volume = sigfunc();
  //printf("volume.high = %i\n", volume->high);
  //printf("volume.max_value = %e\n", volume->max_value);
  //printf("volume.min_value = %e\n", volume->min_value);
  //printf("volume.value = [");
  //for( i = 0; i <= volume->high; i++)
  //  printf("%.5e ", volume->value[i]);
  //printf("]\n");

  // Obtain signatures and write them to output file
  f = fopen( output_name, "wb");
  if( f == NULL) {
    perror("Could not open output file");
    return 1;
  }
  // First write the number of alphas, i.e. the number of
  // elements in each signature
  alpha = sig_alpha();
  no_of_alphas = alpha->high + 1;
  num_written = fwrite( &no_of_alphas, sizeof(no_of_alphas), 1, f);
  if( num_written != 1) {
    perror("Could not write to output file!\n");
    return 2;
  }
  // Then write the number of float signatures
  num_float_sigs = NUM_FLOAT_SIGS;
  num_written = fwrite( &num_float_sigs, sizeof(num_float_sigs), 1, f);
  if( num_written != 1) {
    perror("Could not write to output file!\n");
    return 3;
  }
  // Then write the number of int signatures
  num_int_sigs = NUM_INT_SIGS;
  num_written = fwrite( &num_int_sigs, sizeof(num_int_sigs), 1, f);
  if( num_written != 1) {
    perror("Could not write to output file!\n");
    return 3;
  }
  // Write each float signature
  for( i = 0; i < NUM_FLOAT_SIGS; i++) {
    float_sigfunc = float_sigfunc_array[i];
    floatsig = float_sigfunc();
    num_written = fwrite( floatsig->value, sizeof(Alf_float), no_of_alphas, f);
    if( num_written < no_of_alphas) {
      perror("Could not write to output file!\n");
      return 4;
    }
  }
  // Write each int signature
  for( i = 0; i < NUM_INT_SIGS; i++) {
    int_sigfunc = int_sigfunc_array[i];
    intsig = int_sigfunc();
    num_written = fwrite( intsig->value, sizeof(int), no_of_alphas, f);
    if( num_written < no_of_alphas) {
      perror("Could not write to output file!\n");
      return 4;
    }
  }
  fclose( f);

  return 0;
}
