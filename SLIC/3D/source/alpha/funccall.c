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
void initialize_bc_model(BAYES_CLASSIFIER *bc,
                         /* method for initializing joint density estimator */
                         void (*init_model) (JDE_MODEL *, long unsigned int, lon
g unsigned int *,
                                             long unsigned int, double *),
                         /* method for adding an example to joint density estima
tor */

      bc->init_model(&(bc->jde[i]), m, num_possible_values, n_minus_m, bin_toler
ance_ver_1);

typedef struct {
  /* method for initializing joint density estimator */
  void (*init_model) (JDE_MODEL *, long unsigned int, long unsigned int *,
                      long unsigned int, double *);
  /* method for adding an example to joint density estimator */
  void (*add_to_model) (JDE_MODEL *, DV *, CV *);
  /* method for preparing joint density estimator for prediction */
  void (*compute_model) (JDE_MODEL *);
