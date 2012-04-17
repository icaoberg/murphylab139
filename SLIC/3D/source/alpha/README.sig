_________________________________________________________________
A brief memo on how to access the alpha-shape signature functions

Assume a file <DATA> with the coordinates of input points . The
alpha-shape software will produce two files <DATA>.dt (storing the
Delaunay triangulation, generated via Detri and <DATA>.alf (essentially
storing a set of intervals, for each face of the Delaunay triangulation).
Together, the files represent the whole family of alpha shapes. Eg,

        % detri <DATA>
        % mkalf <DATA>
        % your_code <DATA>

Your code can then load this representation into main memory:

        #include "alf.h"

        char *data_name = "<DATA>";
        char *dt_path = STRDUP (dt_PATH (data_name));
        char *alf_path = STRDUP (alf_PATH (data_name));
        Alf_adt alp = alf_load_all (data_name, dt_path, alf_path)

A floating-point signature, say, the volume, you would load like this:

        Sig_float *volume = sig_volume ();

Here, "loading" means computing the signature function from the
alpha-shape representation, and allocating the approrpiate memory for
it.  The deallocation routine is:

        sig_FREE (volume);

The type Sig_float is defined as:

        typedef struct sig_float_type
        {
          int high;
          Alf_float max_value;
          Alf_float min_value;
          Alf_float *value; /* [0..high] */
        } Sig_float;

An integer signature can be loaded analogously, but has a different
type specification:

        typedef struct sig_int_type
        {
          int high;
          int max_value;
          int min_value;
          int *value; /* [0..high] */
        } Sig_int;

Both, integer and floating-point signatures can be deallocated using
the sig_FREE() macro.

A signature s->value[] array goes from index 0 (rank 0, "empty set"
--- rank 1 denotes the point set for unweighted alpha shapes) through
index s->high (== rank alf_info ()->ranks, denoting "infinity"):

           - alpha[0] == -infinity;
           - alpha[s->high] == +infinity
           - s->high == alf_info ()->ranks, for any signature;

and:

           - s->min_value = Min { s->value[i], 0  < i < s->high };
           - s->max_value = Max { s->value[i], 0  < i < s->high };

Tha alpha value alpha[r] corresponding to the shape of rank r can also
be found in a signature function: sig_alpha().

And now the list of currently implemented signatures:

  Sig_Float* sig_volume();
  Sig_Float* sig_area_f_boundary();
  Sig_Float* sig_area_f_regular();
  Sig_Float* sig_area_f_singular();
  Sig_Float* sig_length_e_singular();
  Sig_Float* sig_w();

  Sig_int*  sig_betti_0();
  Sig_int*  sig_betti_1();
  Sig_int*  sig_betti_2();

  Sig_int*  sig_num_t();
  Sig_int*  sig_num_f_boundary();
  Sig_int*  sig_num_f_singular();
  Sig_int*  sig_num_f_regular();
  Sig_int*  sig_num_f_interior();
  Sig_int*  sig_num_e_boundary();
  Sig_int*  sig_num_e_singular();
  Sig_int*  sig_num_e_regular();
  Sig_int*  sig_num_e_interior();
  Sig_int*  sig_num_v_boundary();
  Sig_int*  sig_num_v_singular();
  Sig_int*  sig_num_v_regular();
  Sig_int*  sig_num_v_interior();

  Sig_float* sig_alpha();
  Sig_float* sig_spectrum();

Refer to the source files in the sig/ subdirectory for the
corresponding "one-liner" of documentation.
