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
/* mkalf/mkalf.c  ---  Make Alpha-shape Family structure... from DT. */

/*---------------------------------------------------------------------------*/

static char version[] = "@(#) Mkalf 2.2";
static char purpose[] = "3D Alpha Shape File and Data Structure";
static char  author[] = "Ernst Mucke";

#include "copyright.h"

static char usage[] = "\n\
USAGE: \n\
       %s [OPTIONS] <DATA> \n\
with:  \n\
       <DATA> ... path name of input data \n\
OPTIONS: \n\
       -t[<LEVEL>] .... test option: <LEVEL> = 0 (default), 1, or 2 \n\
ABC mode:         (All But Computation) batch mode: \n\
       -A ............. read commands from standard input, \n\
       -B <FILE> ...... from given <FILE>, or \n\
       -C <COMMAND> ... from command-line. \n\
";

/*--------------------------------------------------------------------------*/

#include "mkalf.h"

/*---------------------------------------------------------------------------*/

/* program flags/parameters and local procedures (w/ common variables) */

static float load_factor = 0.8;  /* load factor of hash tables */

static int test_level = 0;

static char *cmd_name, *data_name, *abc_file_name;

static FILE *info_file;

static void mkalf(), print_info();
static Alf_adt build_alf();
static void tetrahedra(), triangles(), edges(), vertices();
static void t_hook();
static void e_hook();

static void parse_abc_option(), abc_mode();

static int abort_flag = FALSE;  /* see mkalf() and mkalf_error() */
static void mkalf_error();

/*--------------------------------------------------------------------------*/

main (argc, argv)
     int argc;
     char *argv[];
     /* Start-up and command-line options. */
{
  int c, io, is_wrong = FALSE;
  char abc_option = ' ';
  print ("%s --- %s\n%s\n", basic_strip (version), purpose, copyright, author);
  basic_error_hook (mkalf_error);
  basic_malloc_debug (1);
#ifdef __DEBUG__
  trist_test_flag = TRUE;
  trist_io_check_flag = TRUE;
  dt_test_proto_flag = FALSE;
#endif
  cmd_name = argv[0];
  (void) basic_getarg_init (argc - 1, argv + 1);
  while ((c = basic_getarg ("l:t;SAB:C:")) != 0)
    switch (c)
      {
       case 'l':
        { /* -l <LOAD_FACTOR> ... load factor for hash tables */
          if (basic_getarg_optarg)
            io = sscanf (basic_getarg_optarg, "%f", &load_factor);
          if (io != 1)
            is_wrong = TRUE;
          else
            print ("Load factor changed to: %f\n", load_factor);
          break;
        }
       case 't':
        { /* -t[<LEVEL>] ... test option (<LEVEL> = 0, 1, or 2) */
          if (basic_getarg_optarg == NULL)
            test_level = 1;
          else
            sscan (basic_getarg_optarg, "%d", &test_level);
          break;
        }
       case -1:
        { /* <data_name> */
          data_name = basic_getarg_optarg;
          break;
        }
       case 'A':
       case 'B':
       case 'C':
        { /* -A  or
             -B <COMMAND_FILE>  or
             -C <COMMAND_STRING> ... read binaries and execute commands */
          parse_abc_option (&abc_option, c, basic_getarg_optarg);
          break;
        }
       default:
        is_wrong = TRUE;
      }
  if (is_wrong or (not data_name))
    {
      fprint (stderr, usage, cmd_name);
      exit (1);
    }
  if (not (basic_access (data_name)))
    basic_error ("Data file \"%s\" doesn't exist.", data_name);
  abort_flag = TRUE;  /* from now on we core-dump on errors */
  { /* main part */
    if (abc_option == ' ')
      mkalf ();
    else
      abc_mode ();
  }
  basic_malloc_info_print (stdout);
  if (abc_option == 'C')
    unlink (abc_file_name);
  return (0);
}

/*--------------------------------------------------------------------------*/

static void mkalf ()
     /* Mkalf's "main" program. */
{
  double time_total = basic_utime ();
  double time_terminal = basic_seconds ();
  char   *dt_path = STRDUP (  dt_PATH (data_name));
  char  *alf_path = STRDUP ( alf_PATH (data_name));
  char *info_path = STRDUP (info_PATH (data_name));
  Dt_input_scan *data;
  Trist_num qc;
  Alf_adt alp;
  ;
  data = dt_input_scan (data_name);
  sos_matrix (data->n, 4,
              data->scale, alf_sos_len (data), alf_sos_lenp (data));
  dt_input_load (data_name, data);
  ;
  /* open info file; remove old mkalf output, if any, using 'awk' */
  basic_system ("awk '%s' <%s >%s.TEMP; mv -f %s.TEMP %s",
                "BEGIN {q=0} $1==\"#\" {q++; if (q==4) exit} {print}",
                info_path, info_path, info_path, info_path);
  info_file = basic_fopen (info_path, "a");
  fprint (info_file, "#\n# %s (%s), on %s, %s\n#\n  Data: %s\n  Title: %s\n",
          basic_strip (version), cmd_name, basic_hostname (), basic_date (),
          data_name, data->title);
  ;
  alp = build_alf (dt_path, data->has_weights);
  alf_save (alf_path, alp);
  ;
  time_total = basic_utime () - time_total;
  time_terminal = basic_seconds () - time_terminal;
  print_info (time_total, time_terminal);
  qc = ((Alf *) alp)->dt->num;
  if (test_level > 1)
    check_dt (data->n, (int *) NULL, qc);
  alf_kill (alp);
}

/*--------------------------------------------------------------------------*/

static Alf_adt build_alf (dt_path, has_weights)
     char dt_path[];
     int has_weights;
     /* Refinement of mkalf().
        Reads DT from file dt_path and builds Alf structure: alf_context. */
{
  srandom (basic_seed ());  /* for hash tables */
  ;
  alf_context = alf_malloc ();
  alf_context->dt = dt_load (dt_path);
  trist_set (alf_context->dt->trist);
  ;
  if (alf_context->dt->num.t_flat != 0)
    basic_error ("build_alf: triangulation has degenerate tetrahedra");
  if (has_weights and (abs (alf_context->dt->type) == DT))
    basic_error ("build_alf: wrong triangulation type");
  ;
  alf_context->is_weighted = (abs (alf_context->dt->type) == DT_WEIGHTED);
  alf_begin ();
  ;
  spectrum_open ();
  print ("\nComputing hash tables and rho-values:\n");
  tetrahedra ();
  triangles ();
  edges ();
  vertices ();
  spectrum_close (info_file);
  ;
  if (test_level == 1)
    {
      check_alf_hash ();
      check_ranks ();
    }
  return ((Alf_adt) alf_context);
}

/*--------------------------------------------------------------------------*/

static void tetrahedra ()
     /* Allocate/make hash/rank table for tetrahedra and compute rho-values. */
{
  int m, *a, n = alf_context->dt->num.t, r = sizeof (Alf_tetra_key);
  print ("   Tetrahedra...\n");
  n = (int) (((double) n / load_factor) + 0.5);
  if (n < 2)
    n = 2;  /* hashing breaks if num.t == 1 and n==1 (ie, m==2) */
  Assert_always (n >= alf_context->dt->num.t);  /* otherwise: overflow!  */
  a = MALLOC (int, r);
  basic_uhash_new (n, r, &m, a);
  alf_context->t_hash_a = a;
  alf_context->t_hash_m = m;
  alf_context->t_hash_table = MALLOC (Alf_tetra_key, m);
  BZERO (alf_context->t_hash_table,   Alf_tetra_key, m);
  alf_context->t_rank = MALLOC (Alf_tetra_rank, m + 1);
  BZERO (alf_context->t_rank,   Alf_tetra_rank, m + 1);
  alf_tetra_hash_build (t_hook);
}

static void t_hook (ef, ix)
     int ef, ix;
{ 
  spectrum_tetra (ef, ix, &(alf_context->t_rank[ix].rho));
}

/*--------------------------------------------------------------------------*/

static void triangles ()
     /* Allocate rank table for triangles and compute rho-values. */
{
  print ("   Triangles...\n");
  alf_context->f_rank = MALLOC (Alf_triangle_rank, alf_context->dt->num.f + 1);
  BZERO (alf_context->f_rank,   Alf_triangle_rank, alf_context->dt->num.f + 1); 
  {
    int t;
    trist_forall (t)
      spectrum_triangle (EdFacet (t, 0), t, &(alf_context->f_rank[t].rho));
  }
}

/*--------------------------------------------------------------------------*/

static void edges ()
     /* Allocate/make hash/rank table for tetrahedra and compute rho-values. */
{
  int m, *a, n = alf_context->dt->num.e, r = sizeof (Alf_edge_key);
  print ("   Edges...\n");
  n = (int) (((double) n / load_factor) + 0.5);
  Assert_always (n >= alf_context->dt->num.e);  /* otherwise: overflow! */
  a = MALLOC (int, r);
  basic_uhash_new (n, r, &m, a);
  alf_context->e_hash_a = a;
  alf_context->e_hash_m = m;
  alf_context->e_hash_table = MALLOC (Alf_edge_key, m);
  BZERO (alf_context->e_hash_table,   Alf_edge_key, m);
  alf_context->e_rank = MALLOC (Alf_edge_rank, m + 1);  
  BZERO (alf_context->e_rank,   Alf_edge_rank, m + 1);
  alf_edge_hash_build (e_hook);
}

static void e_hook (ef, ix)
     int ef, ix;
{
  spectrum_edge (ef, ix, &(alf_context->e_rank[ix].rho));
}

/*--------------------------------------------------------------------------*/

static void vertices ()
     /* Allocates rank table for vertices and computes rho-values.
        Unweighted case: rho == 0, rank == 1. */
     /* Must be called after edges(). */
{
  int hx, ef, i, j;
  Alf_vertex_rank *rvi, *rvj;
  print ("   Vertices...\n");
  alf_context->v_rank = MALLOC (Alf_vertex_rank, alf_context->dt->n + 1);
  BZERO (alf_context->v_rank,   Alf_vertex_rank, alf_context->dt->n + 1);
  ;
  /* phase 1: set rho -1 (attached), 0 (redundant/dumped), 1 (otherwise) */
  if (alf_context->is_weighted)
    upfor (hx, 0, alf_context->e_hash_m - 1)                /* weighted case */
      {
        ef = alf_context->e_hash_table[hx].min_ef;
        if (ef)
          {
            i = Org (ef);
            j = Dest (ef);
            rvi = &(alf_context->v_rank[i]);
            rvj = &(alf_context->v_rank[j]);
            /* NOTE: if vertex is already known to be attached
               (ie rho==-1), then no further tests are necessary */         
            if (rvi->rho >= 0)
              if (alf_hidden0 (i, j))
                  rvi->rho = -1;
                else
                  rvi->rho = 1;
            if (rvj->rho >= 0)
              if (alf_hidden0 (j, i))
                rvj->rho = -1;
              else
                rvj->rho = 1;
          }
      }
  else
    upfor (hx, 0, alf_context->e_hash_m - 1)              /* unweighted case */
      {
        ef = alf_context->e_hash_table[hx].min_ef;
        if (ef)
          {
            alf_context->v_rank[Org (ef)].rho = 1;
            alf_context->v_rank[Dest (ef)].rho = 1;
          }
        }
  ;
  /* phase 2: set rho -1 (redundant/dumped), 0 (attached), >0 (otherwsie) */
  upfor (i, 1, alf_context->dt->n)
    {
      rvi = &(alf_context->v_rank[i]);
      switch (rvi->rho)
      {
       case 0:
        rvi->rho = -1;  /* redundant/dumped */
        break;
       case -1: 
        rvi->rho = 0;  /* attached */
        break;
       case 1:
        spectrum_vertex (i, &(rvi->rho));
        break;
       default:
        Assert_always (FALSE);
      }
    }
}

/*--------------------------------------------------------------------------*/

static void print_info (total_time, terminal_time)
     double total_time, terminal_time;
     /* Refinement of mkalf(). */
{
  Lia_info *li = lia_info ();
  Alf_info *ainf = alf_info ();
  Trist_info *ti = trist_info ();
  Basic_prime_info *bpi = basic_prime_info ();
  Basic_malloc_info *mi = basic_malloc_info ();
  int i, d, *rho, *hidden, *minor;
  fprint (info_file, "* Lia counters\n");
  fprint (info_file, "%12s . lia_mul calls\n", basic_counter (li->mul_calls));
  fprint (info_file, "%12s . lia_mul elops\n", basic_counter (li->mul_elops));
  fprint (info_file, "%12s . Lia add calls\n", basic_counter (li->padd_calls));
  fprint (info_file, "%12s . Lia add elops\n", basic_counter (li->padd_elops));
  fprint (info_file, "%12s . Lia sub calls\n", basic_counter (li->psub_calls));
  fprint (info_file, "%12s . Lia sub elops\n", basic_counter (li->psub_elops));
  fprint (info_file, "%12d . Lia maximum digit\n", li->maximum_digit);
  fprint (info_file, "* Primitive calls\n");
  alf_calls (&rho, &hidden, &d);
  upfor (i, 0, d)
    fprint (info_file, "%12d . Rho/Size%d\n", rho[i], i);
  upfor (i, 0, d - 1)
    fprint (info_file, "%12d . Hidden%d\n", hidden[i], i, i);
  fprint (info_file, "* Minor evaluations\n");
  sos_minor_calls (&minor, &d);
  upfor (i, 0, d)
    fprint (info_file, "%12d . %d-by-%d\n", minor[i], i, i);
  fprint (info_file, "* Hashing\n");
  fprint (info_file, "%12.2f . %s, %10d / %10d\n",
          (double) ainf->e_hash_fnexts / ainf->e_hash_queries,    
          " edge fnexts/query",
          ainf->e_hash_fnexts, ainf->e_hash_queries);
  fprint (info_file, "%12.2f . %s, %10d / %10d, load factor %.2f\n",
          (double) ainf->e_hash_probes / ainf->e_hash_queries,
          " edge probes/query",
          ainf->e_hash_probes, ainf->e_hash_queries,
          (double) alf_context->dt->num.e / alf_context->e_hash_m);
  fprint (info_file, "%12.2f . %s, %10d / %10d, load factor %.2f\n",
          (double) ainf->t_hash_probes / ainf->t_hash_queries,
          "tetra probes/query",
          ainf->t_hash_probes, ainf->t_hash_queries,
          (double) alf_context->dt->num.t / alf_context->t_hash_m);
  fprint (info_file, "%12d . prime tests\n", bpi->tests);
  fprint (info_file, "%12d . prime modulos\n", bpi->mods);
  fprint (info_file, "* Miscellaneous\n");
  fprint (info_file, "%12.3f . Mb SoS (parameter matrix)\n",
          basic_mbytes (sos_bytes ()));
  fprint (info_file, "%12.3f . Mb Trist (%d bytes per triangle)\n",
          basic_mbytes (ti->bytes), ti->bpt);
  fprint (info_file, "%12.3f . Mb Alf\n", basic_mbytes (ainf->bytes));
  fprint (info_file, "%12.3f . Mb malloc, total (arena: %.3f MB)\n",
          basic_mbytes (mi->total), basic_mbytes (mi->arena));
  dt_print_info_sec (info_file, total_time,    " CPU secs, mkalf total");
  dt_print_info_sec (info_file, terminal_time, "real secs, mkalf total");
}

/*--------------------------------------------------------------------------*/

static void parse_abc_option (current, option, string)
     char *current;  /* input/output */
     char option;
     char string[];
     /* Refinement of main(). */
{
  static char internal[30];
  FILE *f;
  if ((*current != ' ') and (*current != option))
    basic_error ("Can't use -A, -B, and -C options together: -%c %s",
                 option, string);
  if ((*current == ' ') and (option == 'C'))
    { /* start new internal command file */
      sprintf (internal, ".mkalf.%d.abc", basic_seed ());
      abc_file_name = internal;
      f = basic_fopen (abc_file_name, "w");
      fprint (f, "%s\n", string);
      basic_fclose (f);
    }
  else if (option == 'C')
    { /* continue internal command file */
      f = basic_fopen (abc_file_name, "a");
      fprint (f, "%s\n", string);
      basic_fclose (f);
    }
  else if (option == 'B')
    { /* take given command file */
      abc_file_name = string;      
    }
  else
    { /* read from standard input */
      abc_file_name = ".";
    }
  *current = option;
}

/*--------------------------------------------------------------------------*/

static void abc_mode ()
     /* Mkalf [ABC] batch mode. (All but computation. :) */
{
  int i;
  FILE *f = If ((strcmp (abc_file_name, ".") == 0),
                stdin, basic_fopen (abc_file_name, "r"));
  char *buf, *token[99];
  char *dt_path = STRDUP (dt_PATH (data_name));
  char *alf_path = STRDUP (alf_PATH (data_name));
  Alf_adt alp = alf_load_all (data_name, dt_path, alf_path);
  if (test_level > 0)
    {
      Trist_num qc;
      qc = ((Alf *) alp)->dt->num;
      print ("\nChecking *.alf and *.dt file structure.\n");
      check_dt (trist_max_org (), (int *) NULL, qc);
      check_alf_hash ();
      check_ranks ();
      print ("Okay.\n\n");
    }
  print ("Entering [ABC] batch mode...\n", basic_strip (version));
  while (buf = basic_cb_getline (f))
    {
      buf[strlen (buf) - 1] = '\0';  /* delete '\n' from basic_cb_getline() */
      i = basic_tokenize (buf, token, 99);
      if (basic_arg_find ("help", i, token))
        {
          dt_print_cmd (data_name, alf_dt (), i, token);
          alf_print_cmd (data_name, i, token);
        }
      else if (basic_arg_find ("print", i, token))
        {
          if (alf_is_print_cmd (i, token))
            alf_print_cmd (data_name, i, token);
          else 
            dt_print_cmd (data_name, alf_dt (), i, token);
        }
    }
  if (f and (f != stdin))
    basic_fclose (f);
  alf_kill (alp);
}

/*--------------------------------------------------------------------------*/

static void mkalf_error (message)
     char message[];
{
  fprint (stderr, "\n%s: %s\n", cmd_name, message);
  if (abort_flag)
    abort ();
  else
    exit (1);
}
