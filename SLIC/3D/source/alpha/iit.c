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
/* basic/iit.c --- Integer interval tree. */

/*
 * n intervals;
 *   specified by their indices e[1..n]
 *   and endpoint-access function:
 *                low  (e[i])
 *                high (e[i])
 *        is_contained (e[i], x)
 *   eg:
 *        interval e[i]          ... "[" low (e[i]) "," high (e[i]) ")"
 *        is_contained (e[i], x) ... (    (low (e[i] <= x)
 *                                    and (x < high (e[i]))
 */


/*--------------------------------------------------------------------------*/

#include "basic.h"

/*--------------------------------------------------------------------------*/

#define MAGIC  130862017  /* Magic number! */

typedef struct node_type
{
  int value;
  int a, b;
  struct node_type *left, *right;
} Node;

/*
 * NOTE: Node allocations has to be done with standard malloc() rather than
 * MALLOC(), because my "malloc-paranoia" package is really designed for (few)
 * dynamic arrays than for many small nodes (hint: FREE() gets *really* slow :)
 */

typedef struct tree_type
{
  int magic;
  Node *root;
  int *sigma;
  int *omega;
  int (* low) ();
  int (* high) ();
  int (* is_contained)  ();
  Basic_iit_info info;
} Tree;

static Tree *tree = NULL;
static int path_nodes;

static Node* make_node();
static Node* new_node();     
static void __select (int *index, int *value, int i, int j);
static void query();
static int path_length();
static int sigma_compare(), omega_compare();
static void node_kill();

#ifdef __DEBUG__
 static int is_valid_input(), is_valid_output();
 static int is_right_split();
 static int is_sorted(), is_empty();
#endif

int basic_iit_proto_flag = FALSE;  /* Exported flag! */

#define Proto  if (basic_iit_proto_flag) print /* ... */

/*--------------------------------------------------------------------------*/

Basic_iit basic_iit_build (n, e, low, high, is_contained)
     int n;
     int e[]; /* [1..n], input/output */
     int (* low) ();
     int (* high) ();
     int (* is_contained)  ();
     /* Builds interval tree given as mentioned above. */
{
  double time;
  ;
  tree = MALLOC (Tree, 1);
  tree->magic = MAGIC;
  tree->low = low;
  tree->high = high;
  tree->is_contained = is_contained;
  tree->sigma = e;
  tree->omega = MALLOC (int, n + 1);
  BZERO (tree->omega,   int, n + 1);
  ;
  BZERO (&(tree->info), Basic_iit_info, 1);
  tree->info.n = n;
  ;
  /* sigma[1..n] = e[1..n], sorted wrt. low() */
  time = basic_utime ();
  basic_qsort (tree->sigma, 1, n, sigma_compare);
  tree->info.lowsort_time = basic_utime () - time;
  ;
  /* make first node, and recurse... */
  time = basic_utime ();
  tree->root = make_node (1, n);
  tree->info.buildup_time = basic_utime () - time;
  ;
  return ((Basic_iit) tree);
}

/*--------------------------------------------------------------------------*/

static Node* make_node (i, j)
     int i, j;
     /* Makes node out of tree->sigma[i..j], and recurses. */
{
  Assert (is_valid_input (i, j));
  if (i > j)
    return ((Node *) NULL);
  else
    {
      Node *node = new_node ();
      int lambda, iota;
      int q, r;
      ;
      /* select value & get "right of" intervals sigma[r+1..j] */
      __select (&r, &(node->value), i, j);
      ;
      /* mark "contains" intervals in sigma[i..r] to omega[i<q+1..r] */
      q = r;
      downfor (lambda, r, i)
        if (tree->is_contained (node->value, tree->sigma[lambda]))
          {
            tree->omega[q] = tree->sigma[lambda];
            tree->sigma[lambda] = 0;
            q --;
          }
      ;
      /* move remaining "left of" intervals from sigma[i..r] to sigma[i..q] */
      iota = i;
      upfor (lambda, i, r)
        if (tree->sigma[lambda] != 0)
          {
            tree->sigma[iota] = tree->sigma[lambda];
            iota ++;
          }
      Assert (iota == q + 1);
      ;
      /* copy omega[q+1..r] back to sigma[q+1..r] & sort omega[q+1..r] */
      upfor (lambda, q + 1, r)
        {
          tree->sigma[lambda] = tree->omega[lambda];
        }
      basic_qsort (tree->omega, q + 1, r, omega_compare);
      node->a = q + 1;
      node->b = r;
      ;
#if 0
      print (" NODE=%d [%d..%d], left: %d, cont: %d, right: %d\n",
             node->value, i, j, q - i + 1, r - q, j - r);
#endif
      ;
      Assert (is_valid_output (node, i, j));
      ;
      /* recurse */
      node->left  = make_node (i, q);
      node->right = make_node (r + 1, j);
      ;      ;
      return (node);
    }
}

/*--------------------------------------------------------------------------*/

static Node *new_node ()
{
  tree->info.nodes ++;
  return ((Node *) malloc (sizeof (Node)));  /* Node allocation! */
}

/*--------------------------------------------------------------------------*/
static void __select (index, value, i, j)
     int *index;
     int *value;
     int i, j;
     /* Refinement of make_node().
        Select proper value and split off "right of" triangles. */
{
  int r = j - (j - i) / 3;
  int k = tree->low (tree->sigma[r]);
  ;
  while ((r < j) and (tree->low (tree->sigma[r + 1]) == k))
    r ++;
  ;
  if (not (tree->is_contained (k, tree->sigma[r])))
    { /* adjust r to the left for "open" intervals */
      while ((r > i) and (not (tree->is_contained (k, tree->sigma[r - 1]))))
        {
          r --;
          print (" basic_iit: (-)\n");
        }
      if (not tree->is_contained (k, tree->sigma[r]))
        { 
          r --;
          print (" basic_iit: [-]\n");
          Assert_always (r == i - 1);
          print (" basic_iit WARNING: empty NODE!?!\n");
        }
    }
  Assert (is_right_split (i, j, r, k));
  *index = r;
  *value = k;
}

/*--------------------------------------------------------------------------*/

void basic_iit_query (t, x, report)
     Basic_iit t;
     int x;
     void (* report) ();
     /* Query ... */
{
  double time;
  if (basic_iit_proto_flag)
    time = basic_utime ();
  Assert_always (t and (((Tree *) t)->magic == MAGIC));
  tree = (Tree *) t;
  Proto ("<<basic_iit_query: ");
  query (tree->root, x, report);
  Proto ("= %.2fs>>\n", basic_utime () - time);
}

static void query (node, x, report)
     Node *node;
     int x;
     void (* report) ();
{
  int lambda;
  if (node == NULL)
    return;
  else if (x == node->value)
    {
      Proto ("%dD ", node->value);
      upfor (lambda, node->a, node->b)
        report (tree->sigma[lambda]);
    }
  else if (x < node->value)
    {
      Proto ("%dL ", node->value);
      query (node->left, x, report);
      upfor (lambda, node->a, node->b)
        {
          if (tree->is_contained (x, tree->sigma[lambda]))
            report (tree->sigma[lambda]);
          else
            return;
        }
      return;
    }
  else /* (node->value < x) */
    {
      Proto ("%dR ", node->value);
      query (node->right, x, report);
      downfor (lambda, node->b, node->a)
        {
          if (tree->is_contained (x, tree->omega[lambda]))
            report (tree->omega[lambda]);
          else
            return;
        }
      return;
    }
}

/*--------------------------------------------------------------------------*/

void basic_iit_kill (t)
     Basic_iit t;
{
  if (t)
    {
      Tree *tree = (Tree *) t;
      Assert_always (tree->magic == MAGIC);
      FREE (tree->sigma);
      FREE (tree->omega);
      node_kill (tree->root);
      FREE (tree);
    }
}

static void node_kill (node)
     Node *node;
{
  if (node)
    {
      node_kill (node->left);
      node_kill (node->right);
      free (node);  /* Node de-allocation! */
    }
}

/*--------------------------------------------------------------------------*/

Basic_iit_info* basic_iit_info (t)
     Basic_iit t;
{
  Basic_iit_info *inf;
  Assert_always (t and (((Tree *) t)->magic == MAGIC));
  tree = (Tree *) t;
  inf = &(tree->info);
  inf->bytes_array = (inf->n + 1) * sizeof (int) * 2;
  inf->bytes_tree = (inf->nodes) * sizeof (Node);
  inf->bytes = sizeof (Tree) + inf->bytes_tree + inf->bytes_array;
  path_nodes = 0;
  inf->max_depth = 0;
  inf->filled_depth = MAXINT;
  inf->avg_depth = (double) path_length (tree->root, 1) / inf->n;
  inf->avg_nodesize = (double) inf->n / inf->nodes;
  Assert_always (path_nodes == inf->nodes);
  return (inf);
}

static int path_length (node, depth)
     Node *node;
     int depth;
{
  if (node == NULL)
    {
      if (depth  < tree->info.filled_depth)
        tree->info.filled_depth = depth - 1;  /* side effect! */
      return (0);
    }
  else
    {
      path_nodes ++;
      if (depth > tree->info.max_depth)
        tree->info.max_depth = depth;  /* side effect! */
      return (  depth * (node->b - node->a + 1)
              + path_length (node->left,  depth + 1)
              + path_length (node->right, depth + 1));
    }
}

/*--------------------------------------------------------------------------*/

static int sigma_compare (i, j)
     int *i, *j;
     /* Compare tree->sigma[i],[j] wrt. tree->low(); cf: basic_iit_build(). */
{
  int a = tree->low (*i);
  int b = tree->low (*j);
  if (a == b)
    return (0);
  else if (a < b)
    return (-1);
  else
    return (1);
}

static int omega_compare (i, j)
     int *i, *j;
     /* Compare tree->omega[i],[j] wrt. tree->high(); cf: make_node(). */
{
  int a = tree->high (*i);
  int b = tree->high (*j);
  if (a == b)
    return (0);
  else if (a < b)
    return (-1);
  else
    return (1);
}

/*--------------------------------------------------------------------------*/

#ifdef __DEBUG__

static int is_valid_input (i, j)
     int i, j;
{
  int lambda, iota;
  Assert (is_sorted (tree->sigma, i, j, tree->low));
  Assert (is_empty (tree->omega, i, j));
  upfor (lambda, i, j)
    {
      iota = tree->sigma[lambda];
      Assert (   tree->is_contained (tree->low  (iota), iota)
              or tree->is_contained (tree->high (iota), iota));
    }
  return (TRUE);
}

static int is_valid_output (node, i, j)
     Node *node;
     int i, j;
{
  int lambda;
  ;
  Assert ((i <= node->a) and (node->b <= j));
  Assert (is_sorted (tree->sigma, node->a, node->b, tree->low));
  Assert (is_sorted (tree->omega, node->a, node->b, tree->high));
  Assert (is_right_split (i, j, node->b, node->value));
  ;
  upfor (lambda, node->a, node->b)
    Assert (    tree->is_contained (node->value, tree->sigma[lambda])
            and tree->is_contained (node->value, tree->omega[lambda]));
  ;
  upfor (lambda, i, node->a - 1)
    Assert (not tree->is_contained (node->value, tree->sigma[lambda]));
  upfor (lambda, node->b + 1, j)
    Assert (not tree->is_contained (node->value, tree->sigma[lambda]));
  ;
  return (TRUE);
}

static int is_right_split (i, j, r, value)
     int i, j, r;
     int value;
{
  int iota, lambda;
  iota = 0;
  upfor (lambda, i, j)
    if (tree->is_contained (value, tree->sigma[lambda]))
      iota = lambda;
  return (iota == r);
}  

static int is_sorted (array, i, j, endpoint)
     int array[];  /* [i..j] */
     int i, j;
     int (* endpoint) ();
{
  int lambda;
  upfor (lambda, i, j - 1)
    if (endpoint (array[lambda]) > endpoint (array[lambda + 1]))
      return (FALSE);
  return (TRUE);
}

static int is_empty (array, i, j)
     int array[];  /* [i..j] */
     int i, j;
{
  int lambda;
  upfor (lambda, i, j)
    if (array[lambda] != 0)
      return (FALSE);
  return (TRUE);
}

#endif
