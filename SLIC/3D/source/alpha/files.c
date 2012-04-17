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
/* basic/files.c --- Encapsulated file utilities.  Cf, man fopen. */

/*--------------------------------------------------------------------------*/

#define AUTOMATIC_UNCOMPRESS  /* Uncomment line to turn feature off. */

/* If AUTOMATIC_UNCOMPRESS is on, than basic_fopen (pathname, "r") will first
   look for a compressed file, uncompress it (via a Unix pipe) and read it
   (from the pipe).  Is the mode == "w", basic_fopen() will first delete the
   compressed file, and then create a new file pathname.

   The first call of basic_fopen() will decide on what compression scheme is
   used. */

/*--------------------------------------------------------------------------*/

#include "basic.h"
#include <unistd.h>
#include <sys/stat.h>

/*--------------------------------------------------------------------------*/

#ifdef AUTOMATIC_UNCOMPRESS

typedef struct plist_record
{
  int id;
  struct plist_record *next;
} Plist;
static Plist *plist = NULL;  /* local linked list remembers files that were
                                open with popen(); NOTE: no de-allocation! */

static char *format[]  = { "", "%s.Z", "%s.z",      "%s.gz",      0 };
static char *command[] = { "", "zcat", "gunzip -c", "gzip -d -c", 0 };
static int scheme = 0;
static void select_scheme();  /* select from supported compression schemes */

#endif

char *basic_fopen_zpath = NULL;  /* exports actual file path of last
                                    basic_fopen() call, *if* file was
                                    compressed & user wants to know */

/*--------------------------------------------------------------------------*/

FILE* basic_fopen (path_name, mode)
     char path_name[], mode[];
     /* Opens file with given path_name in given mode. Cf, man fopen. */
{
  FILE *file;
  char fname[MAXPATHLEN];
  char fmode[20];
  sprint (fname, "%s", path_name);
  sprint (fmode, "%s", mode);
#ifndef AUTOMATIC_UNCOMPRESS
  {
    file = fopen (fname, fmode);
    if (not file)
      basic_error ("basic_fopen: Can't open file \"%s\" with mode \"%s\".", 
                   fname, fmode);
  }
#else
  { /* check for compressed file first */
    static char zname[MAXPATHLEN];
    basic_fopen_zpath = NULL;
    if (scheme == 0)
      select_scheme (&scheme, fname);
    sprint (zname, format[scheme], fname);
    if (access (zname, F_OK) != -1)
      {
        if (strpbrk (fmode, "wW"))
          {
            if (access (zname, W_OK) == -1)
              basic_error ("basic_fopen: Can't remove file \"%s\".",
                           zname);
            print ("Removing compressed file \"%s\" ...\n", zname);
            basic_system ("rm -f %s", zname);
            file = fopen (fname, fmode);
          }
        else
          { /* open compressed file for reading */
            char cmd[20+MAXPATHLEN];
            sprintf (cmd, "%s %s", command[scheme], zname);
            file = popen (cmd, "r");
            if (file)
              { /* push file id onto plist, use malloc() to hide it */
                Plist *new = (Plist *) malloc (sizeof (Plist));
                new->id = fileno (file);
                new->next = plist;
                plist = new;
              }
            basic_fopen_zpath = zname;  /* export compressed file name,
                                           in case caller wants to know */
          }
      }
    else
      file = fopen (fname, fmode);
    if (not file)
      basic_error ("basic_fopen: Can't open file \"%s\", %s \"%s\".",
                   fname, "with mode", fmode);
  }
#endif
  return (file);
}

/*--------------------------------------------------------------------------*/

void basic_fclose (file)
     FILE *file;
     /* Closes the given file.  Cf, man fclose. */
{
#ifndef AUTOMATIC_UNCOMPRESS
  if (fclose (file) != 0)
    fprint (stderr, "WARNING: basic_fclose: Unsucessful\n");
#else
  Plist *p = plist;
  while (p != NULL)
    { /* if file is in plist, use pclose() */
      if (p->id == fileno (file))
        { 
          if (pclose (file) != 0)
            fprint (stderr, "WARNING: basic_fclose: Unsucessful pclose()\n");
          p->id = 0;  /* forget about the file id */
          return;
        }
      p = p->next;
    }
  /* otherwise, use fclose() */
  if (fclose (file) != 0)
    fprint (stderr, "WARNING: basic_fclose: Unsucessful fclose()\n");
#endif
}

/*--------------------------------------------------------------------------*/

#ifdef AUTOMATIC_UNCOMPRESS

static void select_scheme (scheme, fname)
     int *scheme; /* output */
     char fname[];
{
  char zname[MAXPATHLEN];
  int i = 1;
  while (format[i])
    {
      sprint (zname, format[i], fname);
      if (access (zname, F_OK) != -1)
        {
          *scheme = i;
#if 1
          print ("(Using \"%s %s\" for decompression.)\n",
                 command[i], format[i]);
#endif    
          break;
        }
      i ++;
    }
}

#endif

/*--------------------------------------------------------------------------*/

int basic_access (path_name)
     char path_name[];
     /* Checks if file with given path_name exists. Cf, man access. */
{
#ifndef AUTOMATIC_UNCOMPRESS
  return (access (path_name, F_OK) != -1);
#else
  if (access (path_name, F_OK) != -1)
    return (TRUE);
  else
    { /* check for path_name.Z */
      char zname[MAXPATHLEN];
      if (scheme == 0)
        select_scheme (&scheme, path_name);
      sprint (zname, format[scheme], path_name);
      return (access (zname, F_OK) != -1); 
    }
#endif
}

/*--------------------------------------------------------------------------*/

time_t basic_modtime (path_name)
     char path_name[];
     /* Returns last modifaction time of file, measured in (time_t) seconds
        since 00:00:00, Jan 1, 1970.  NOTE: Figuring out whether the file
        is compressed or not, is a hack! */
{
  time_t no_access = (time_t) 0;
  if (not basic_access (path_name))
    return (no_access);
  else
    {
      struct stat buf;
      if (stat (path_name, &buf) == 0)
        return (buf.st_mtime);
      else
        { /* check for compressed file; spaghetti code... */
          char zname[MAXPATHLEN];
          sprint (zname, format[scheme], path_name);
          if (stat (zname, &buf) == 0)
            return (buf.st_mtime);
          else
            return (no_access);
        }
    }
}
