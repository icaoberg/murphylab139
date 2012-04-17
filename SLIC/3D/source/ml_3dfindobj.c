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
/*///////////////////////////////////////////////////////////////////////////////
//
// ml_3dfindobj.c - finds objects in a 3D image, where objects are 26-connected
//                regions of above-threshold voxels.
//
// Copyright (c) 2000-2002 Meel Velliste, Carnegie Mellon University
//
// created:	13-14 Jan 2000
//
// modification history:
//              22 Jun 2000 MV: Added code to increase stack size to cope
//                              with larger images
//              05 Sept 2000 MV: Separated feature calculations from obj-
//                               ect finding. Created a MEX interface and
//                               removed standalone interface
//              Spring 2002 MV: Implemented initial background detection, so as
//                              to avoid building up a huge stack when finding
//                              the background as a hole.
//              17-24 May 2002 MV: Implemented a standalone interface again
//              05 July 2002 MV: Added the option of excluding objects below
//                               a certain size. Improved memory efficiency.
//                               Improved source code layout by grouping 
//                               functions for easier readability.
//
// Version 1.8
//
////////////////////////////////////////////////////////////////////////////////
// Every voxel is an unsigned char or a signed/unsigned short.
// Image is gray level type (0..255) or (0..4095) or (0..16383) or (0..65535) etc.
// Thresholding: WHITE = GrayLevel > THRESHOLD, BLACK = GrayLevel < THRESHOLD.
// GrayLevel = THRESHOLD is used to encode voxels that have already been considered
// so that they will not be counted twice.
// Therefore initially all voxels which = THRESHOLD are set to THRESHOLD+1; */

#define inline
#define MEGABYTE 0x100000
#define STACK_MEGABYTES 3000

#ifdef _MEX_
#include "mex.h"
#include "matrix.h"
#endif /*/ #ifdef _MEX_ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <memory.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/types.h>

/*/ global variables (DO NOT FORGET TO INITIALIZE THESE) */
int TotalObjects;
int TotalSizeOfObjects;
int TotalHoles;
int TotalSizeOfHoles;
int findholes; /*/ This global flag indicates if holes are to be found */
int min_obj_size;

#define COORD_type signed short
#define SIZE COORD_type
#define DATA_8bit unsigned char
#define HEADER_SIZE 16

#define EXIT	exit(1);

/*/ MACROS */
#define THRESHOLD 3
#define WHITE 1
#define BLACK 0
#define BKGND 2
#define idx(x,y,z)              (y) + (x)*Im->height + (z)*Im->height*Im->width
#define is_white(x,y,z)         Im->p[idx((x),(y),(z))] == WHITE
#define is_black(x,y,z)         Im->p[idx((x),(y),(z))] == BLACK
#define is_bkgnd(x,y,z)         Im->p[idx((x),(y),(z))] == BKGND
#define set_voxel(x,y,z,val)    Im->p[idx((x),(y),(z))] = (val)

struct HOLE {
  unsigned long size; /*/ Number of voxels that make up the hole */
  struct VOXEL *pVoxels; /*/ List of voxel coordinates in this hole */
  struct HOLE *pNext;
  int part_of_bkgnd; /*/ Boolean: Whether this hole is actually  */
  /*/ part of the background */
};

#define REAL_TYPE double
struct COORD_REAL {
	REAL_TYPE x;
	REAL_TYPE y;
	REAL_TYPE z;
};

struct OBJECT {
  unsigned long size; /*/ Number of voxels in the object*/
  struct VOXEL *pVoxels; /*/ List of voxel coordinates in this object */
  unsigned long n_holes; /*/ Number of holes in the object */
  struct HOLE *pHoleList;
  struct HOLE *pLastHole;
  struct OBJECT *pNext;
};

struct IMAGE {
  long height; /*/ y-size */
  long width;  /*/ x-size */
  long depth;  /*/ z-size */
  long bytes_per_voxel;
  long threshold; /* The graylevel value at which image is to be thresholded */
  struct COORD_REAL COF; /*/ Center of fluorescence of the image */
  DATA_8bit *p; /*/ pointer to image data*/
};

/*/ Declare a global image pointer. This is used by the recursive routines 
// so that it does not have to be passed on stack */
struct IMAGE *Im;

struct VOXEL {
        COORD_type x;
        COORD_type y;
        COORD_type z;
        struct VOXEL *pNext;
};

/*
struct COORD {
	COORD_type x;
	COORD_type y;
	COORD_type z;
};
*/

/*/////////////////////////////////////////////////////////////////////
// FUNCTIONS FOR ADDING ELEMENTS TO THE VARIOUS LISTS*/

void add_voxel( struct VOXEL **v, COORD_type x, COORD_type y, COORD_type z)
{ /*/ Add voxel coordinates to the list */
  struct VOXEL *pNewVoxel;

  /*/ Allocate memory */
	pNewVoxel = (struct VOXEL*) malloc( sizeof( struct VOXEL));
	if( pNewVoxel == NULL){ printf( "Out of memory!\n"); EXIT;}
	/*/ Fill in the coordinates */
	pNewVoxel->x = x;
	pNewVoxel->y = y;
	pNewVoxel->z = z;
	/*/ Link new member to the beginning of the list */
	pNewVoxel->pNext = *v;
	*v = pNewVoxel;
}

void add_voxel_to_object( struct OBJECT *o, COORD_type x, COORD_type y, COORD_type z)
{
	o->size++;
	add_voxel( &(o->pVoxels), x, y, z);
}

void add_voxel_to_hole( struct HOLE *h, COORD_type x, COORD_type y, COORD_type z)
{
	h->size++;
	add_voxel( &(h->pVoxels), x, y, z);
}

struct HOLE *add_hole_to_object( struct OBJECT *pObject)
{
  /*/ Allocate memory */
	struct HOLE *pNewHole = (struct HOLE*) malloc( sizeof( struct HOLE));
	if( pNewHole == NULL){ printf( "Out of memory!\n"); EXIT;}
	/*/ Initialise new hole */
	memset( pNewHole, 0, sizeof( struct HOLE));
	/*/ Link new member to the end of the list */
	if( pObject->pLastHole == NULL)
		pObject->pHoleList = pNewHole;
	else {
		pObject->pLastHole->pNext = pNewHole;
	}
	pObject->pLastHole = pNewHole;
	pObject->n_holes++;
	return pNewHole;
}

struct OBJECT *new_object( void)
{
  /*/ Allocate memory */
	struct OBJECT *pNewObject = malloc( sizeof( struct OBJECT));
	if( pNewObject == NULL){ printf( "Out of memory!\n"); EXIT;}
	/*/ Initialise new object */
	memset( pNewObject, 0, sizeof( struct OBJECT));
	return pNewObject;
}

/*///////////////////////////////////////////////////////////////////
// FUNCTIONS FOR REMOVING ELEMENTS FROM THE VARIOUS LISTS */

void free_voxels( struct VOXEL *pVoxels)
{
	struct VOXEL *pThis, *pNext;

	pThis = pVoxels;
	while( pThis != NULL) {
		pNext = pThis->pNext;
		free( pThis);
		pThis = pNext;
	}
}

void free_holes( struct HOLE *pHoles)
{
	struct HOLE *pThis, *pNext;

	pThis = pHoles;
	while( pThis != NULL) {
	        free_voxels( pThis->pVoxels);
		pNext = pThis->pNext;
		free( pThis);
		pThis = pNext;
	}
}

void free_objects( struct OBJECT *pObjects)
{
	struct OBJECT *pThis, *pNext;

	pThis = pObjects;
	while( pThis != NULL) {
		pNext = pThis->pNext;
		free_voxels( pThis->pVoxels);
		free_holes( pThis->pHoleList);
		free( pThis);
		pThis = pNext;
	}
}

void free_image( struct IMAGE *Image)
{
	free( Image->p);
}

/*/////////////////////////////////////////////////////////////////////
// FUNCTIONS FOR OBJECT AND HOLE FINDING */

struct HOLE *h;
void findhole( COORD_type x, COORD_type y, COORD_type z)
     /*/ Recursively finds hole at given voxel (x,y,z) */
{
	static COORD_type X, Y, Z;

	add_voxel_to_hole( h, x, y, z);
	/*/ make sure the current voxel wont be done again*/
	set_voxel(x, y, z, THRESHOLD);
	TotalSizeOfHoles++;
	/*/ Scan 6-neighborhood for connected voxels */
	X = x - 1;
	if( X >= 0) 
		if( is_black( X, y, z)) findhole( X, y, z);
		else if( is_bkgnd( X,y,z)) h->part_of_bkgnd = 1;
	X = x + 1;
	if( X < Im->width) 
		if( is_black( X, y, z)) findhole( X, y, z);
		else if( is_bkgnd( X,y,z)) h->part_of_bkgnd = 1;
	Y = y - 1;
	if( Y >= 0)
		if( is_black( x, Y, z)) findhole( x, Y, z);
		else if( is_bkgnd( x,Y,z)) h->part_of_bkgnd = 1;
	Y = y + 1;
	if( Y < Im->height)
		if( is_black( x, Y, z)) findhole( x, Y, z);
		else if( is_bkgnd( x,Y,z)) h->part_of_bkgnd = 1;
	Z = z - 1;
	if( Z >= 0)
		if( is_black( x, y, Z)) findhole( x, y, Z);
		else if( is_bkgnd( x,y,Z)) h->part_of_bkgnd = 1;
	Z = z + 1;
	if( Z < Im->depth)
		if( is_black( x, y, Z)) findhole( x, y, Z);
		else if( is_bkgnd( x,y,Z)) h->part_of_bkgnd = 1;
}

struct OBJECT *o;
void findobj( COORD_type x, COORD_type y, COORD_type z)
     /*/ Recursively finds object at given voxel (x,y,z) */
{
	COORD_type X, Y, Z;

	TotalSizeOfObjects++;
	add_voxel_to_object( o, x, y, z);
	/*/ make sure the current voxel wont be done again */
	set_voxel( x, y, z, THRESHOLD);
	/*/ Scan 26-neighborhood for connected voxels */
	for( Z = z-1; Z <= z+1; Z++) {
	        if( Z < 0 || Z >= Im->depth) continue;
		for( X = x-1; X <= x+1; X++) {
		  if( X < 0 || X >= Im->width) continue; /*/ avoid going over edge */
			for( Y = y-1; Y <= y+1; Y++) {
			        if( Y < 0 || Y >= Im->height) continue;
				if( is_white( X, Y, Z)) {
				    findobj( X, Y, Z);
				} else if( findholes) {
				    if( is_black( X, Y, Z)) {
					h = add_hole_to_object( o);
					findhole( X, Y, Z);
					TotalHoles++;
				    }
				}
			}
		}
	}
}

struct OBJECT *find_objects( void)
{
	COORD_type x, y, z;
	struct OBJECT *pLastObject, *pDummyHead, *pNewObject, *pHead;
	/*/ Create dummy list head */
	pDummyHead = new_object();
	pLastObject = pDummyHead;

	/*/ Scan the whole image for white voxels */
	for( z = 0; z < Im->depth; z++) {
	    for( x = 0; x < Im->width; x++) {
	        for( y = 0; y < Im->height; y++) {
		    if( is_white( x, y, z)) {
		        pNewObject = new_object();
			/*/ find the object that contains the given voxel (x,y,z) */
			o = pNewObject;
			findobj( x, y, z);
			if( pNewObject->size >= min_obj_size) {
			  /*/ Add object to list only if it meets minimum size */
			    TotalObjects++;
			    /*/ Link new member to the end of the list */
			    pLastObject->pNext = pNewObject;
			    pLastObject = pNewObject;
			} else {
			  /*/ Otherwise let the new object go */
			    free_objects( pNewObject);
			}
		    }
		}
	    }
	}
	/*/ Drop the dummy head */
	pHead = pDummyHead->pNext;
	free( pDummyHead);
	/*/ Return the head of the object list */
	return pHead;
}

/*/////////////////////////////////////////////////////////////////////
// MISCELLANEOUS FUNCTIONS */

void remove_background( struct OBJECT *pFirstObject) 
{ /*/ Removes all the "bogus" holes that are actualy part of the background */
        struct HOLE *pHole, *pPrevHole, *pBogusHole;
	struct OBJECT *pThisObj;
	pThisObj = pFirstObject;
	while( pThisObj != NULL) {
	      pHole = pPrevHole = pThisObj->pHoleList;
	      /*/ Search for holes labeled as background */
	      while( pHole != NULL) {
		  if( pHole->part_of_bkgnd) {
		    /*/ Remove any such hole */
		    pThisObj->n_holes--;
		    pBogusHole = pHole;
		    if( pHole == pThisObj->pHoleList) {/*/ If it's the first object*/
		      pThisObj->pHoleList = pHole->pNext;
		      pPrevHole = pHole->pNext;
		    } else
		      pPrevHole->pNext = pHole->pNext;
		    pHole = pHole->pNext;
		    free_voxels( pBogusHole->pVoxels);
		    free( pBogusHole);
		  } else {
		    pPrevHole = pHole;
		    pHole = pHole->pNext;
		  }
	      }
	      pThisObj = pThisObj->pNext;
	}
}

void ensure_binary( struct IMAGE *Im)
{
	COORD_type x, y, z;
	for( z = 0; z < Im->depth; z++)
	        for(  x = 0; x < Im->width; x++)
		        for( y = 0; y < Im->height; y++)
				if( Im->p[idx(x,y,z)] > WHITE)
					set_voxel( x, y, z, WHITE);
}

void mark_background_pixels( struct IMAGE *Im)
{
	COORD_type x, y, z;
	int straight_line;
	/*/ Mark the contiguous below-threshold pixels starting from */
	  /*/ each side of the image*/
	for( z = 0; z < Im->depth; z++) {
	  for( x = 0; x < Im->width; x++) {
	    straight_line = 1;
	    /*/ From north*/
	    for( y = 0; y < Im->height; y++) {
	      if( is_black(x,y,z))
		set_voxel( x, y, z, BKGND);
	      else {
		straight_line = 0;
		break;
	      }
	    }
	    if( !straight_line) {
	      /*/ From south*/
	      for( y = Im->height-1; y >= 0; y--) {
		if( is_black(x,y,z))
		  set_voxel( x, y, z, BKGND);
		else
		  break;
	      }
	    }
	  }
	}
	for( z = 0; z < Im->depth; z++) {
	  for( y = 0; y < Im->height; y++) {
	    straight_line = 1;
	    /*/ From west*/
	    for( x = 0; x < Im->width; x++) {
	      if( is_black(x,y,z))
		set_voxel( x, y, z, BKGND);
	      else {
		straight_line = 0;
		break;
	      }
	    }
	    if( !straight_line) {
	      /*/ From east*/
	      for( x = Im->width-1; x >= 0; x--) {
		if( is_black(x,y,z))
		  set_voxel( x, y, z, BKGND);
		else
		  break;
	      }
	    }
	  }
	}
	if( Im->depth > 1) {
	  for( x = 0; x < Im->width; x++) {
	    for( y = 0; y < Im->height; y++) {
	      straight_line = 1;
	      /*/ From top*/
	      for( z = 0; z < Im->depth; z++) {
		if( is_black(x,y,z))
		  set_voxel( x, y, z, BKGND);
		else {
		  straight_line = 0;
		  break;
		}
	      }
	      /*/ From bottom*/
	      if( !straight_line) {
		for( z = Im->depth-1; z >= 0; z--) {
		  if( is_black(x,y,z))
		    set_voxel( x, y, z, BKGND);
		  else
		    break;
		}
	      }
	    }
	  }
	}
}

/*/////////////////////////////////////////////////////////////////////*/
#ifdef _STANDALONE_

#define HEADER_SIZE 16

struct IMAGE loadimage( FILE *f)
{
	struct IMAGE i;
	unsigned long file_length;
	unsigned long n_voxels;

	/*/ Read Header
	fread( &i, HEADER_SIZE, 1, f); /*/ Read width, height and depth*/
	if( ferror( f)){ printf( "Error reading image file!\n"); EXIT;}
	n_voxels = i.width*i.height*i.depth;
	fseek( f, 0, SEEK_END);
	file_length = ftell( f);
	if( HEADER_SIZE + n_voxels*i.bytes_per_voxel != file_length) {
		printf("Image file corrupted!\n"); EXIT; }
	/*/printf( "Image is %ix%ix%i, %i byte%s per voxel.\n",i.width,i.height,i.depth,i.bytes_per_voxel,(i.bytes_per_voxel>1)?"s":""); */
	/*/ Allocate memory for image */
	i.p = malloc( n_voxels*i.bytes_per_voxel);
	if( i.p == NULL){ printf( "Cannot allocate memory for image!\n"); EXIT;}
	/*/printf("n_voxels = %i\n",n_voxels);*/
	/*/ Read in data*/
	fseek( f, HEADER_SIZE, SEEK_SET);
	fread( i.p, i.bytes_per_voxel, n_voxels, f);
	if( ferror( f)){ printf( "Error reading image file!\n"); EXIT;}
	return i;
}

/*/#define OUTPUT_FILE "/tmp/ml_3dobjfind.out"*/
int output_objects_to_file( struct OBJECT *pObj, char *filename)
{
	struct OBJECT *pThis;
	struct VOXEL *pThisVoxel;
	FILE *f;
	int items_written;
	unsigned long n_objects_written;
	unsigned long n_voxels_written;
	unsigned long buffer;
	struct {
	  unsigned short y;
	  unsigned short x;
	  unsigned short z;
	} voxel_buffer;
	
	/*/ Open the file*/
	f = fopen( filename, "wb");
	if( f == NULL){ printf( "Could not open output file!\n"); EXIT;}
	/*/ Output the number of objects*/
	buffer = TotalObjects;
	items_written = fwrite( &buffer, 4, 1, f);
	if( items_written != 1){ printf( "Unable to write to output file!\n"); EXIT;}
	/*/ Output info about each object*/
	pThis = pObj;
	n_objects_written = 0;
	while( pThis != NULL) {
	  /*/ output object size*/
		buffer = pThis->size;
		items_written = fwrite( &buffer, 4, 1, f);
		if( items_written != 1){ printf( "Unable to write to output file!\n"); EXIT;}
		/*/ output the number of holes*/
		buffer = pThis->n_holes;
		items_written = fwrite( &buffer, 4, 1, f);
		if( items_written != 1){ printf( "Unable to write to output file!\n"); EXIT;}
		/*/ output the list of voxels*/
		pThisVoxel = pThis->pVoxels;
		n_voxels_written = 0;
		while( pThisVoxel != NULL) {
		  voxel_buffer.y = pThisVoxel->y + 1;
		  voxel_buffer.x = pThisVoxel->x + 1;
		  voxel_buffer.z = pThisVoxel->z + 1;
		  items_written = fwrite( &voxel_buffer, sizeof(voxel_buffer), 1, f);
		  if( items_written != 1){ printf( "Unable to write to output file!\n"); EXIT;}
		  pThisVoxel = pThisVoxel->pNext;
		  n_voxels_written++;
		}
		/*/ Check that the correct number of voxels was written*/
		if( pThis->size != n_voxels_written){
                    printf("Not all voxels could be written to output file!\n");
                    EXIT;
                }
		n_objects_written++;
		pThis = pThis->pNext;
	}
	/*/ Check that the correct number of objects was written*/
	if( TotalObjects != n_objects_written){
            printf("Not all objects could be written to output file!\n");
            EXIT;
        }
}

#define MIN_ARGS 3
#define MAX_ARGS 4
void print_usage( int argc, char* argv[])
{
	printf( "\n  Usage: %s <imagefile> <outputfile> <findholes> <min_obj_size>\n\n",argv[0]);
        printf("  The <findholes> argument must be 0 or 1\n");
}

int main( int argc, char* argv[])
{

    FILE *imagefile;
    struct IMAGE Image;
    struct OBJECT *pObjects;
    struct rlimit rlim;
    int stack_limit = 0; /*/ Stack limit in MegaBytes*/
    char *findholes_arg;
    char *min_objsize_arg;
    
    /*/ Initialize global variables*/
    TotalObjects = 0;
    TotalSizeOfObjects = 0;
    TotalHoles = 0;
    TotalSizeOfHoles = 0;
    findholes = 0;
    min_obj_size = 1;

    /*/ Process command line arguments*/
    if( argc < MIN_ARGS+1 || argc > MAX_ARGS+1){ print_usage(argc,argv); EXIT;}
    findholes_arg = argv[3];
    findholes = (int) findholes_arg[0] - '0';
    if( strlen(findholes_arg)!=1 || (findholes >> 1)){
        print_usage(argc,argv);
        EXIT;
    }
    /*/ If minimum object size specified, then get it*/
    if( argc > MIN_ARGS+1) {
        min_objsize_arg = argv[4];
	min_obj_size = (int) atoi( min_objsize_arg);
	if( min_obj_size < 1){
	    print_usage(argc,argv);
	    printf("Minimum object size cannot be smaller than 1");
	    EXIT;
	}
    }

    /*/ Load the image*/
    imagefile = fopen( argv[1], "rb");
    if( imagefile == NULL){ printf( "Cannot open file: %s\n", argv[1]); EXIT;}
    Image = loadimage( imagefile);
    fclose( imagefile);
    /*/ Initialize image attributes*/
    Image.threshold = 1;
    Image.COF.x = 0.0;
    Image.COF.y = 0.0;
    Image.COF.z = 0.0;

    /*/printf("(x y z) = (%i %i %i)\n", Image.width, Image.height, Image.depth);*/
    /*/printf("Finished loading image\n");*/

    /*/ Increase stack size*/
    getrlimit( RLIMIT_STACK, &rlim);
    rlim.rlim_cur = STACK_MEGABYTES*MEGABYTE;
    setrlimit( RLIMIT_STACK, &rlim);
    getrlimit( RLIMIT_STACK, &rlim);
    stack_limit = rlim.rlim_cur/MEGABYTE;
    if( stack_limit < STACK_MEGABYTES) {
      printf( "Unable to set stack limit to %i MB.\n", STACK_MEGABYTES);
      printf( "Setting Stack limit to the maximum available %i MB\n", stack_limit);
    }

    /*/ Initialize global Image pointer*/
    Im = &Image;
    /*/ Make sure image is binary (0-s and 1-s)*/
    ensure_binary( Im);
    /*/ Then mark background at the edges of image*/
    if( findholes)
      mark_background_pixels( Im);
    /*/ Find objects*/
    pObjects = find_objects();

    /*/ Remove the background*/
    if( findholes)
      remove_background( pObjects);

    /*/ Output list of objects to a file*/
	/*/printf("Total Objects = %i\n",TotalObjects);*/
    output_objects_to_file( pObjects, argv[2]);

    /*/ release memory*/
    free_objects( pObjects);
    free_image( Im);
    return 0;
}

#endif /*/ #ifdef _STANDALONE_*/

/*/////////////////////////////////////////////////////////////////////*/
#ifdef _MEX_

mxArray *make_voxellist( struct VOXEL *pVoxelList, long n_voxels)
{
  struct VOXEL *pThis;
  mxArray *newObject;
  int NDims = 2;
  int Dims[2];
  COORD_type *data;
  long i, idx;
  int x, y, z;

        Dims[0] = 3;
        Dims[1] = n_voxels;
        newObject = mxCreateNumericArray(NDims, Dims, mxUINT16_CLASS, mxREAL);
	data = (COORD_type*) mxGetData( newObject);
	pThis = pVoxelList;
	i = 0;
	while( pThis != NULL) {
	    data[i] = pThis->y + 1;
	    data[i+1] = pThis->x + 1;
	    data[i+2] = pThis->z + 1;
	    pThis = pThis->pNext;
	    i = i + 3;
	}
	return newObject;
}

#define N_HOLE_FIELDS 2
mxArray *make_holelist( struct HOLE *pHoleList, long n_holes)
{
  struct HOLE *pThis;
  int nfields = N_HOLE_FIELDS;
  const char *field_names[N_HOLE_FIELDS] = {"size","voxels"};
  mxArray *newHoleList, *newHole, *voxelList;
  mxArray *Size;
  double *buffer;
  long n = 0;

  newHoleList = mxCreateCellMatrix( 1, n_holes);
  pThis = pHoleList;
  while( pThis != NULL) {
    mxSetN( newHoleList, ++n); /*/ append cell array*/
      /*/ Create new object structure for the hole*/
      newHole = mxCreateStructMatrix( 1, 1, nfields, field_names);
      /*/ fill in the structure*/
	/*/ set hole size*/
      Size = mxCreateDoubleMatrix( 1, 1, mxREAL);
      buffer = mxGetPr( Size);
      *buffer = (double) pThis->size;
      mxSetField( newHole, 0, "size", Size);
      /*/ list of voxels constituting the hole*/
      voxelList = make_voxellist( pThis->pVoxels, pThis->size);
      mxSetField( newHole, 0, "voxels", voxelList);
      /*/ Put the new hole struct in the list*/
      mxSetCell( newHoleList, n-1, newHole);
      /*/ Get next hole*/
      pThis = pThis->pNext;
  }
  return newHoleList;
}

#define NFIELDS 4
void output_objects( struct OBJECT *pObj, mxArray *output)
{
	struct OBJECT *pThis;
	int nfields = NFIELDS;
	const char *field_names[NFIELDS] = {"size","voxels","n_holes","holes"};
	mxArray *newObject;
	long n_objects;
	mxArray *Size, *NHoles;
	double *buffer;
	mxArray *voxelList, *holeList;
	
	pThis = pObj;
	n_objects = 0;
	while( pThis != NULL) {
	  mxSetM( output, ++n_objects); /*/ append cell array*/
	    /*/ Create new object structure*/
		newObject = mxCreateStructMatrix( 1, 1, nfields, field_names);
		/*/ fill in the structure*/
		  /*/ output object size*/
		Size = mxCreateDoubleMatrix( 1, 1, mxREAL);
		buffer = mxGetPr( Size);
		*buffer = (double) pThis->size;
		mxSetField( newObject, 0, "size", Size);
		/*/ output the list of voxels*/
		voxelList = make_voxellist( pThis->pVoxels, pThis->size);
		mxSetField( newObject, 0, "voxels", voxelList);
		/*/ Output the number of holes*/
		NHoles = mxCreateDoubleMatrix( 1, 1, mxREAL);
		buffer = mxGetPr( NHoles);
		*buffer = (double) pThis->n_holes;
		mxSetField( newObject, 0, "n_holes", NHoles);
		/*/ Output the list of holes*/
		holeList = make_holelist( pThis->pHoleList, pThis->n_holes);
		mxSetField( newObject, 0, "holes", holeList);
		/*/ put the object in the output cell array*/
		mxSetCell( output, n_objects-1, newObject);
		pThis = pThis->pNext;
	}
}

void mexFunction(
		 int nlhs,
		 mxArray *plhs[], 
		 int nrhs,
		 const mxArray *prhs[])
{
  char* buf ;                 /*/ buffer for string value*/
    int buflen ;                /*/ Lentgh of string buffer*/
      int status ;                /*/ Storage for return values*/
	mxArray *output ;           /*/ Return value to Matlab*/

    struct IMAGE Image;
    struct OBJECT *pObjects;
    struct rlimit rlim;
    double* data; /*/ Data that matlab provides us*/
      int i; /*/ temporary counter*/
    mxArray *mxcell;
    int x, y, z;
    long n_voxels; /*/ The total number of voxels in the image*/
    long NDims;
    const int *dims;
    int Dims[3];
    unsigned char *img;
    int stack_limit = 0; /*/ Stack limit in MegaBytes*/
      /*/struct HOLE *pBackground;*/
    min_obj_size = 1;


    if (nrhs < 2 || nrhs > 3) {
        mexErrMsgTxt("ml_3dfindobj() requires 2 or 3 input arguments.");
    }

    /*/ Initialize global variables*/
    TotalObjects = 0;
    TotalSizeOfObjects = 0;
    TotalHoles = 0;
    TotalSizeOfHoles = 0;

    /*/ Get the option argument to determine if holes have to be found*/
    data = mxGetPr( prhs[1]);
    if( data == NULL)
	mexErrMsgTxt("Second argument to ml_3dfindobj() must be FINDHOLES");
    findholes = (int) (*data);

    /*/ Get the option argument to specify a minimum object size*/
    if( nrhs > 2) {
        data = mxGetPr( prhs[2]);
	if( data == NULL)
	    mexErrMsgTxt("Third argument to ml_3dfindobj() must be MIN_OBJ_SIZE");
	min_obj_size = (int) (*data);
	if( min_obj_size < 1)
	    mexErrMsgTxt("MIN_OBJ_SIZE must be at least one");
    }

    /*/printf("Starting to load image\n");*/
    /*/ Get Image size info*/
    NDims = mxGetNumberOfDimensions( prhs[0]);
    if( NDims > 3) mexErrMsgTxt("IMAGE should not have more than 3 dimensions");
    dims = mxGetDimensions( prhs[0]);
    Dims[0] = dims[0];
    Dims[1] = dims[1];
    if( NDims < 3) {
        Dims[2] = 1;
        NDims = 3;
    } else {
        Dims[2] = dims[2];
    }
    n_voxels = Dims[0] * Dims[1] * Dims[2]; 

    /*/ Get a hold of the input image*/
    if( !mxIsUint8( prhs[0])) 
        mexErrMsgTxt("Input IMAGE must be of class uint8");
    img = (unsigned char*) mxGetData( prhs[0]);
    Image.bytes_per_voxel = 1;
    Image.width = Dims[1];
    Image.height = Dims[0];
    Image.depth = Dims[2];
    Image.threshold = 1;

    /*/ Make a copy of the image so we do not have to modify the original*/
    Image.p = (unsigned char*) malloc( n_voxels);
    if( Image.p == NULL) mexErrMsgTxt("Could not allocate memory for temporary image");
    memcpy( Image.p, img, n_voxels);

    /*/printf("(x y z) = (%i %i %i)\n", Image.width, Image.height, Image.depth);*/


    /*/printf("Finished loading image\n");*/

    /*/ Increase stack size	getrlimit( RLIMIT_STACK, &rlim);*/
    rlim.rlim_cur = STACK_MEGABYTES*MEGABYTE;
    setrlimit( RLIMIT_STACK, &rlim);
    getrlimit( RLIMIT_STACK, &rlim);
    stack_limit = rlim.rlim_cur/MEGABYTE;
    if( stack_limit < STACK_MEGABYTES) {
        printf( "Unable to set stack limit to %i MB.\n", STACK_MEGABYTES);
        printf( "Setting Stack limit to the maximum available %i MB\n", stack_limit);
    }

    /*/ Initialize global Image pointer*/
    Im = &Image;
    /*/ Make sure image is binary (0-s and 1-s)*/
    ensure_binary( Im);
    /*/ Then mark background at the edges of image*/
    if( findholes)
        mark_background_pixels( Im);
    /*/ Find objects*/
    pObjects = find_objects();

    /*/ Remove the background*/
    if( findholes)
        remove_background( pObjects);

    /*/ OUTPUT LIST OF OBJECTS TO  MATLAB*/
	  /*/printf("Total Objects = %i\n",TotalObjects);*/
    plhs[0] = mxCreateCellMatrix( TotalObjects, 1) ;
    output = plhs[0];
    output_objects( pObjects, output);

    /*/ release memory*/
    free_objects( pObjects);
    free_image( Im);
    return;
}

#endif /*/ #ifdef _MEX_*/
