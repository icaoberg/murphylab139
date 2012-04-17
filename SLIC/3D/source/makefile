MEX='/usr/local/matlab6/bin/mex'
CC='/usr/bin/gcc'

cvip_pgm3Dtexture.o : cvip_pgm3Dtexture.c
	gcc -c -IInclude -ansi cvip_pgm3Dtexture.c
ml_3dbgsub: ml_3dbgsub.c
	${MEX} -D_MEX_ ml_3dbgsub.c
ml_binarize: ml_binarize.c
	${MEX} -D_MEX_ ml_binarize.c
ml_3dfindobj: ml_3dfindobj.c
	${MEX} -D_MEX_ ml_3dfindobj.c
ml_3dfindobj_sa: ml_3dfindobj.c
	${CC} -D_STANDALONE_ -o ml_3dfindobj_sa ml_3dfindobj.c
ml_3Dtexture: ml_3Dtexture.c cvip_pgm3Dtexture.o
	${MEX}  -D_MEX_ ml_3Dtexture.c cvip_pgm3Dtexture.o
ml_3Dtexture_all: ml_3Dtexture_all.c
	${MEX} -D_MEX_ ml_3Dtexture_all.c
