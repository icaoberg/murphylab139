MEX='/usr/local/matlab/bin/mex'
CC='/usr/bin/gcc'

ml_Har_Temporal_Texture: ml_Har_Temporal_Texture.c ml_Extract_Temporal_Texture.o 
	${MEX} ml_Har_Temporal_Texture.c ml_Extract_Temporal_Texture.o

ml_Extract_Temporal_Texture.o: ml_Extract_Temporal_Texture.c
	gcc -c -IInclude -ansi ml_Extract_Temporal_Texture.c
