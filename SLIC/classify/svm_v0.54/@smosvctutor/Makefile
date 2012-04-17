#
# uncomment for windows systems
# 
#mex=mex.bat

/smotrain.mexglx	: smosvctrain.cpp \
			  Cache.h     \
			  LrrCache.h  \
			  LrrCache.cpp \
			  InfCache.h  \
			  InfCache.cpp \
			  SmoTutor.h  \
		 	 SmoTutor.cpp
		$(mex) -g smosvctrain.cpp InfCache.cpp LrrCache.cpp SmoTutor.cpp
