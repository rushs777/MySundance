EXTRA_OBJS = ../../DenseMatrixMatrixProduct/PlayaDenseMatrixMatrixProduct.o \
             ../../SourceFiles/PlayaSVD.o \
	     ../../SourceFiles/denseSerialMatrixIO.o \
	     ../../SourceFiles/ODERHSBase.o \
	     ../../SourceFiles/QuadraticODERHSBase.o \
	     ../../SourceFiles/MyNLO.o \
	     ../../SourceFiles/velocityROM.o \
	     ../../SourceFiles/meshAddOns.o \
	     ../../SourceFiles/NSEProjectedODE.o

CLIENT_EXTRA_INCLUDES = -I../../SourceFiles  -I../../DenseMatrixMatrixProduct $(SITE_EXTRA_INCLUDES)

CLIENT_EXTRA_LIBS = 

include ../../Makefile.site
include ../../Makefile.common

all: $(EXTRA_OBJS)
	make uRO.exe
