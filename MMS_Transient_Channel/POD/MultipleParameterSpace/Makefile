EXTRA_OBJS = ../../../DenseMatrixMatrixProduct/PlayaDenseMatrixMatrixProduct.o \
	     ../../../SourceFiles/POD_SVD.o \
	     ../../../SourceFiles/denseSerialMatrixIO.o \
	     ../../../SourceFiles/ODERHSBase.o \
	     ../../../SourceFiles/QuadraticODERHSBase.o \
	     ../../../SourceFiles/NSEProjectedODE.o \
	     ../../../SourceFiles/meshAddOns.o 



CLIENT_EXTRA_INCLUDES = -I../../../SourceFiles -I../../../DenseMatrixMatrixProduct $(SITE_EXTRA_INCLUDES)

CLIENT_EXTRA_LIBS = 

include ../../../Makefile.site
include ../../../Makefile.common

all: $(EXTRA_OBJS)
	make -j4 POD_Generator.exe
	make -j4 stage_KKT_required_files.exe

