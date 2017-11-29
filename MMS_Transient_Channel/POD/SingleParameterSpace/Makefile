#

include ../../../Makefile.site
include ../../../Makefile.common

EXTRA_OBJS = ../../../DenseMatrixMatrixProduct/PlayaDenseMatrixMatrixProduct.o \
../../../SourceFiles/POD_SVD.o \
../../../SourceFiles/denseSerialMatrixIO.o \
	     ../../../SourceFiles/ODERHSBase.o \
	     ../../../SourceFiles/QuadraticODERHSBase.o \
	     ../../../SourceFiles/NSEProjectedODE.o \
	     ../../../SourceFiles/meshAddOns.o 



CLIENT_EXTRA_INCLUDES = -I../../../SourceFiles -I../../../DenseMatrixMatrixProduct $(SITE_EXTRA_INCLUDES)

CLIENT_EXTRA_LIBS = 

all: $(EXTRA_OBJS)
	make POD_Generator.exe
	make stage_KKT_required_files.exe
