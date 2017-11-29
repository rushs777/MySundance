

include ../../Makefile.site
include ../../Makefile.common

# Relative Path
EXTRA_OBJS = ../../DenseMatrixMatrixProduct/PlayaDenseMatrixMatrixProduct.o \
../../SourceFiles/PlayaSVD.o \
	     ../../SourceFiles/denseSerialMatrixIO.o \
	     ../../SourceFiles/ODERHSBase.o \
	     ../../SourceFiles/QuadraticODERHSBase.o \
	     ../../SourceFiles/MyNLO.o \
	     ../../SourceFiles/velocityROM.o \
	     ../../SourceFiles/meshAddOns.o \
	     ../../SourceFiles/NSEProjectedODE.o


# Relative Path
CLIENT_EXTRA_INCLUDES = -I../../SourceFiles  -I../../DenseMatrixMatrixProduct $(SITE_EXTRA_INCLUDES)

CLIENT_EXTRA_LIBS = 


all: $(EXTRA_OBJS)
	make uRO.exe
