EXTRA_OBJS =  ../../DenseMatrixMatrixProduct/PlayaDenseMatrixMatrixProduct.o\
	      integralOperator.o \
	      KKTBase.o \
	      KKT_Transient_Channel.o\
	      sensorData.o \
	      ../../SourceFiles/denseSerialMatrixIO.o

CLIENT_EXTRA_INCLUDES = -I../../SourceFiles -I. -I../../DenseMatrixMatrixProduct $(SITE_EXTRA_INCLUDES)

CLIENT_EXTRA_LIBS = 

include ../../Makefile.site
include ../../Makefile.common

all: $(EXTRA_OBJS)
	make ODECO.exe
