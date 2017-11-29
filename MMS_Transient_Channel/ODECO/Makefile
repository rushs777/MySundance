include ../../Makefile.site
include ../../Makefile.common

EXTRA_OBJS =  ../../DenseMatrixMatrixProduct/PlayaDenseMatrixMatrixProduct.o\
	integralOperator.o \
	     KKTBase.o \
	     KKT_Transient_Channel.o\
	     sensorData.o \
	     ../../SourceFiles/denseSerialMatrixIO.o


# Absolute Path
#CLIENT_EXTRA_INCLUDES = -I/home/sirush/PhDResearch/SourceFiles -I/home/sirush/PhDResearch/MMS_Transient_Channel/ODECO
# Relative Path
CLIENT_EXTRA_INCLUDES = -I../../SourceFiles -I. -I../../DenseMatrixMatrixProduct $(SITE_EXTRA_INCLUDES)

CLIENT_EXTRA_LIBS = 


all: $(EXTRA_OBJS)
	make ODECO.exe
