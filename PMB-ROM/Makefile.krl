# ------------------------------------------------------------------------

EXTRA_OBJS = ../SourceFiles/PlayaSVD.o PMB.o  ../DenseMatrixMatrixProduct/PlayaDenseMatrixMatrixProduct.o ../SourceFiles/ODERHSBase.cpp ../SourceFiles/QuadraticODERHSBase.o ../SourceFiles/denseSerialMatrixIO.o

#EXTRA_OBJS =

# ***
# *** End copy section
# ***



#--------------------------------------------------------------------------
# Set the value of TRILINOS_INSTALL_DIR to point to the directory
# into which you have installed Trilinos. For example, if you have
# installed into /usr/local/trilinos-3.14.59, this variable will be set to
#
# TRILINOS_INSTALL_DIR = /usr/local/trilinos-3.14.59 
#
# automatically by CMake. If you're using a single application makefile to
# build against several versions of Trilinos you may need to edit
# TRILINOS_INSTALL_DIR manually.
#
#--------------------------------------------------------------------------

#TRILINOS_INSTALL_DIR = /usr/local/trilinos-devel-opt

# Kevin's system
TRILINOS_INSTALL_DIR = ${HOME}/Code/BUILDS/OPT
VIENTO_INSTALL_DIR = ${HOME}/Projects/CFD/BUILDS/OPT
# Simon's system
#TRILINOS_INSTALL_DIR = /usr/local/Trilinos/BUILDS/SERIAL-OPT
#VIENTO_INSTALL_DIR = ${HOME}/PhDResearch/BUILDS/Viento


#--------------------------------------------------------------------------
# If needed, set extra include and lib macros
#--------------------------------------------------------------------------
#CLIENT_EXTRA_INCLUDES = -I/home/sirush/PhDResearch/SourceFiles
CLIENT_EXTRA_INCLUDES =-I../SourceFiles -I../DenseMatrixMatrixProduct 

CLIENT_EXTRA_LIBS = 








#--------------------------------------------------------------------------
# The remaining lines will not normally need to be changed. 
#
# Cases where you'll need to edit them include 
# (1) Your application needs some compiler flags not used in the Trilinos
#     build. Be careful to avoid inconsistency: for instance, mixing C++ 
#     libraries built with and without STL checking can cause segfaults.
# (2) Your application must link to some 3rd party libraries not among those 
#     specified when building Trilinos. You'll need to add them to the
#     linker command line. You may need to add their locations
#     to the library search path, and to the rpath if you're using shared
#     libraries.  
#--------------------------------------------------------------------------

# Include the Trilinos export makefile from package=Sundance.
#include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.Sundance
include $(VIENTO_INSTALL_DIR)/include/Makefile.export.Viento

# Add the Trilinos installation directory to the search paths
# for libraries and headers
LIB_PATH = $(TRILINOS_INSTALL_DIR)/lib $(Viento_LIBRARY_DIRS)

INCLUDE_PATH = $(TRILINOS_INSTALL_DIR)/include $(Viento_INCLUDE_DIRS) $(CLIENT_EXTRA_INCLUDES)

# Set the C++ compiler and flags to those specified in the export makefile
CXX = $(Viento_CXX_COMPILER)

CXXFLAGS = $(Viento_CXX_FLAGS)

# Add the Trilinos libraries, search path, and rpath to the 
# linker command line arguments 
LIBS = $(CLIENT_EXTRA_LIBS) $(Viento_SHARED_LIB_RPATH_COMMAND) \
$(Viento_LIBRARIES) \
$(Viento_TPL_LIBRARIES) $(Viento_EXTRA_LD_FLAGS) 

# Rules for building executables and objects. 
%.exe : %.o $(EXTRA_OBJS)
	$(CXX) -o $@ $(LDFLAGS) -I$(INCLUDE_PATH) $(Viento_TPL_INCLUDES) $(CXXFLAGS) $< $(EXTRA_OBJS) -L$(LIB_PATH) $(LIBS)

%.o : %.cpp
	$(CXX) -c -o $@ $(CXXFLAGS) -I$(INCLUDE_PATH) $(Viento_TPL_INCLUDES) $<

all:

clean:
	rm -f *.o

spotless:
	rm -f *.o *.exe

redo:
	make spotless;	\
	make testVMB.exe; \
	./testVMB.exe
