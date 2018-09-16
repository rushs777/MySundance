# ------------------------------------------------------------------------
#
# Makefile.client.Sundance is an example of a simple makefile that will let
# a user build a Sundance-based application in an arbitrary directory.
# It will be configured by cmake to refer to the Trilinos installation
# directory. If you want to use the same client makefile with several
# Trilinos installations (debug and optimized, for instance) you can do 
# so: you'll merely need to hand-edit the TRILINOS_INSTALL_DIR that has been
# set automatically during installation.
#
# Author: Kevin Long, Texas Tech University
#         kevin.long@ttu.edu
# 
# To use this makefile, you can to do three steps:
#
# (1) Create a new Makefile that includes the installed Makefile.client.
#
# (2) In your Makefile, add the definition of EXTRA_OBJS to list names of
#     any object files required by your application other than that corresponding
#     to the application name.
#
# (3) (optional) Define CLIENT_EXTRA_LIBS and CLIENT_EXTRA_INCLUDES to
#     if you need to augment the header and/or library search path and the
#     list of libraries. 
#
# Example: Suppose you are building an application "RANS.exe" whose code is in 
# three source files: RANS.cpp, TurbulenceModel.cpp, and BC.cpp. You have
# installed Trilinos in /home/bob/Soft/Trilinos/INSTALL. You are building your 
# application in /home/bob/FlowModels/RANS. You would do the following:
#
# (1) Create a new makefile /home/bob/FlowModels/RANS/Makefile and include
#     /home/bob/Soft/Trilinos/INSTALL/Makefile.client.<PACKAGE> 
#
# (2) Fill in /home/bob/FlowModels/RANS/Makefile to set 
#     EXTRA_OBJS = TurbulenceModel.o BC.o
#
# (3) Do "make RANS.exe" to build the executable.
#
# If your usage is more complex (see comments below), then you may want to
# just copy this installed makefile to your own directory and then make the
# needed modifications.  However, if you do this you might have to make
# adjustments with new versions of Trilinos in case there are changes that
# need to be made to the export makefile system.  However, such changes are
# unlikely.
# 
# ------------------------------------------------------------------------


# ***
# *** Copy this below section into your own Makefile that includes this
# *** installed makefile Makefile.client.Sundance
# ***


#--------------------------------------------------------------------------
# If your application contains more than one source file, list any additional
# object files in the variable EXTRA_OBJS, e.g.,
# 
# EXTRA_OBJS = SomeFile.o SomeOtherFile.o
#
# In this section you can also add any dependency specifications for these
# extra files, e.g.,
# 
# SomeFile.o: SomeFile.hpp
# SomeOtherFile.o: SomeOtherFile.hpp SomeFile.hpp
# 
#--------------------------------------------------------------------------
# Sandbox Makefile

# Turn this one on for meshBoundaries.cpp
#EXTRA_OBJS = /home/sirush/PhDResearch/SourceFiles/meshAddOns.o

# Turn this one on for expandingPOD.cpp
#EXTRA_OBJS = /home/sirush/PhDResearch/SourceFiles/POD_SVD.o

# Turn this one on for ODE_POD_SVD.cpp
#EXTRA_OBJS = /home/sirush/PhDResearch/ODEOptimization/MoreAdvanced/integralOperator.o \
             /home/sirush/PhDResearch/ODEOptimization/MoreAdvanced/KKTBase.o \
             /home/sirush/PhDResearch/ODEOptimization/MoreAdvanced/KKT_Transient_Channel.o\
             /home/sirush/PhDResearch/ODEOptimization/MoreAdvanced/sensorData.o \
             /home/sirush/PhDResearch/SourceFiles/denseSerialMatrixIO.o \
	     /home/sirush/PhDResearch/SourceFiles/POD_SVD.o

# Turn this one for velocityROMRework.cpp, vstar_list_expr, exeToBeCalled
EXTRA_OBJS =

# Turn this on for snapshotMatrixBuilder.cpp, LinearOperatorMemoryStorage
#EXTRA_OBJS = /home/sirush/PhDResearch/SourceFiles/denseSerialMatrixIO.o


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
#TRILINOS_INSTALL_DIR = ${HOME}/Code/BUILDS/OPT
#VIENTO_INSTALL_DIR = ${HOME}/Projects/CFD/BUILDS/OPT

# Simon's system
#TRILINOS_INSTALL_DIR = /usr/local/Trilinos/BUILDS/SERIAL-OPT
VIENTO_INSTALL_DIR = ${HOME}/PhDResearch/BUILDS/Viento

#--------------------------------------------------------------------------
# If needed, set extra include and lib macros
#--------------------------------------------------------------------------
CLIENT_EXTRA_INCLUDES = -I/home/sirush/PhDResearch/SourceFiles -I/home/sirush/PhDResearch/Sandbox #-I/home/sirush/PhDResearch/ODEOptimization/MoreAdvanced



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
	$(CXX) -o $@ $(LDFLAGS) $(CXXFLAGS) $< $(EXTRA_OBJS) -L$(LIB_PATH) $(LIBS)

%.o : %.cpp
	$(CXX) -c -o $@ $(CXXFLAGS) -I$(INCLUDE_PATH) $(Viento_TPL_INCLUDES) $<

all:

clean:
	rm -f *.o

spotless:
	rm -f *.o *.exe

