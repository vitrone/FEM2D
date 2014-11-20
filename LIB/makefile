#=============================================================================#
#                                                                             #
# Author: Vishal Vaibhav                                                      #
#                                                                             #
# Description: makefile for creating shared object libraries                  #
#                                                                             #
# History:                                                                    #
# 19 July 2013: created                                                       #
#                                                                             #
#=============================================================================#
SHELL = /bin/sh

# Basics of linking:
# (1) Static linking: a static library or statically-linked library is a set of 
# routines, external functions and variables which are resolved in a caller at 
# compile-time and copied into a target application by a compiler, linker, or 
# binder, producing an object file and a stand-alone executable.
# (2)

CC = gcc
# Default locations where gcc searches for header files 
# /usr/local/include
# /usr/include
INC_DIR = ../INCLUDE
BUILDDIR = ../BUILD

# Implicit Rules: Some rules are known to "make". They need not be stated
# explicitly, for example, how to create .o files from .c files.
# For compiling C programs, following implicit rule is used 
# $(CC) -c $(CFLAGS)
# Therefore, CFLAG must be defined correctly.

# # COMPILER OPTIONS: (Some may be GCC spepcific and needs to be signled out)
#
# -c <file list>: used to generate object files (default extension is '.o')
# -o <name>: used to specify name of the output
#
# -fPIC : generate position independent code (PIC)
#
# -fexceptions: enable exception handling
#
# -fno-omit-frame-pointer: 
#
# -shared: produce shared object
#
# -D_GNU_SOURCE: defines the macro _GNU_SOURCE which enables all GNU C extensions
#
# -Wl, <options>: passes <option> to the linker, comma separated options are
#  split into multiple options. 
#  example: '-Wl, -Map, output.map' passes '-Map output.map' to the linker
#
# # WARNINGS: 
# -Wall: enables all warnings
#
# # STATIC LINKING: (linking of .a files)
#  -lNAME: Used to link a library which has the name libNAME.a in the standard
#          library locations. GCC gives preference to shared objects over static
#          if both are present.
#
# EVIRONMENT VARIABLES:
#  -I : for specifying include path
#  -L : for library path
# Both of these above option can be configured by using environment variables 
#  C_INCLUDE_PATH, LIBRARY_PATH
#  Load path can be added by using LD_LIBRARY_PATH
#
# # TYPES OF LIBRARIES: 
# (1) Static type: extension '.a', requires static linking
# (2) Shared type: extension '.so', requires dynamic linking

EXT = so

#ifeq ($(SD),static)
#  EXT=a
#  RES_EXT=lib
#  S_GRP=-Wl,--start-group
#  E_GRP=-Wl,--end-group
#else
#  EXT=so
#  RES_EXT=so
#  S_GRP=
#  E_GRP=
#endif
#
#ifeq ($(EXT),so)
#    MKL_LD_PATH=-L"$(MKL_PATH)"
#    MKL_PREFIX=-l
#    MKL_SUFFIX=
#else
#    MKL_LD_PATH=
#    MKL_PREFIX="$(MKL_PATH)/lib
#    MKL_SUFFIX=.a"
#endif
#
# MKL options
MKL_PATH   = "$(MKLROOT)/lib/intel64"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/intel64"

# OpenMP library : OMP_LIB for parallelism in MKL
IFACE_LIB      = $(MKL_PATH)/libmkl_intel_ilp64.$(EXT)
THREADING_LIB  = $(MKL_PATH)/libmkl_intel_thread.$(EXT)
SEQUENTIAL_LIB = $(MKL_PATH)/libmkl_sequential.$(EXT)
OMP_LIB        = -L$(CMPLR_PATH) -liomp5
CORE_LIB       = $(MKL_PATH)/libmkl_core.$(EXT)

MKL_LIBS = -Wl,--start-group $(IFACE_LIB) $(THREADING_LIB) $(CORE_LIB) -Wl,--end-group $(OMP_LIB)
#MKL_LIBS = -Wl,--start-group $(IFACE_LIB) $(SEQUENTIAL_LIB) $(CORE_LIB) -Wl,--end-group

# Machine dependent options
MACH_DEP_OPT = -march=native
# Optimization options, -Ofast enables all -03 level options 
OPTIMIZE = -Ofast -funroll-all-loops

INCLUDES = -I$(INC_DIR) -I$(MKLROOT)/include

CFLAGS = $(INCLUDES) -g -ansi -D_GNU_SOURCE -fexceptions -fPIC   \
         -fno-omit-frame-pointer -std=c99                     \
	 $(MACH_DEP_OPT) $(OPTIMIZE) -DMKL_ILP64 -m64                        
          

LDFLAGS =  -shared -L$(MKL_PATH) $(MKL_LIBS) -lpthread -lm -lmkl_rt

# Default locations where gcc searches for libraries 
# /usr/local/lib
# /usr/lib

TARGET = $(BUILDDIR)/libfem2d.so

SOURCES = triangle.c      \
	  triangle_data.c \
	  triangle_mem.c  \
	  triangle_check.c\
	  triangle_debug.c\
	  fem2d.c         \
	  fem1d.c         \
	  matlib.c        \
	  matlib_solver.c \
	  matlib_io.c     \
	  jacobi.c        \
	  legendre.c      \
	  sym_quad.c

OBJECTS = $(SOURCES:%.c=$(BUILDDIR)/%.o)

HLIST = triangle.h \
	matlib.h   \
	basic.h    \
	debug.h    \
	ehandler.h \
	legendre.h \
	fem2d.h    \
	fem1d.h

HEADERS = $(addprefix $(INC_DIR)/, $(HLIST))

.PHONY: all clean

#.SECONDARY: $(OBJECTS)

all : $(TARGET) 

$(TARGET): $(OBJECTS) $(HEADERS)
	@echo "Creating shared object..."	
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -o $(TARGET)

$(BUILDDIR)/%.o: %.c
	@echo "Building object files..."
	$(CC) $(CFLAGS) -c $< -o $@ 	

clean :
	rm -f $(OBJECTS)

