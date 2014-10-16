#=============================================================================#
#                                                                             #
# Author: Vishal Vaibhav                                                      #
#                                                                             #
# Description: makefile for creating .mexa64 binaries                         #
#                                                                             #
#                                                                             #
#=============================================================================#
SHELL = /bin/sh

MEX = mex
CC = gcc
# MKL options
MKL_PATH   = "$(MKLROOT)/lib/intel64"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/intel64"
EXT = so

IFACE_LIB     = $(MKL_PATH)/libmkl_intel_lp64.$(EXT)
THREADING_LIB = $(MKL_PATH)/libmkl_intel_thread.$(EXT)
OMP_LIB       = -L$(CMPLR_PATH) -liomp5
CORE_LIB      = $(MKL_PATH)/libmkl_core.$(EXT)

MKL_LIBS = -Wl,--start-group $(IFACE_LIB) $(THREADING_LIB) $(CORE_LIB) -Wl,--end-group $(OMP_LIB)

INCLUDES = -I/usr/local/MATLAB/R2013a/extern/include   \
           -I/usr/local/MATLAB/R2013a/simulink/include \
           -I../INCLUDE -I$(MKLROOT)/include


# CFLAG is used for the implicit rule to generate object files for C
# $(CC) -c $(CFLAGS)
CFLAGS = -ansi -D_GNU_SOURCE -DMATLAB_MEX_FILE -fexceptions -fPIC   \
         -fno-omit-frame-pointer -std=c99                           \
         $(INCLUDES)  -DMKL_ILP64 -m64

LDFLAGS = -shared -L$(MKL_PATH) $(MKL_LIBS) -L../BUILD                                     \
          -Wl,--version-script,/usr/local/MATLAB/R2013a/extern/lib/glnxa64/mexFunction.map \
	  -Wl,-rpath-link,/usr/local/MATLAB/R2013a/bin/glnxa64                             \
          -L/usr/local/MATLAB/R2013a/bin/glnxa64                                           \
	  -lmx -lmex -lmat -lstdc++ -lpthread -lm -lfem2d -lmkl_rt

CLIBS =               

BUILDDIR = ../MEXBUILD
INC_DIR  = ../INCLUDE
LIB_DIR  = ../LIB
# add the .c mex files here
SOURCES = tricall.c          \
	  mxfem2d_delaunay.c \
	  mxfem2d_remesh.c   \
	  mxfem2d_checkmesh.c \
	  mxmatlib_gemm.c


TARGETS = $(SOURCES:%.c=$(BUILDDIR)/%.mexa64)
OBJECTS = $(SOURCES:%.c=$(BUILDDIR)/%.o)

LIB_DEP = ../BUILD/libfem2d.so

.PHONY: all clean
# Automatic Variables
# $@ - The file name of the target of the rule.
# $< - The name of the first prerequisite
#
# If .o files are created automatically using implicit rules then are assumed
# to be intermediate and are deleted automatically by make. Making it an 
# explicit prerequisite of a special target .SECONDARY prevents this.
# 
# .SECONDARY: $(OBJECTS)

all : $(TARGETS) 

# Pattern Rules: A pattern rule contains the character '%' 
# (exactly one of them) in the target; otherwise, it looks exactly like an 
# ordinary rule. The target is a pattern for matching file names; the '%' 
# matches any nonempty substring, while other characters match only 
# themselves. For example, '%.c' as a pattern matches any file name that 
# ends in '.c'. 's.%.c' as a pattern matches any file name that starts with
# 's.', ends in '.c' and is at least five characters long. 
# (There must be at least one character to match the '%'.) 
# The substring that the '%' matches is called the stem. 
#
# '%' in a prerequisite of a pattern rule stands for the same stem that 
# was matched by the '%' in the target. In order for the pattern rule to 
# apply, its target pattern must match the file name under consideration. 
# All of its prerequisites (after pattern substitution) must name 
# files that exist or can be made. These files become prerequisites of 
# the target.

#$(BUILDDIR)/%.mexa64: %.mexa64
#	mv $< $@

$(BUILDDIR)/%.mexa64 : $(BUILDDIR)/%.o $(LIB_DEP)
	@echo "Building mex-binary files..."
	$(CC) $(LDFLAGS) $< -o $@ $(CLIBS) 

$(BUILDDIR)/%.o : %.c
	@echo "Building object files..."
	$(CC) $(CFLAGS) -c $< -o $@ 

# Note: The following contruct gives error while running the MATLAB unittester.
# In fact it is incorrect because it lists all the elements in OBJECTS as
# prerequisite for all the elements in TARGETS. It is better to be explicit
# here.
# $(TARGETS): $(OBJECTS)
# 	$(CC) $(LDFLAGS) $< -o $(BUILDDIR)/$@ $(CLIBS) 
# 

clean :
	rm -f $(OBJECTS)

