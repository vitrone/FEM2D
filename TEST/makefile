#=============================================================================#
#                                                                             #
# Author: Vishal Vaibhav                                                      #
#                                                                             #
# Description: makefile for creating test targets                             #
#                                                                             #
# History:                                                                    #
# 22 Feb 2014 : created                                                       #
#                                                                             #
#=============================================================================#
SHELL = /bin/sh

CC = gcc
#
# MKL options
MKL_PATH   = "$(MKLROOT)/lib/intel64"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/intel64"
EXT = so

IFACE_LIB     = $(MKL_PATH)/libmkl_intel_ilp64.$(EXT)
THREADING_LIB = $(MKL_PATH)/libmkl_intel_thread.$(EXT)
OMP_LIB       = -L$(CMPLR_PATH) -liomp5
CORE_LIB      = $(MKL_PATH)/libmkl_core.$(EXT)

MKL_LIBS = -Wl,--start-group $(IFACE_LIB) $(THREADING_LIB) $(CORE_LIB) -Wl,--end-group $(OMP_LIB)

INCLUDES = -I../INCLUDE  -I$(MKLROOT)/include

# CFLAG is used for the implicit rule to generate object files for C
# $(CC) -c $(CFLAGS)
# use -save-temps to see the preprocess effects in .i files
#
# Machine dependent options
MACH_DEP_OPT = -march=native
# Optimization options, -Ofast enables all -03 level options 
OPTIMIZE = -Ofast -funroll-all-loops

CFLAGS = $(INCLUDES) -g -ansi -D_GNU_SOURCE -fexceptions -fPIC  \
	 $(MACH_DEP_OPT) $(OPTIMIZE)                          \
         -fno-omit-frame-pointer -std=c99 -DMKL_ILP64 -m64

LDFLAGS = -L../BUILD -L$(MKL_PATH) $(MKL_LIBS) -lpthread -lm -lcunit -lfem2d -lmkl_rt

SOURCES = CUnit_fem2d.c   \
	  CUnit_matrix.c  \
	  matlib_io_print.c

TSTBIN = BIN
#OBJECTS = $(SOURCES:%.c=$(TSTBIN)/%.o)
TARGETS = $(SOURCES:%.c=$(TSTBIN)/%)
FEMLIB = ../BUILD/libfem2d.so


.PHONY: bin_dir all clean


all : bin_dir $(TARGETS) 

bin_dir:
	@mkdir -p $(TSTBIN)

$(TSTBIN)/% : $(TSTBIN)/%.o $(FEMLIB)
	@echo "Building binary files..."
	$(CC) $(CFLAGS) $(LDFLAGS) $(CLIBS) $< -o $@ 

$(TSTBIN)/%.o : %.c
	@echo "Building object files..."
	$(CC) $(CFLAGS) -c $< -o $@ 


clean :
	@rm -f $(OBJECTS)

