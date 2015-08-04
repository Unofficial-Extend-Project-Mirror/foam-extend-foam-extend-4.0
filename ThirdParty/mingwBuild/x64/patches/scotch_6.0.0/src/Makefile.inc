# This make include file is intended for building under MinGW32. As is,
# it relies on :
#   - pthread-win32 (http://sourceware.org/pthreads-win32/),
#   - zlib (http://www.zlib.net/)
#   - mpich2 (http://www.mcs.anl.gov/research/projects/mpich2/).
# It must be adapted to reflect your environment, in particular
# installation root directories, include path and library name.
# Since all files are linked with the MPI libraries, this file is
# for compiling PT-Scotch only. To compile Scotch, remove the
# -DSCOTCH_PTSCOTCH flag, as well as all references to MPI in the
# CFLAGS and CLIBFLAGS	=
#LDFLAGS variables.

#MPI_ROOTDIR		= /c/Programs/OpenMPI_v1.6.1-x64

EXE			=
LIB			= .a
OBJ			= .o

MAKE			= make
AR			= ar
ARFLAGS			= -ruv
CAT			= cat
CCS			= gcc
CCP			= mpicc
CCD			= mpicc -I$(MPI_ROOTDIR)/include

#--- Compiler/loader flags
CFLAGS_CPL		= -O0 -g3 --std=c99
CFLAGS_INC		=
CFLAGS_DEF		= -DCOMMON_RANDOM_FIXED_SEED -DCOMMON_STUB_FORK -DSCOTCH_PTSCOTCH -DSCOTCH_RENAME -D'pipe(pfds)=_pipe(pfds,1024,0x8000)' -DHAVE_STDINT_H=0 -DHAVE_UINT_T=1
CLIBFLAGS	=
LDFLAGS   		=

PGMFILES=$(PROGRAMFILES)

#--- MPI
CFLAGS_INC	+= -I$(MPI_ROOTDIR)/include
#--- Comment/Uncomment for threaded MPI
CLIBFLAGS	=
LDFLAGS    	+= -L$(MPI_ROOTDIR)/bin -lm -lmpi -lmpid
#CLIBFLAGS	=
#LDFLAGS    	+= -L$(MPI_ROOTDIR)/lib -lm -lmpich2mt

#--- Pthread : Uncomment for pthread support
#PTHREAD_ROOTDIR		= $(PGMFILES)/pthread-win32
#CFLAGS_INC		+= -I$(PTHREAD_ROOTDIR)/include
#CLIBFLAGS	=
LDFLAGS		+= -L$(PTHREADS_HOME)/Pre-built.2/lib/x64 -lpthreadGC2

#--- zlib: Uncomment for compressed files
#ZLIB_ROOTDIR		= $(PGMFILES)/zlib-1.2.3
#CFLAGS_INC	+= -I$(ZLIB_ROOTDIR)/include
#CLIBFLAGS	=
LDFLAGS	+= -lz

#--- COMMON_PTHREAD: Uncomment for compressed files
#CFLAGS_DEF	+= -DCOMMON_PTHREAD -DCOMMON_FILE_COMPRESS_GZ

#--- SCOTCH_PTHREAD: Uncomment for threaded MPI
#CFLAGS_DEF	+= -DSCOTCH_PTHREAD

CFLAGS		= $(CFLAGS_CPL) $(CFLAGS_INC) $(CFLAGS_DEF)
CP			= cp
LEX			= flex -Pscotchyy -olex.yy.c
LN			= cp
MKDIR			= mkdir
MV			= mv
RANLIB			= ranlib
YACC			= bison -pscotchyy -y -b y
