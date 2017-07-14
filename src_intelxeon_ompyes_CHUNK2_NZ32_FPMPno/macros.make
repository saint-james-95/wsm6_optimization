# Beginning of script-generated settings
ARCH      = intelxeon
OPENMP    = yes
CHUNK     = 2
TRANSLATEI= auto
FPMP      = no
IINSIDE   = no
NZ        = 32
DEFINES   = 
CPP_FLAGS = 
DEFINES += -DALIGN_OK
CPP_FLAGS += -I/home/keirouz/Pearls2_Chapter02_Copy/src_intelxeon_ompyes_CHUNK2_NZ32_FPMPno/tools
CHUNKCMD = /home/keirouz/Pearls2_Chapter02_Copy/src_intelxeon_ompyes_CHUNK2_NZ32_FPMPno/tools/chunkfilter.pl
LITERALKCMD = /home/keirouz/Pearls2_Chapter02_Copy/src_intelxeon_ompyes_CHUNK2_NZ32_FPMPno/tools/literalkfilter.pl
# End of script-generated settings

# Beginning of macros.make.intelxeon
# For parallel make
GMAKEMINUSJ = -j8

# Compiler settings. Add build-specific cpp definitions here if needed.  
BUILD_DEFINES = 
LDFLAGS   = $(INSTRFLAGS) $(OMPFLAG)
CPP       = /lib/cpp
CPP_OPTS  = -traditional -P
CPP_FLAGS += $(CPP_OPTS) $(GPTLINCLUDE) $(DEFINES) $(BUILD_DEFINES)
FC        = ifort
FCserial  = ifort
FCfrontend = ifort
FCgfort   = gfortran
FCncarg   = ncargf90

# Flag to enable "-fp-model precise" ("makenim" automatically sets this to 
# yes or no)
ifeq ($(FPMP),yes)
  PRECISEFLAG = -fp-model precise
endif

# Fortran flags for various subdirs. If DEBUG=yes, set for debugging
# fpe0 will cause the code to abort on overflow, divide-by-zero, or invalid
# fpe0 does not change the output but adds about 15% to the run time
DEBUG = no
ifeq ($(DEBUG),yes)
  PHYSFLAGS  = -g -O0 -ftz -traceback -I ../include -fpe0 -fp-model precise
else
  PHYSFLAGS  = -g -O3 -ftz -traceback -I ../include -opt-report-phase=hlo -vec-report6 -align array64byte -xHost $(PRECISEFLAG)
endif

# Whether to enable auto-profiling
PROFILE = no
ifeq ($(PROFILE),yes)
  INSTRFLAGS = -finstrument-functions -rdynamic
  PHYSFLAGS += $(INSTRFLAGS)
endif

# Flag to target native mode on Xeon-Phi (empty for everyone except xeonphi)
MICFLAGS =

# Implicit real*8 flag needed in physics
R8FLAG = -r8

# Flag to enable OpenMP ("makenim" automatically sets this to yes or no)
ifeq ($(OPENMP),yes)
  OMPFLAG = -qopenmp
  BUILD_DEFINES += -D_OPENMP
endif

# Flag to switch to i-on-inside versions of nislfv_rain_plm*() routines
ifeq ($(IINSIDE),yes)
  BUILD_DEFINES += -DIINSIDE
endif

# GPTL timing library: To shut off entirely, set USE_GPTL = no
USE_GPTL = yes
ifeq ($(USE_GPTL),yes)
  GPTL = ../gptl_install
  GPTLINCLUDE = -I$(GPTL)/include
  GPTLLIB     = -L$(GPTL)/lib -lgptl
endif
# End of macros.make.intelxeon

# Beginning of macros.make.all
# USE_GPTL=yes means we'll link to the library rather than use dummy objects
# In this case GPTLINCLUDE was set in macros.make

RM = rm -f

ifeq ($(USE_GPTL),yes)
  GPTLOBJS =
else
  GPTLINCLUDE = -I$(SRCDIR)/dummygptl
  GPTLOBJS = $(shell ls $(SRCDIR)/dummygptl/*.o)
endif

FCX    = $(FCserial)

WSM6EXE = $(BINDIR)/wsm6kernel
# End of macros.make.all

