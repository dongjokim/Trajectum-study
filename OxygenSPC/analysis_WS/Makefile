PROGRAM         = JCorrSPCRun3

LIBDIRARCH      = lib

# Detect operating system
UNAME_S := $(shell uname -s)

# Set compiler and flags based on platform
ifeq ($(UNAME_S),Darwin)
    # macOS settings
    CXX             = clang++
    CXXFLAGS        = -std=c++20 -O2 -g -arch arm64 \
               -isysroot $(shell xcrun --show-sdk-path) \
               -mmacosx-version-min=14.0 \
               -Wno-unknown-warning-option \
               -fPIC \
               -stdlib=libc++
    LD              = clang++
else
    # Linux settings
    CXX             = g++
    CXXFLAGS        = -std=c++20 -O2 -g -fPIC
    LD              = g++
endif

SOFLAGS         = -shared -Wno-deprecated
CXXFLAGS        += $(shell root-config --cflags | sed 's/-std=c++17/-std=c++20/')
LIBS            = $(shell root-config --libs)

#INCS += $(shell $(FASTJET)/bin/fastjet-config --cxxflags )
#LIBS += $(shell $(FASTJET)/bin/fastjet-config --libs --plugins ) $(INCS)

# LIBS	       += -L$(PYTHIA8)/$(LIBDIRARCH) -lpythia8 -llhapdfdummy
# LIBS	       += -L$(PYTHIA8)/lib -lpythia8
# INCS           += -I$(PYTHIA8)/include

# Platform-specific O2 Framework paths
ifeq ($(UNAME_S),Darwin)
    # macOS O2 paths
    O2_ROOT ?= /Users/djkim/alice/sw/osx_arm64/O2/daily-20260420-0000-local1
    O2PHYSICS_ROOT ?= /Users/djkim/alice/sw/osx_arm64/O2Physics/master-local1
    FAIRLOGGER_ROOT ?= /Users/djkim/alice/sw/osx_arm64/FairLogger/latest
    FMT_ROOT ?= /Users/djkim/alice/sw/osx_arm64/fmt/latest
    FAIRMQ_ROOT ?= /Users/djkim/alice/sw/osx_arm64/FairMQ/latest
    BOOST_ROOT ?= /Users/djkim/alice/sw/osx_arm64/Boost/latest
    ARROW_ROOT ?= /Users/djkim/alice/sw/osx_arm64/arrow/latest
    MS_GSL_ROOT ?= /Users/djkim/alice/sw/osx_arm64/ms_gsl/latest
else
    # Linux O2 paths - use environment variables if set, otherwise use default paths
    # O2_ROOT = /home/alice/sw/ubuntu2204_x86-64/O2/latest
    # O2PHYSICS_ROOT = /home/alice/sw/ubuntu2204_x86-64/O2Physics/latest
    # FAIRLOGGER_ROOT = /home/alice/sw/ubuntu2204_x86-64/FairLogger/latest
    # FMT_ROOT = /home/alice/sw/ubuntu2204_x86-64/fmt/latest
    # FAIRMQ_ROOT = /home/alice/sw/ubuntu2204_x86-64/FairMQ/latest
    # BOOST_ROOT = /home/alice/sw/ubuntu2204_x86-64/boost/latest
    # ARROW_ROOT = /home/alice/sw/ubuntu2204_x86-64/arrow/latest
    # MS_GSL_ROOT = /home/alice/sw/ubuntu2204_x86-64/ms_gsl/latest
    UBUNTU_VERSION := $(shell lsb_release -rs | tr -d '.' | cut -c1-4)
    O2_ROOT ?= $(shell if [ -n "$$O2_ROOT" ]; then echo $$O2_ROOT; else echo ${HOME}/alice/sw/ubuntu$(UBUNTU_VERSION)_x86-64/O2/latest; fi)
    O2PHYSICS_ROOT ?= $(shell if [ -n "$$O2PHYSICS_ROOT" ]; then echo $$O2PHYSICS_ROOT; else echo ${HOME}/alice/sw/ubuntu$(UBUNTU_VERSION)_x86-64/O2Physics/latest; fi)
    FAIRLOGGER_ROOT ?= $(shell if [ -n "$$FAIRLOGGER_ROOT" ]; then echo $$FAIRLOGGER_ROOT; else echo ${HOME}/alice/sw/ubuntu$(UBUNTU_VERSION)_x86-64/FairLogger/latest; fi)
    FMT_ROOT ?= $(shell if [ -n "$$FMT_ROOT" ]; then echo $$FMT_ROOT; else echo ${HOME}/alice/sw/ubuntu$(UBUNTU_VERSION)_x86-64/fmt/latest; fi)
    FAIRMQ_ROOT ?= $(shell if [ -n "$$FAIRMQ_ROOT" ]; then echo $$FAIRMQ_ROOT; else echo ${HOME}/alice/sw/ubuntu$(UBUNTU_VERSION)_x86-64/FairMQ/latest; fi)
    BOOST_ROOT ?= $(shell if [ -n "$$BOOST_ROOT" ]; then echo $$BOOST_ROOT; else echo ${HOME}/alice/sw/ubuntu$(UBUNTU_VERSION)_x86-64/boost/latest; fi)
    ARROW_ROOT ?= $(shell if [ -n "$$ARROW_ROOT" ]; then echo $$ARROW_ROOT; else echo ${HOME}/alice/sw/ubuntu$(UBUNTU_VERSION)_x86-64/arrow/latest; fi)
    MS_GSL_ROOT ?= $(shell if [ -n "$$MS_GSL_ROOT" ]; then echo $$MS_GSL_ROOT; else echo ${HOME}/alice/sw/ubuntu$(UBUNTU_VERSION)_x86-64/ms_gsl/latest; fi)
endif

# Use alice sw environment when set (e.g. after 'alienv enter O2Physics/latest') so all
# preinstalled dependency includes and lib paths are picked up automatically.
ifneq ($(CPLUS_INCLUDE_PATH),)
INCS += $(addprefix -I,$(subst :, ,$(CPLUS_INCLUDE_PATH)))
endif
ifneq ($(LIBRARY_PATH),)
LIBS += $(addprefix -L,$(subst :, ,$(LIBRARY_PATH)))
endif

INCS += -I$(O2_ROOT)/include \
        -I$(O2_ROOT)/include/Framework \
        -I$(O2PHYSICS_ROOT)/include \
        -I$(FAIRLOGGER_ROOT)/include \
        -I$(FMT_ROOT)/include \
        -I$(FAIRMQ_ROOT)/include \
        -I$(BOOST_ROOT)/include \
        -I$(ARROW_ROOT)/include \
        -I$(MS_GSL_ROOT)/include
# PWGCF/DataModel is not installed by O2Physics (DataModel subdir commented out in build).
# Default to O2Physics source tree in alice dir; override with O2PHYSICS_SRC= if elsewhere.
O2PHYSICS_SRC ?= $(HOME)/alice/O2Physics
INCS += -I$(O2PHYSICS_SRC)

LIBS += -L$(O2_ROOT)/lib \
        -L$(O2PHYSICS_ROOT)/lib \
        -L$(FAIRLOGGER_ROOT)/lib \
        -L$(FMT_ROOT)/lib \
        -L$(ROOTSYS)/lib \
        -lO2Framework \
        -lO2FrameworkFoundation \
        -lO2PhysicsPWGCFCore \
        -lO2PhysicsJCorran \
        -lO2PhysicsAnalysisCore \
        -lFairLogger \
        -lfmt

CXXFLAGS       += $(INCS)

HDRSDICT = ../src/AliJBaseTrack.h ../src/AliJBaseEventHeader.h
HDRS += $(HDRSDICT) nanoDict.h

SRCS 		= $(HDRS:.h=.cxx)
OBJS		= $(HDRS:.h=.o)

%.o: %.cxx %.h
	#$(COMPILE.C) $(OUTPUT_OPTION) $(INCS) -c $<
	$(COMPILE.C) $(OUTPUT_OPTION) $<


$(PROGRAM):     $(OBJS) $(PROGRAM).C 
				$(CXX) -L$(PWD) $(PROGRAM).C $(CXXFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
				@echo "$(PROGRAM) done"

.PHONY : clean debug

clean : 
	@echo cleaning up
	rm -f $(OBJS) core *Dict*  $(PROGRAM).o $(PROGRAM) || echo "already cleaned"

cl:  clean $(PROGRAM)

debug:
	echo $(OBJS)

nanoDict.cc: $(HDRSDICT)
		@echo "Generating dictionary ..."
		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
		@rootcling nanoDict.cc -c -Done $(HDRSDICT)

# Platform-specific include paths
ifeq ($(UNAME_S),Darwin)
    INCLUDES = -I/Users/djkim/alice/sw/osx_arm64/O2/dev-local1/include \
              -I/Users/djkim/alice/sw/osx_arm64/ROOT/v6-32-06-alice1-local1/include \
              -I/Users/djkim/alice/sw/osx_arm64/fmt/10.1.1-local1/include
    LDFLAGS = -L/Users/djkim/alice/sw/osx_arm64/O2/dev-local1/lib \
              -L/Users/djkim/alice/sw/osx_arm64/ROOT/v6-32-06-alice1-local1/lib
else
    INCLUDES = -I$(O2_ROOT)/include \
              -I$(ROOTSYS)/include \
              -I$(FMT_ROOT)/include
    LDFLAGS = -L$(O2_ROOT)/lib \
              -L$(ROOTSYS)/lib
endif
