# Makefile for the DMRG project
#-------------------------------------------------------------
CC = /opt/local/bin/g++

CPPFLAGS =#-DWITH_LAPACK
OPTFLAGS = -O3 -Wall
DBGFLAGS = -g -O0 -DDEBUG_MODE
MKLPATH = /opt/intel/mkl/lib
MKLINCLUDE = /opt/intel/mkl/include/intel64/lp64
#CCINCLUDE = $(PREFIX)/include
INCLUDE  = -I/usr/local/include -I$(MKLINCLUDE) 
CCFLAGS = -std=c++11 $(OPTFLAGS) $(INCLUDE) $(CPPFLAGS) 
CCGFLAGS = $(DBGFLAGS) $(INCLUDE) $(CPPFLAGS)  
LDFLAGS = $(OPTFLAGS) -L/usr/local/lib  -L$(MKLPATH) #-L$(LIBGALAHAD)   
LDGFLAGS = $(DBGFLAGS) -L/usr/local/lib -L$(MKLPATH) #-L$(LIBGALAHAD)   
LIBS = -lboost_system -lboost_filesystem -lmkl_intel_lp64 -lmkl_sequential \
	-lmkl_core -lmkl_blas95_lp64 -lmkl_lapack95_lp64 

#-------------------------------------------------------------
TAGT = a.out
GTAGT = dbg.out
SRCS = cmdargs.cc 
SRCS+= inputparams.cc 
SRCS+= taskparams.cc 
SRCS+= task.cc 
SRCS+= master_scheduler.cc
SRCS+= scheduler.cc
SRCS+= lattice.cc
SRCS+= latticelibrary.cc
SRCS+= main.cc

HDRS = optionparser.h cmdargs.h inputparams.h task.h scheduler.h \
       constants.h lattice.h \
       mytask.h

AUX = Makefile input.parm changelog
#-------------------------------------------------------------
# compilation and linking
VPATH = ./:./scheduler:./lattice

OBJS = $(patsubst %.cc,%.o, $(SRCS))
GOBJS= $(patsubst %.cc,debug_objs/%.o, $(SRCS))

all: $(TAGT)
$(TAGT) : $(OBJS)
	$(CC) -o $(TAGT) $(OBJS) $(LDFLAGS) $(LIBS)

$(GTAGT) : $(GOBJS)
	$(CC) -o $(GTAGT) $(GOBJS) $(LDFLAGS) $(LIBS)

%.o: %.cc
	$(CC) -c $(CCFLAGS) -o $@ $<

debug_objs/%.o: %.cc mkdebugdir
	$(CC) -c $(CCGFLAGS) -o $@ $<

mkdebugdir:
	mkdir -p debug_objs
	
#-------------------------------------------------------
# Detailed dependencies
DEPHDRS = optionparser.h cmdargs.h
cmdargs.o: $(DEPHDRS)
debug_objs/cmdargs.o: $(DEPHDRS)
DEPHDRS += inputparams.h
inputparams.o: inputparams.h
debug_objs/inputparams.o: $(DEPHDRS)
taskparams.o: inputparams.h
debug_objs/taskparams.o: $(DEPHDRS)
DEPHDRS += task.h
task.o: $(DEPHDRS)
debug_objs/task.o: $(DEPHDRS)
DEPHDRS += scheduler.h 
master_scheduler.o: $(DEPHDRS) 
scheduler.o: $(DEPHDRS) 
debug_objs/scheduler.o: $(DEPHDRS)
DEPHDRS += scheduler.h 
DEPHDRS += constants.h 
DEPHDRS += lattice.h 
lattice.o: lattice.h 
latticelibrary.o: lattice.h 
DEPHDRS += mytask.h 
main.o: $(DEPHDRS)

#-------------------------------------------------------
.PHONY: clean
clean:
	rm -f *.o debug_objs/*.o 
	
.PHONY: xclean
xclean:
	rm -f jobinfo* *.data*

