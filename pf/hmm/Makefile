CC = gcc
CPP = g++

BIN = ./main
C_SRCS = kdtree.c
CPP_SRCS = single_integrator.cpp Logger.cpp

C_OBJS = $(C_SRCS:.c=.o)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
OBJS = $(C_OBJS) $(CPP_OBJS)

INC = -I./
LIBS = -L./

CFLAGS = -g -pg -O0
LDFLAGS = -lfigtree -lann_figtree_version

all: $(BIN)
	rm -rf *.dat
	touch rrg.dat rrgp.dat traj.dat

$(BIN): $(OBJS)
	@echo "  [LD]    $@"
	@$(CPP) $(LDFLAGS) -o $(BIN) $(OBJS) $(LIBS)

%.o: %.c
	@echo "  [CC]    $@"
	@$(CC) $(CFLAGS) $(INCS) -c $< -o $@

%.o: %.cpp
	@echo "  [CPP]   $@"
	@$(CPP) $(CFLAGS) $(INCS) -c $< -o $@

clean:
	rm -rf $(BIN) $(C_OBJS) $(CPP_OBJS) *.d gmon.out

ctags:
	ctags -R $(C_SRCS) $(CPP_SRCS) *.h

depend: Makefile $(CPP_SRCS)
	makedepend $(CPP_SRCS)


# DO NOT DELETE

single_integrator.o: common.h /usr/include/stdio.h /usr/include/_types.h
single_integrator.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
single_integrator.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
single_integrator.o: /usr/include/secure/_stdio.h
single_integrator.o: /usr/include/secure/_common.h /usr/include/math.h
single_integrator.o: /usr/include/architecture/i386/math.h
single_integrator.o: /usr/include/stdlib.h /usr/include/Availability.h
single_integrator.o: /usr/include/AvailabilityInternal.h
single_integrator.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
single_integrator.o: /usr/include/sys/appleapiopts.h
single_integrator.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
single_integrator.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
single_integrator.o: /usr/include/machine/_structs.h
single_integrator.o: /usr/include/sys/resource.h
single_integrator.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
single_integrator.o: /usr/include/sys/_endian.h
single_integrator.o: /usr/include/libkern/_OSByteOrder.h
single_integrator.o: /usr/include/libkern/i386/_OSByteOrder.h
single_integrator.o: /usr/include/alloca.h /usr/include/machine/types.h
single_integrator.o: /usr/include/i386/types.h /usr/include/sys/time.h
single_integrator.o: /usr/include/time.h /usr/include/_structs.h
single_integrator.o: /usr/include/sys/_select.h Logger.h kdtree.h figtree.h
Logger.o: Logger.h
