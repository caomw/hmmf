CC = gcc
CPP = g++

BIN = ./main
C_SRCS = kdtree.c
CPP_SRCS = common.cpp systems/singleint.cpp hmmf.cpp main.cpp

C_OBJS = $(C_SRCS:.c=.o)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
OBJS = $(C_OBJS) $(CPP_OBJS)

C_DEPS = $(C_SRCS:.c=.d)
CPP_DEPS = $(CPP_SRCS:.cpp=.d)
DEPS = $(C_DEPS) $(CPP_DEPS)

INC = -I./ -I systems/
LIBS = -L./

CFLAGS = -g -pg -O3
LDFLAGS =

all: $(BIN)
	rm -rf data/*.dat
	touch data/rrg.dat data/rrgp.dat data/traj.dat

$(BIN): $(OBJS)
	@echo "  [LD]    $@"
	@$(CPP) $(LDFLAGS) -o $(BIN) $(OBJS) $(LIBS)

%.o: %.c
	@echo "  [CPP]   $@"
	@$(CC) $(CFLAGS) $(INCS) -c $< -o $@

%.o: %.cpp
	@echo "  [CC]    $@"
	@$(CPP) $(CFLAGS) $(INCS) -c $< -o $@

clean:
	rm -rf $(BIN) $(C_OBJS) $(CPP_OBJS) $(DEPS) gmon.out

ctags:
	ctags -R $(C_SRCS) $(CPP_SRCS) *.h


depend:
	makedepend -- $(CFLAGS) -- $(C_SRCS) -- $(CPP_SRCS)


# DO NOT DELETE

kdtree.o: /usr/include/stdio.h /usr/include/_types.h
kdtree.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
kdtree.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
kdtree.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
kdtree.o: /usr/include/stdlib.h /usr/include/Availability.h
kdtree.o: /usr/include/AvailabilityInternal.h /usr/include/sys/wait.h
kdtree.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
kdtree.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
kdtree.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
kdtree.o: /usr/include/machine/_structs.h /usr/include/sys/resource.h
kdtree.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
kdtree.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
kdtree.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
kdtree.o: /usr/include/machine/types.h /usr/include/i386/types.h
kdtree.o: /usr/include/string.h /usr/include/secure/_string.h
kdtree.o: /usr/include/math.h /usr/include/architecture/i386/math.h kdtree.h
common.o: common.h /usr/include/math.h /usr/include/architecture/i386/math.h
common.o: /usr/include/sys/cdefs.h /usr/include/sys/time.h
common.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
common.o: /usr/include/i386/_types.h /usr/include/sys/_structs.h
common.o: /usr/include/machine/_structs.h /usr/include/i386/_structs.h
common.o: /usr/include/sys/appleapiopts.h /usr/include/time.h
common.o: /usr/include/_types.h /usr/include/_structs.h
common.o: /usr/include/sys/_select.h kdtree.h
systems/singleint.o: systems/singleint.h common.h /usr/include/math.h
systems/singleint.o: /usr/include/architecture/i386/math.h
systems/singleint.o: /usr/include/sys/cdefs.h /usr/include/sys/time.h
systems/singleint.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
systems/singleint.o: /usr/include/i386/_types.h /usr/include/sys/_structs.h
systems/singleint.o: /usr/include/machine/_structs.h
systems/singleint.o: /usr/include/i386/_structs.h
systems/singleint.o: /usr/include/sys/appleapiopts.h /usr/include/time.h
systems/singleint.o: /usr/include/_types.h /usr/include/_structs.h
systems/singleint.o: /usr/include/sys/_select.h kdtree.h
hmmf.o: hmmf.h common.h /usr/include/math.h
hmmf.o: /usr/include/architecture/i386/math.h /usr/include/sys/cdefs.h
hmmf.o: /usr/include/sys/time.h /usr/include/sys/_types.h
hmmf.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
hmmf.o: /usr/include/sys/_structs.h /usr/include/machine/_structs.h
hmmf.o: /usr/include/i386/_structs.h /usr/include/sys/appleapiopts.h
hmmf.o: /usr/include/time.h /usr/include/_types.h /usr/include/_structs.h
hmmf.o: /usr/include/sys/_select.h kdtree.h systems/singleint.h common.h
main.o: hmmf.h common.h /usr/include/math.h
main.o: /usr/include/architecture/i386/math.h /usr/include/sys/cdefs.h
main.o: /usr/include/sys/time.h /usr/include/sys/_types.h
main.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
main.o: /usr/include/sys/_structs.h /usr/include/machine/_structs.h
main.o: /usr/include/i386/_structs.h /usr/include/sys/appleapiopts.h
main.o: /usr/include/time.h /usr/include/_types.h /usr/include/_structs.h
main.o: /usr/include/sys/_select.h kdtree.h systems/singleint.h common.h
