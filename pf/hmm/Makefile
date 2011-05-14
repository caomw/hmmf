CC = gcc
CPP = g++

BIN = ./main
C_SRCS = kdtree.c
CPP_SRCS = main.cpp halton.cpp

C_OBJS = $(C_SRCS:.c=.o)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
OBJS = $(C_OBJS) $(CPP_OBJS)

INC = -I./
LIBS = -L./

CFLAGS = -fPIC -g -pg -O0 -Wall -pthread
LDFLAGS = -lfigtree -lann_figtree_version

all: $(BIN)
	rm -rf *.dat
	touch rrg.dat rrgp.dat miniout.dat traj.dat

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
	ctags -R *.c *.cpp *.h

depend: Makefile $(CPP_SRCS)
	makedepend $(CPP_SRCS)


# DO NOT DELETE

main.o: common.h /usr/include/stdio.h /usr/include/_types.h
main.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
main.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
main.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
main.o: /usr/include/math.h /usr/include/architecture/i386/math.h
main.o: /usr/include/stdlib.h /usr/include/Availability.h
main.o: /usr/include/AvailabilityInternal.h /usr/include/sys/wait.h
main.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
main.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
main.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
main.o: /usr/include/machine/_structs.h /usr/include/sys/resource.h
main.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
main.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
main.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
main.o: /usr/include/machine/types.h /usr/include/i386/types.h
main.o: /usr/include/sys/time.h /usr/include/time.h /usr/include/_structs.h
main.o: /usr/include/sys/_select.h halton.h kdtree.h
halton.o: halton.h
