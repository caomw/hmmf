CC = gcc
CPP = g++

BIN = ./main
C_SRCS = utils/kdtree.c
CPP_SRCS = utils/common.cpp systems/ship.cpp filter.cpp main.cpp

C_OBJS = $(C_SRCS:.c=.o)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
OBJS = $(C_OBJS) $(CPP_OBJS)

INCS = -I./ -I systems/ -I utils/ -I /home/pratik/apps/eigen/
LIBS = -L./

CFLAGS = -Wall -g -pg -O3
LDFLAGS = -g -pg

all: $(BIN)
	touch data/monte_carlo.dat data/rrg.dat data/rrgp.dat data/traj.dat

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
	rm -rf $(BIN) $(C_OBJS) $(CPP_OBJS) gmon.out

ctags:
	ctags -R $(C_SRCS) $(CPP_SRCS) *.h


depend:
	makedepend -f deps -- $(CFLAGS) -- $(C_SRCS) -- $(CPP_SRCS)

include deps
