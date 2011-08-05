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

CFLAGS = -g -pg -O0
LDFLAGS = -g -pg

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
	makedepend -f deps -- $(CFLAGS) -- $(C_SRCS) -- $(CPP_SRCS)

include deps
