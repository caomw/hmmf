CC = gcc
CPP = g++

BIN = ./temp
C_SRCS = kdtree.c
CPP_SRCS = temp.cpp halton.cpp

C_OBJS = $(C_SRCS:.c=.o)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
OBJS = $(C_OBJS) $(CPP_OBJS)

INC = -I./
LIBS =

CFLAGS = -fPIC -g -pg -O3 -Wall -pthread
LDFLAGS =

all: $(BIN)
	rm -rf *.dat
	touch rrg.dat rrgp.dat miniout.dat traj.dat

$(BIN): $(OBJS)
	@echo "  [LD]    $@"
	@$(CPP) $(CFLAGS) -o $(BIN) $(OBJS) $(LIBS)

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

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif
