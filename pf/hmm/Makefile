CC = gcc
CPP = g++

BIN = ./main
C_SRCS =
CPP_SRCS = main.cpp halton.cpp

C_OBJS = $(C_SRCS:.c=.o)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
OBJS = $(C_OBJS) $(CPP_OBJS)

INC = -I./
LIBS =
kdlib = libkdtree.a

CFLAGS = -fPIC -g -O0 -Wall -pthread
LDFLAGS = $(kdlib) -lm

all: $(BIN)

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

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif
