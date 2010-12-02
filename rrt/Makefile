CC = gcc
CPP = g++

BIN = ./rrt
C_SRCS = kdtree.c
CPP_SRCS = random.cpp rrt.cpp

C_OBJS = $(C_SRCS:.c=.o)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
OBJS = $(C_OBJS) $(CPP_OBJS)

C_DEPS=$(C_SRCS:.c=.d)
CPP_DEPS=$(CPP_SRCS:.cpp=.d)
DEPS=$(C_DEPS) $(CPP_DEPS)

LIBS =  
INCS = 

CFLAGS = -g -O3 -Wall -lpthread -DUSE_LIST_NODE_ALLOCATOR -lm

all: $(BIN)

$(BIN): $(DEPS) $(OBJS)
	@echo "  [LD]    $@"
	@$(CPP) $(CFLAGS) -o $(BIN) $(OBJS) $(LIBS)

%.d: %.c
	@echo "  [DEP]   $@"
	@$(CC) $(INCS) $(CFLAGS) -MM $< -MF $@
	@sed -i 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\)\n\1 : \2,g' $*.d

%.d: %.cpp
	@echo "  [DEP]   $@"
	@$(CPP) $(INCS) $(CFLAGS) -MM $< -MF $@
	@sed -i 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\)\n\1 : \2,g' $*.d

%.o: %.c
	@echo "  [CC]    $@"
	@$(CC) $(CFLAGS) $(INCS) -c $< -o $@

%.o: %.cpp
	@echo "  [CPP]   $@"
	@$(CPP) $(CFLAGS) $(INCS) -c $< -o $@

clean:
	rm -rf $(BIN) $(BIN).exe $(C_OBJS) $(CPP_OBJS) $(DEPS)

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif
