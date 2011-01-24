CC = gcc
CPP = g++

BIN = main
C_SRCS = 
CPP_SRCS = main.cpp quickshift_cpu.cpp

C_OBJS = $(C_SRCS:.c=.o)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
OBJS = $(C_OBJS) $(CPP_OBJS)

C_DEPS = $(C_SRCS:.c=.d)
CPP_DEPS = $(CPP_SRCS:.cpp=.d)
DEPS = $(C_DEPS) $(CPP_DEPS)

LIBS = `pkg-config --libs opencv` -lm
INCS = `pkg-config --cflags opencv`

CFLAGS = -O2 -g

all: $(BIN)

$(BIN): $(DEPS) $(OBJS)
	@echo "  [LD]	 $@"
	@$(CPP) $(CFLAGS) -o $(BIN) $(OBJS) $(LIBS)
	@echo

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
	-rm -rf $(BIN) $(C_OBJS) $(CPP_OBJS) $(DEPS)

run:
	$(BIN)

