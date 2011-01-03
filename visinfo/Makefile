CC = gcc
CPP = g++

BIN = sift
C_SRCS =
CPP_SRCS = sift.cpp

C_OBJS = $(C_SRCS:.c=.o)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
OBJS = $(C_OBJS) $(CPP_OBJS)

C_DEPS = $(C_SRCS:.c=.d)
CPP_DEPS = $(CPP_SRCS:.cpp=.d)
DEPS = $(C_DEPS) $(CPP_DEPS)

VLFEAT_LIB_DIR = lib/
VLFEAT_INC_DIR = .

LIBS = `pkg-config --libs opencv` -L$(VLFEAT_LIB_DIR) -lvl -lm
INCS = `pkg-config --cflags opencv` -I$(VLFEAT_INC_DIR)

CFLAGS = -O2 -g

all: $(BIN)

$(BIN): $(DEPS) $(OBJS)
	@echo "  [LD]    $@"
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

