CC = gcc
CPP = g++

BIN = testvl
C_SRCS = 
CPP_SRCS = testvl.cpp

C_OBJS = $(C_SRCS:.c=.o)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
OBJS = $(C_OBJS) $(CPP_OBJS)

C_DEPS = $(C_SRCS:.c=.d)
CPP_DEPS = $(CPP_SRCS:.cpp=.d)
DEPS = $(C_DEPS) $(CPP_DEPS)

OPENCV_LIB_DIR = /usr/lib
OPENCV_INC_DIR = /usr/include/opencv
PLAYER_LIB_DIR = /home/mit/apps/player/lib
PLAYER_INC_DIR = /home/mit/apps/player/include/player-3.0

VLFEAT_LIB_DIR = lib/
VLFEAT_INC_DIR = .

LIBS = -L$(OPENCV_LIB_DIR) -L$(VLFEAT_LIB_DIR) -lcxcore -lcv -lhighgui -lvl -lm
INCS = -I$(OPENCV_INC_DIR) -I$(VLFEAT_INC_DIR)

CFLAGS = -O2 -Wall -Wextra

all: $(BIN)

#@$(CPP) $(CFLAGS) -o $(BIN) $(OBJS) $(LIBS) `pkg-config --libs playerc++`
$(BIN): $(DEPS) $(OBJS)
	@echo "  [LD]	 $@"
	@$(CPP) $(CFLAGS) -o $(BIN) $(OBJS) $(LIBS)
	@echo
	@echo "Add \"$(OPENCV_LIB_DIR)\" to your LD_LIBRARY_PATH env variable (if not added already) and export it"
	@echo "Add \"$(PLAYER_LIB_DIR)\" to your LD_LIBRARY_PATH env variable (if not added already) and export it"
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

