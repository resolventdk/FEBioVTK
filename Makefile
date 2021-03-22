CC = g++  # C compiler
CFLAGS = -fPIC -g # C flags
LDFLAGS = -shared -fopenmp  # linking flags
RM = rm -f   # rm command
TARGET = libvtk_gcc64.so  # target lib
SDK_INC = -I../../  # path to febio header files
SDK_LIB = -L../../build/lib/ -lfecore_gcc64  # path to febio libraries
SRC_DIR = ./src

SRCS := $(wildcard $(SRC_DIR)/*.cpp)  # source files
OBJS = $(SRCS:.cpp=.o)  # objectives

.PHONY: all
all: ${TARGET}

# compile the sources
# must be done this way for one at a time compilation
# cannot use, $OBJS: $SRCS
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(SDK_INC)

# then link together as library remember to include ext libraries
$(TARGET): $(OBJS)
	$(CC) $(INC) ${LDFLAGS} -o $@ $^ $(SDK_LIB)

.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)
