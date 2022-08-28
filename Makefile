CC := g++
CFLAGS := #-Wall
EXEC := 

SRC_DIR := src
OBJ_DIR := bin
SRC_FILES := $(shell find . -name '*.cpp') #$(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst ./$(SRC_DIR)/%.cpp,./$(OBJ_DIR)/%.o,$(SRC_FILES))

STATICFLAGS := -static-libgcc -static-libstdc++ -static -s

ifeq ($(OS),Windows_NT)
	EXEC := .exe
endif

all:	clear compile

compile:	$(OBJ_FILES)
	$(CC) $(CFLAGS) -o main$(EXEC) $(OBJ_FILES) $(STATICFLAGS) -fopenmp

$(OBJ_DIR)/%.o:	$(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -I$(SRC_DIR) -c -o $@ $< -O2 -fopenmp

clear:
	rm -r -f bin/*