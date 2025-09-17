CXX      := clang++
CXXFLAGS := -pedantic-errors -Wall -Wextra -fPIC `root-config --cflags `
LDFLAGS	 := `root-config --glibs ` -shared
BUILD    := build
OBJ_DIR  := $(BUILD)/objects
SO_DIR  := $(BUILD)/shared
TARGET   := libReadout.so
INCLUDE  := -Iinclude/
SRC      := $(wildcard src/*.cpp)
LIB_DIR	 := ../libs

ifneq ($(MODE),)
  ifeq ($(MODE), debug)
    CXXFLAGS += -DDEBUG -g -fsanitize=address,undefined,leak
  endif
  ifeq ($(MODE), release)
    CXXFLAGS += -O3 -march=native -flto=auto -funroll-loops -fno-plt
  endif
  ifeq ($(MODE), dev)
    CXXFLAGS += -O2
  endif
endif

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(OBJECTS:.o=.d)

all: build $(SO_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -MMD -o $@

$(SO_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(shell pwd)/$(SO_DIR)/$(TARGET) $^ $(LDFLAGS)

-include $(DEPENDENCIES)

.PHONY: all build clean debug release info

build:
	@mkdir -p $(SO_DIR)
	@mkdir -p $(OBJ_DIR)

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(SO_DIR)/*

info:
	@echo "[*] Libraries: "
	@$(foreach item, ${SO_DIR}, echo "\t$(item)";)
	@echo "[*] Object dir: "
	@$(foreach item, ${OBJ_DIR}, echo "\t$(item)";)
	@echo "[*] Sources: "
	@$(foreach item, ${SRC}, echo "\t$(item)";)
	@echo "[*] Objects: "
	@$(foreach item, ${OBJECTS}, echo "\t$(item)";)
	@echo "[*] Dependencies: "
	@$(foreach item, ${DEPENDENCIES}, echo "\t$(item)";)
