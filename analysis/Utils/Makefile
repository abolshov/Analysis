CXX      := clang++
CXXFLAGS := -pedantic-errors -Wall -Wextra -fPIC `root-config --cflags `
LDFLAGS	 := `root-config --glibs ` -shared
BUILD    := build
OBJ_DIR  := $(BUILD)/objects
LIB_DIR  := $(BUILD)/libs
TARGET   := libUtils.so
INCLUDE  := -Iinclude/
SRC      := $(wildcard src/*.cpp)

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

all: build $(LIB_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -MMD -o $@

$(LIB_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(shell pwd)/$(LIB_DIR)/$(TARGET) $^ $(LDFLAGS)

-include $(DEPENDENCIES)

.PHONY: all build clean debug release info

build:
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(OBJ_DIR)

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(LIB_DIR)/*

info:
	@echo "[*] Libraries: "
	@$(foreach item, ${LIB_DIR}, echo "\t$(item)";)
	@echo "[*] Object dir: "
	@$(foreach item, ${OBJ_DIR}, echo "\t$(item)";)
	@echo "[*] Sources: "
	@$(foreach item, ${SRC}, echo "\t$(item)";)
	@echo "[*] Objects: "
	@$(foreach item, ${OBJECTS}, echo "\t$(item)";)
	@echo "[*] Dependencies: "
	@$(foreach item, ${DEPENDENCIES}, echo "\t$(item)";)