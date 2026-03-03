CXX      := clang++
CXXFLAGS := -pedantic-errors -Wall -Wextra `root-config --cflags `
LDFLAGS	 := `root-config --glibs `
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps
TARGET   := analysis
INCLUDE  := -Iinclude/
SRC      := $(wildcard src/*.cpp)

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES := $(OBJECTS:.o=.d)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -MMD -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

-include $(DEPENDENCIES)

.PHONY: all build clean debug dev release info all

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g -fsanitize=address,undefined,leak
debug: all

dev: CXXFLAGS += -O2 -DDEV
dev: all

experimental: CXXFLAGS += -O2 -DDEV -DEXPERIMENTAL
experimental: all

release: CXXFLAGS += -O3 -march=native -flto=auto -funroll-loops -fno-plt
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

info:
	@echo "[*] Applications: "
	@$(foreach item, ${APP_DIR}, echo "\t$(item)";)
	@echo "[*] Object dir: "
	@$(foreach item, ${OBJ_DIR}, echo "\t$(item)";)
	@echo "[*] Sources: "
	@$(foreach item, ${SRC}, echo "\t$(item)";)
	@echo "[*] Objects: "
	@$(foreach item, ${OBJECTS}, echo "\t$(item)";)
	@echo "[*] Dependencies: "
	@$(foreach item, ${DEPENDENCIES}, echo "\t$(item)";)
