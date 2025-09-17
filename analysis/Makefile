MODE ?= dev

all: main

debug:
	$(MAKE) MODE=debug all

release:
	$(MAKE) MODE=release all

dev:
	$(MAKE) MODE=dev all

main: readout selection utils
	@echo "Building main in $(MODE) mode..."
	@make -C main

utils:
	@echo "Building utils in $(MODE) mode..."
	@make -C utils

readout:
	@echo "Building readout in $(MODE) mode..."
	@make -C readout

selection: readout utils
	@echo "Building selection in $(MODE) mode..."
	@make -C selection

.PHONY: selection utils readout clean main debug dev release

clean:
	@make -C utils clean
	@make -C selection clean
	@make -C readout clean
	@make -C main clean