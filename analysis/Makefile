main: readout selection utils
	@make -C main

utils: 
	@make -C utils

readout:
	@make -C readout

selection: readout utils
	@make -C selection

.PHONY: selection utils readout clean main

clean:
	@make -C utils clean
	@make -C selection clean
	@make -C readout clean
	@make -C main clean