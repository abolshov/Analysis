main: readout selection
	@make -C main

readout:
	@make -C readout

selection: readout
	@make -C selection

.PHONY: selection readout clean main

clean:
	@make -C selection clean
	@make -C readout clean
	@make -C main clean