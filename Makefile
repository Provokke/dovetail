FC = gfortran
FFLAGS ?= -O2 -std=f2018 -Wall -Wextra -fcheck=all
BUILD := build

.PHONY: test clean

test: $(BUILD)/check
	./$(BUILD)/check

$(BUILD)/check: src/dovetail.f90 test/check.f90
	mkdir -p $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -o $@ src/dovetail.f90 test/check.f90

clean:
	rm -rf $(BUILD)
