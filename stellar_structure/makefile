MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# edit with path to shtawa
SHTAWA_DIR := ../../shtawa


lane_emden:
	gfortran -I${SHTAWA_DIR}/bin ${SHTAWA_DIR}/bin/shtawa.o $@.f90 -o $@.bin

clean:
	rm *.bin