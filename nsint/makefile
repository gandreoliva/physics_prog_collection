all:
	echo "Choose: make metric or make structure (output binaries are *.bin)"

metric:
	gfortran tov_metric.f90 -o tov_metric.bin

structure:
	gfortran -Wall tov_structure.f90 -o tov_structure.bin

clean:
	rm *.bin