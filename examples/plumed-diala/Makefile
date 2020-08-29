default: run fes.dat

run:
	miniomm --input acemd_input

fes.dat:
	plumed sum_hills --hills HILLS

fes.png: fes.dat
	gnuplot fes.gp


clean:
	-rm bck.* restart.* HILLS colvar fes.dat plumed.log output.* miniomm_*

