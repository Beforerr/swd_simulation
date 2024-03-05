# Define shell
SHELL=bash

# Directory of precompiled VPIC libraries:

PROJECTDIR=/glade/u/home/xinan/vpic-kokkos-dev/build/bin

h1d:
	${PROJECTDIR}/vpic low2high_1d.cxx	

h2d:
	${PROJECTDIR}/vpic low2high_2d.cxx	

h2d_large:
	${PROJECTDIR}/vpic low2high_2d_large.cxx	

h2d_slowmode:
	${PROJECTDIR}/vpic low2high_2d_slowmode.cxx	

h2d_alfven:
	${PROJECTDIR}/vpic low2high_2d_alfven.cxx	

h2d_mpause:
	${PROJECTDIR}/vpic low2high_2d_alfven_mpause.cxx	

h2d_wind:
	${PROJECTDIR}/vpic low2high_2d_alfven_wind.cxx	

cleanall: 	
#	cd translate; make clean 
	-rm data/*
	-rm -rf hydro
	-rm -rf particle
	-rm -rf restart*
	-rm *.Linux
	-rm *.bin
	-rm *.vpc
	-rm *.btr
	-rm -rf rundata
	-rm -rf window
	-rm -rf fields
	-rm -rf spectrum
	-rm energy_spectrum
	-rm *.cy
	-rm *.op
	-rm *.ppu
	-rm outfile
	-rm info
	-rm *.eps
	-rm *.mpg
	-rm core.*
	-rm *~









