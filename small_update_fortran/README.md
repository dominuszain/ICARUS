ICARUS v1.21

A small update to the previous version. Now the output file also shows the nuclear mass number specified in the input file.

gfortran main.f90 -o icarus
./icarus < input

in the input file, you have to specify the address to the file that contains the neutron capture cross sections, and the number of lines in that file. here, you also specify the mass number of the target nucleus. 

note that the units in the neutron capture cross section file are mev for the energy, and mb for the cross sections.
in the output file, the units are kev for energy (kt), gk for temperature (t9), mb for macs, and cm3 s-1 mol-1 for the reaction rates. 
