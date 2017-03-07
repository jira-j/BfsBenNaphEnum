enumchem: bfsenum.cpp bfsenum.hpp
	g++ -O3 -I/usr/local -L/usr/local/lib bfsenum.cpp -lboost_program_options-mt -std=gnu++11  -o bfsenum
