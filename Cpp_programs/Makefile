CC = g++
VER = c++0x
TLIB = pthread
CFLAGS =  -lm -std=$(VER) -l$(TLIB)
MPFRFLAGS = -lgmp -lmpfr
GMPFLAGS = -lgmpxx -lgmp 
INCLUDES = -I../Include/
DEPENDENCIES = ../Include/ ../Include/Geometric_types/ ../Include/Geometric_structures_and_algorithms/Well_separated_pair_decomposition_types/  ../Include/Geometric_functions/  ../Include/Measures/ ../Include/Polynomial_related_types/ ../Include/Measures/Measure_base/ ../Include/Exception_related_types/ ../Include/Numeric_traits_types/

experiment_BBV_vs_sampling_variable_n: experiment_BBV_vs_sampling_variable_n.o
	$(CC)  $(CFLAGS) $(INCLUDES) -o experiment_BBV_vs_sampling_variable_n experiment_BBV_vs_sampling_variable_n.o; $(RM) experiment_MPD_vs_sampling_variable_n.o *~
experiment_BBV_vs_sampling_variable_n.o: experiment_BBV_vs_sampling_variable_n.cpp $(DEPENDENCIES)
	$(CC)  $(CFLAGS) $(INCLUDES) -c experiment_BBV_vs_sampling_variable_n.cpp

experiment_MPD_vs_sampling_variable_n: experiment_MPD_vs_sampling_variable_n.o
	$(CC)  $(CFLAGS) $(INCLUDES) -o experiment_MPD_vs_sampling_variable_n experiment_MPD_vs_sampling_variable_n.o; $(RM) experiment_MPD_vs_sampling_variable_n.o *~
experiment_MPD_vs_sampling_variable_n.o: experiment_MPD_vs_sampling_variable_n.cpp $(DEPENDENCIES)
	$(CC)  $(CFLAGS) $(INCLUDES) -c experiment_MPD_vs_sampling_variable_n.cpp

experiment_MPD_vs_sampling_variable_s: experiment_MPD_vs_sampling_variable_s.o
	$(CC)  $(CFLAGS) $(INCLUDES) -o experiment_MPD_vs_sampling_variable_s experiment_MPD_vs_sampling_variable_s.o; $(RM) experiment_MPD_vs_sampling_variable_s.o *~
experiment_MPD_vs_sampling_variable_s.o: experiment_MPD_vs_sampling_variable_s.cpp $(DEPENDENCIES)
	$(CC)  $(CFLAGS) $(INCLUDES) -c experiment_MPD_vs_sampling_variable_s.cpp

experiment_MPD_vs_sampling_variable_d: experiment_MPD_vs_sampling_variable_d.o
	$(CC)  $(CFLAGS) $(INCLUDES) -o experiment_MPD_vs_sampling_variable_d experiment_MPD_vs_sampling_variable_d.o; $(RM) experiment_MPD_vs_sampling_variable_d.o *~
experiment_MPD_vs_sampling_variable_d.o: experiment_MPD_vs_sampling_variable_d.cpp $(DEPENDENCIES)
	$(CC)  $(CFLAGS) $(INCLUDES) -c experiment_MPD_vs_sampling_variable_d.cpp

experiment_MPD_vs_WSPD_variable_n: experiment_MPD_vs_WSPD_variable_n.o
	$(CC)  $(CFLAGS) $(INCLUDES) -o experiment_MPD_vs_WSPD_variable_n experiment_MPD_vs_WSPD_variable_n.o; $(RM) experiment_MPD_vs_WSPD_variable_n.o *~
experiment_MPD_vs_WSPD_variable_n.o: experiment_MPD_vs_WSPD_variable_n.cpp $(DEPENDENCIES)
	$(CC)  $(CFLAGS) $(INCLUDES) -c experiment_MPD_vs_WSPD_variable_n.cpp

experiment_MPD_vs_WSPD_variable_s: experiment_MPD_vs_WSPD_variable_s.o
	$(CC)  $(CFLAGS) $(INCLUDES) -o experiment_MPD_vs_WSPD_variable_s experiment_MPD_vs_WSPD_variable_s.o; $(RM) experiment_MPD_vs_WSPD_variable_s.o *~
experiment_MPD_vs_WSPD_variable_s.o: experiment_MPD_vs_WSPD_variable_s.cpp $(DEPENDENCIES)
	$(CC)  $(CFLAGS) $(INCLUDES) -c experiment_MPD_vs_WSPD_variable_s.cpp

experiment_MPD_vs_WSPD_variable_d: experiment_MPD_vs_WSPD_variable_d.o
	$(CC)  $(CFLAGS) $(INCLUDES) -o experiment_MPD_vs_WSPD_variable_d experiment_MPD_vs_WSPD_variable_d.o; $(RM) experiment_MPD_vs_WSPD_variable_d.o *~
experiment_MPD_vs_WSPD_variable_d.o: experiment_MPD_vs_WSPD_variable_d.cpp $(DEPENDENCIES)
	$(CC)  $(CFLAGS) $(INCLUDES) -c experiment_MPD_vs_WSPD_variable_d.cpp

experiment_MPD_vs_WSPD_variable_epsilon: experiment_MPD_vs_WSPD_variable_epsilon.o
	$(CC)  $(CFLAGS) $(INCLUDES) -o experiment_MPD_vs_WSPD_variable_epsilon experiment_MPD_vs_WSPD_variable_epsilon.o; $(RM) experiment_MPD_vs_WSPD_variable_epsilon.o *~
experiment_MPD_vs_WSPD_variable_epsilon.o: experiment_MPD_vs_WSPD_variable_epsilon.cpp $(DEPENDENCIES)
	$(CC)  $(CFLAGS) $(INCLUDES) -c experiment_MPD_vs_WSPD_variable_epsilon.cpp

clean: 
	$(RM) *.o *~
