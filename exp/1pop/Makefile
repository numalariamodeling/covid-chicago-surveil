CFLAGS=--ansi --pedantic -O2 -std=c++11 
#CFLAGS=--ansi --pedantic -g 
INCLUDE= -I../../src/
LDFLAGS=  ../../src/*.o

all: epifire gsl ex8_event_driven_SEIRS

epifire: 
	$(MAKE) -C ../src/

gsl:
	$(MAKE) -C ../gsl_subset/

model1g: ../../src/Event_Driven_NUCOVID.h ../../src/Utility.h ../../src/Time_Series.h single_modwave_nosmdet.cpp
	g++ $(CFLAGS) single_modwave_nosmdet.cpp $(INCLUDE) $(LDFLAGS) -o model1g

test: test.cpp
	g++ $(CFLAGS) test.cpp $(INCLUDE) $(LDFLAGS) -o test

model1g_apv: ../../src/Event_Driven_NUCOVID.h ../../src/Utility.h ../../src/Time_Series.h single_modwave_nosmdet_apvary.cpp
	g++ $(CFLAGS) single_modwave_nosmdet_apvary.cpp $(INCLUDE) $(LDFLAGS) -o model1g_apv

model1g_wdet: ../../src/Event_Driven_NUCOVID.h ../../src/Utility.h ../../src/Time_Series.h single_modwave_wsmdet.cpp
	g++ $(CFLAGS) single_modwave_wsmdet.cpp $(INCLUDE) $(LDFLAGS) -o model1g_wdet

clean:
	rm -f test_network ex1_mass_action ex2_percolation ex3_chain_binomial ex4_dynamic_net ex5_diff_eq ex6_network_diff_eq ex7_gillespie_network_SEIRS
