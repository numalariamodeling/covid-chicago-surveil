CFLAGS=--ansi --pedantic -O2 -std=c++11 
#CFLAGS=--ansi --pedantic -g 
INCLUDE= -I../../src/
LDFLAGS=  ../../src/*.o

all: epifire gsl ex8_event_driven_SEIRS

epifire: 
	$(MAKE) -C ../src/

gsl:
	$(MAKE) -C ../gsl_subset/

model2g: ../../src/Event_Driven_NUCOVID.h ../../src/Utility.h ../../src/Time_Series.h agegrp_modwave_nosmdet.cpp
	g++ $(CFLAGS) agegrp_modwave_nosmdet.cpp $(INCLUDE) $(LDFLAGS) -o model2g

