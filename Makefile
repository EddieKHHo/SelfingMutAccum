main: main.o InputPar.o World.o Popul.o Indiv.o SharedFunctions.o
	g++ main.o InputPar.o World.o Popul.o Indiv.o SharedFunctions.o -lgsl -lgslcblas -o main

main.o: main.cpp InputPar.h World.h Popul.h Indiv.h SharedFunctions.h IncludeFiles.h GlobalConstants.h
	g++ -c main.cpp
	
InputPar.o: InputPar.cpp InputPar.h IncludeFiles.h GlobalConstants.h
	g++ -c InputPar.cpp
	
World.o: World.cpp World.h Popul.h Indiv.h SharedFunctions.h IncludeFiles.h GlobalConstants.h
	g++ -c World.cpp
	
Popul.o: Popul.cpp Popul.h Indiv.h SharedFunctions.h IncludeFiles.h GlobalConstants.h
	g++ -c Popul.cpp

Indiv.o: Indiv.cpp Indiv.h SharedFunctions.h IncludeFiles.h GlobalConstants.h
	g++ -c Indiv.cpp
	
SharedFunctions.o: SharedFunctions.cpp SharedFunctions.h IncludeFiles.h GlobalConstants.h
	g++ -c SharedFunctions.cpp

clean:
	rm -r -f *.o main ExampleParameterFile.txt
