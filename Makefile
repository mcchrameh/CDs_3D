CC= mpic++   -fopenmp
CFLAGS=-c -Wall  -std=c++11 -O2  -g -Wextra  #-Werror  ### -Wno-write-strings 
LDFLAGS= -lm -pg
SOURCES=testCahnH3D.cpp CahnH3D.cpp #vtk_export.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=testCahn

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@  $(LDFLAGS)

.cpp.o:
	$(CC)  $(CFLAGS) $< -o $@  $(LDFLAGS)


clean:	 
	rm -f $(OBJECTS)
	rm -f $(EXECUTABLE)

run:
	mpirun -np 2 ./$(EXECUTABLE) 
