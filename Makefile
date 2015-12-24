CC=g++
CFLAGS=
LDFLAGS=`root-config --glibs --cflags`
SOURCES=readTree7414.cc
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=make

all: 
		$(CC) $(CFLAGS) $(SOURCES) $(LDFLAGS) -o $(EXECUTABLE)
			
clean:
		rm -rf *o $(EXECUTABLE)
