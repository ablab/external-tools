PROGRAM   := opera

OBJS := $(wildcard *.o)
CPP := $(wildcard *.cpp)

install : $(PROGRAM)

$(PROGRAM) : $(OBJS)     
	g++ -c -O3 -g -w $(CPP) $<
#g++ -c -g -fno-inline -w $(CPP) $<
	g++ -o ../bin/$(PROGRAM) *.o
	rm *.o
	
clean:
	rm *.o
