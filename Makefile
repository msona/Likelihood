
cflags:=-I. -I$(shell root-config --incdir)

iflags:=$(shell root-config --ldflags) $(shell root-config --libs) -lReflex -lGenVector -lCintex  -L/usr/lib/x86_64-linux-gnu/root5.34/  -lMinuit2 -lMinuit

Likelihood: Likelihood.cc 
	$(CXX) $+ $(cflags) $(iflags) -o $@

clean:
	rm Likelihood

