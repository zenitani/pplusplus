### C++ compiler
# CPP = g++
CPP = clang++

### files
HEADERS = vector3.h particle.h RK.h

###

all: ExB lorenz rossler poincare

ExB: sample_ExB.cpp $(HEADERS)
	$(CPP) -I. sample_ExB.cpp -o sample_ExB -lm
	./sample_ExB > data/ExB.dat

lorenz: sample_lorenz.cpp $(HEADERS)
	$(CPP) -I. sample_lorenz.cpp -o sample_lorenz -lm
	./sample_lorenz > data/lorenz.dat

rossler: sample_rossler.cpp $(HEADERS)
	$(CPP) -I. sample_rossler.cpp -o sample_rossler -lm
	./sample_rossler > data/rossler.dat

poincare: sample_poincare.cpp $(HEADERS)
	$(CPP) -I. sample_poincare.cpp -o sample_poincare -lm
	./sample_poincare > data/poincare.dat

clean: 
	rm sample_{ExB,lorenz,rossler,poincare}
	rm data/*.dat

# end
