CPP=g++
CFLAGS="-Wall"
# c++ -g -Wall -o uniwig uniwig.cpp -lm -lz -O2


uniwig:
	$(CPP) $(CFLAGS) uniwig.cpp -o uniwig -lm -lz  -O2

clean:
	-rm uniwig

test:
	./uniwig test.bed 1 5 0

