CPP=g++
CFLAGS="-Wall"
# c++ -g -Wall -o uniwig uniwig.cpp -lm -lz -O2

BUILD_DIR := ./bin

uniwig: $(BUILD_DIR)/uniwig

opt: $(BUILD_DIR)/opt

$(BUILD_DIR)/uniwig: 
	$(CPP) $(CFLAGS) src/uniwig.cpp -o $(BUILD_DIR)/uniwig -lm -lz -O2

$(BUILD_DIR)/opt:
	$(CPP) $(CFLAGS) src/opt.cpp -o bin/opt

clean:
	-rm $(BUILD_DIR)/uniwig
	-rm $(BUILD_DIR)/opt

test:
	./bin/uniwig test/test.bed 1 5 0
