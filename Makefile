CPP=g++
CFLAGS="-Wall"
# c++ -g -Wall -o uniwig uniwig.cpp -lm -lz -O2

BUILD_DIR := ./bin
LIB_DIR := /home/ys4aj/research/hmm/libBigWig/lib

uniwig: $(BUILD_DIR)/uniwig

opt: $(BUILD_DIR)/opt

$(BUILD_DIR)/uniwig:
	$(CPP) $(CFLAGS) -L$(LIB_DIR) -lBigWig -Wl,-R$(LIB_DIR) src/uniwig.cpp -o $(BUILD_DIR)/uniwig -lm -lz -O2

$(BUILD_DIR)/opt:
	$(CPP) $(CFLAGS) src/opt.cpp -o bin/opt

clean:
	-rm $(BUILD_DIR)/uniwig
	# -rm $(BUILD_DIR)/opt

test:
	./bin/uniwig test/test.bed 1 5 0
