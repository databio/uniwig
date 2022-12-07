CPP=g++
CFLAGS="-Wall"
# c++ -g -Wall -o uniwig uniwig.cpp -lm -lz -O2

BUILD_DIR := ./bin
LIB_DIR := libBigWig/lib/lib

uniwig: $(BUILD_DIR)/uniwig

opt: $(BUILD_DIR)/opt

$(BUILD_DIR)/uniwig: src/uniwig.cpp
	$(CPP) $(CFLAGS) src/uniwig.cpp -L$(LIB_DIR) -lBigWig -Wl,-R$(LIB_DIR)  -o $(BUILD_DIR)/uniwig -lm -lz -O2

$(BUILD_DIR)/opt:
	$(CPP) $(CFLAGS) src/opt.cpp -o bin/opt

clean:
	-rm $(BUILD_DIR)/uniwig

test_old:  # remove?
	./bin/uniwig test/test.bed 1 5 0

tests:
	./bin/uniwig -m 5 test/test.bed test/hg38.chrom.sizes

rebuild:
	$(MAKE) clean
	$(MAKE) uniwig