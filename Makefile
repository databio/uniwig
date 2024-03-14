CPP=g++
CFLAGS="-Wall"
# c++ -g -Wall -o uniwig uniwig.cpp -lm -lz -O2

BUILD_DIR := ./bin
LIB_DIR := ../libBigWig

uniwig: $(BUILD_DIR)/uniwig

opt: $(BUILD_DIR)/opt

$(BUILD_DIR)/uniwig: src/uniwig.cpp
	$(CPP) $(CFLAGS) src/uniwig.cpp -L$(LIB_DIR) -lBigWig -Wl,-R$(LIB_DIR)  -o $(BUILD_DIR)/uniwig -lm -lz -O2 -g

$(BUILD_DIR)/opt:
	$(CPP) $(CFLAGS) src/opt.cpp -o bin/opt

clean:
	-rm $(BUILD_DIR)/uniwig

rebuild:
	$(MAKE) clean
	$(MAKE) uniwig

tests:
	$(MAKE) clean
	$(MAKE) uniwig
	time ./bin/uniwig -m 5 -w 10 test/test4.bed test/hg38.chrom.sizes ./data/bw/test

testbig:
	$(MAKE) clean
	$(MAKE) uniwig
	./bin/uniwig -m 5 -w 10000 test/test_big_sorted.bed test/hg38.chrom.sizes ./data/bw/test_big

testleak:
	$(MAKE) clean
	$(MAKE) uniwig
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./bin/uniwig -s -m 5 -w 10000 test/test_big_sorted.bed test/hg38.chrom.sizes ./data/bw/test_big
