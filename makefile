CC = gcc

COMPILER_FLAGS = -Wall -std=c99 -I/opt/local/include
LINKER_FLAGS = -lm -L/opt/local/lib

DEPS = SyndromeDecoding.c Matrix.c Queue.c

BUILD_FILES = SyndromeDecoding

.PHONY: build 
build:
	make SyndromeDecoding

SyndromeDecoding: $(DEPS)
	$(CC) $(DEPS) -o SyndromeDecoding $(COMPILER_FLAGS) $(LINKER_FLAGS)

.PHONY: clean
clean:
	rm -f *.o $(BUILD_FILES)