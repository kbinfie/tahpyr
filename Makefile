.PHONY: all test clean

all: clean test

test:
	nextflow run main.nf 

clean:
	rm -rf .nextflow* results* work*