working_dir = $(shell pwd)


.PHONY: all
all: bin/mKS bin/view bin/test

.PHONY: clean
clean:
	rm -f bin/mKS bin/view bin/test


bin/mKS: source/*
	@mkdir -p bin
	g++ -Ofast -march=native -mtune=native -std=c++17 -Isource -obin/mKS source/mKS.cpp -llapack -lblas

bin/view: scripts/*
	@mkdir -p bin
	ln -sf $(working_dir)/scripts/view bin/view

bin/test: test/* source/*
	@mkdir -p bin
	g++ -O0 -g3 -std=c++17 -Itest -Isource -obin/test test/test.cpp -llapack -lblas
