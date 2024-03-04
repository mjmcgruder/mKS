.PHONY: all
all: bin/mKS bin/test

.PHONY: clean
clean:
	rm -f bin/mKS bin/test

bin/mKS: source/*
	@mkdir -p bin
	g++ -O0 -g3 -std=c++17 -Isource -obin/mKS source/mKS.cpp -llapack -lblas

bin/test: test/* source/*
	@mkdir -p bin
	g++ -O0 -g3 -std=c++17 -Itest -Isource -obin/test test/test.cpp -llapack -lblas
