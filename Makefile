bin/test: test/* source/*
	@mkdir -p bin
	g++ -O0 -g3 -std=c++17 -Itest -Isource -obin/test test/test.cpp -llapack -lblas
