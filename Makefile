CXXFLAGS = -std=c++17
LDFLAGS  =
LDLIBS   =

include local.mk

working_dir = $(shell pwd)


.PHONY: all
all: bin/mKS bin/view bin/test

.PHONY: clean
clean:
	rm -f bin/mKS bin/view bin/test


bin/mKS: source/*
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -Isource -o$@ source/mKS.cpp $(LDFLAGS) $(LDLIBS)

bin/view: scripts/*
	@mkdir -p bin
	ln -sf $(working_dir)/scripts/view bin/view

bin/test: source/* test/*
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -Isource -Itest -o$@ test/test.cpp $(LDFLAGS) $(LDLIBS)
