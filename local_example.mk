CXX        := g++
WARN_FLAGS := -Wall -Wextra -Wshadow -Wno-unused-parameter -Wno-unused-function -Wno-unused-result -Wno-cast-function-type
SAN_FLAGS  := -fsanitize=address,undefined

CXXFLAGS   += -O0 -g3 $(WARN_FLAGS)
LDLIBS     += -llapack -lblas
