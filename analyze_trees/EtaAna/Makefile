#
CXX = g++
CXXFLAGS = -Iinclude `root-config --cflags --glibs`

SRC  := $(wildcard src/*.cc)
OBJS := $(patsubst src/%.cc, objects/%.o, $(SRC))
INCL := $(wildcard include/*.h)

all: $(OBJS)
	@echo build of EtaAna class complete

objects/%.o: src/%.cc ${INCL} | objects
	$(CXX) $(CXXFLAGS) -c $< -o $@
objects:
	mkdir -p objects
clean:
	rm -r objects/*.o
