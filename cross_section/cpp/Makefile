#
CXX = g++
CXXFLAGS = -Iinclude `root-config --cflags --glibs`

SRC  := $(wildcard src/*.cc)
OBJS := $(patsubst src/%.cc, objects/%.o, $(SRC))
HEADERS := $(wildcard include/*.h)

EXEC = CrossSection

$(EXEC): $(OBJS)
	@echo Building Executable; $(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC)

objects/%.o: src/%.cc ${HEADERS} | objects
	$(CC) $(CXXFLAGS) -c $< -o $@
objects:
	mkdir -p objects
clean:
	rm -f objects/*.o $(EXEC)
