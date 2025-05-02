#
CLASSDIR = /work/halld/home/andrsmit/primex_eta_analysis/analyze_trees/EtaAna
OBJS     = $(wildcard ${CLASSDIR}/objects/*.o)

CXX = g++
CXXFLAGS = `root-config --cflags --glibs` -I${CLASSDIR}/include

all: ana_trees

ana_trees: analyze_trees.o Makefile $(OBJS)
	$(CXX) analyze_trees.o $(OBJS) -o ana_trees $(CXXFLAGS)
#
analyze_trees.o: analyze_trees.cc Makefile $(wildcard ./*.h)
	@echo Compiling analyze_trees.cc; $(CXX) -c analyze_trees.cc $(CXXFLAGS)
#
clean:
	rm -f ana_trees analyze_trees.o
