#
CERN_ROOT = /apps/cernlib/x86_64_rhel6/2005
CLASSDIR = ${PRIMEXDIR}/analyze_trees/EtaAna
OBJS = $(wildcard ${CLASSDIR}/*.o)

all: ana_trees

ana_trees: analyze_trees.o Makefile $(OBJS)
	g++ analyze_trees.o $(OBJS) -o ana_trees -L/usr/lib64/root \
	`root-config --cflags --glibs` -I${HALLD_RECON_HOME}/src/libraries/include -I${CLASSDIR}
#
analyze_trees.o: analyze_trees.cc Makefile $(wildcard ./*.h)
	@echo Compiling analyze_trees.cc; g++ -c analyze_trees.cc \
	`root-config --cflags --glibs` -I${HALLD_RECON_HOME}/src/libraries/include -I${CLASSDIR}
#
clean:
	rm -f ana_trees analyze_trees.o
