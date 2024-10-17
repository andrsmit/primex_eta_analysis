#
CERN_ROOT = /apps/cernlib/x86_64_rhel6/2005

all: ana_trees

ana_trees: analyze_trees.o Makefile EtaAna.o
	g++ analyze_trees.o EtaAna.o -o ana_trees -L/usr/lib64/root \
	`root-config --cflags --glibs` -I${HALLD_RECON_HOME}/src/libraries/include
#
EtaAna.o: EtaAna.cc Makefile $(wildcard ./*.h)
	@echo Compiling EtaAna.cc; g++ -c EtaAna.cc \
	`root-config --cflags --glibs` -I${HALLD_RECON_HOME}/src/libraries/include
#
analyze_trees.o: analyze_trees.cc Makefile $(wildcard ./*.h)
	@echo Compiling analyze_trees.cc; g++ -c analyze_trees.cc \
	`root-config --cflags --glibs` -I${HALLD_RECON_HOME}/src/libraries/include
#
clean:
	rm -f ana_trees analyze_trees.o EtaAna.o
