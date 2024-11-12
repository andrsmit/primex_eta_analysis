#
CERN_ROOT = /apps/cernlib/x86_64_rhel6/2005

all: clean MyReadConfig.o EtaAna.o

MyReadConfig.o: MyReadConfig.cc Makefile $(wildcard ./*.h)
	@echo Compiling MyReadConfig.cc;
	g++ -c MyReadConfig.cc `root-config --cflags --glibs` -I${HALLD_RECON_HOME}/src/libraries/include
#
EtaAna.o: EtaAna.cc MyReadConfig.o Makefile $(wildcard ./*.h)
	@echo Compiling EtaAna.cc; 
	g++ -c EtaAna.cc `root-config --cflags --glibs` -I${HALLD_RECON_HOME}/src/libraries/include
#
clean:
	rm -f EtaAna.o
	rm -f MyReadConfig.o
