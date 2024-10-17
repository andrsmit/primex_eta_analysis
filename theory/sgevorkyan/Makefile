#
OBJS = PrimExCS.o CoulombFormFactor.o StrongFormFactor.o
CERN_ROOT = /apps/cernlib/x86_64_rhel6/2005

all: clean calc_cs

calc_cs: calc_cs.o Makefile $(OBJS)
	g++ calc_cs.o -o calc_cs $(OBJS) -L/usr/lib64/root \
	`root-config --cflags --glibs` -lMathMore -I${HALLD_RECON_HOME}/src/libraries/include -I${HALLD_SIM_HOME}/src/libraries/UTILITIES
#
PrimExCS.o: PrimExCS.cc Makefile $(wildcard ./*.h )
	@echo Compiling PrimExCS.cc; g++ -c PrimExCS.cc \
	`root-config --cflags --glibs` -lMathMore -I${HALLD_RECON_HOME}/src/libraries/include -I${HALLD_SIM_HOME}/src/libraries/UTILITIES
#
CoulombFormFactor.o: CoulombFormFactor.cc Makefile $(wildcard ./*.h)
	@echo Compiling CoulombFormFactor.cc; g++ -c CoulombFormFactor.cc \
	`root-config --cflags --glibs` -lMathMore -I${HALLD_RECON_HOME}/src/libraries/include -I${HALLD_SIM_HOME}/src/libraries/UTILITIES
#
StrongFormFactor.o: StrongFormFactor.cc Makefile $(wildcard ./*.h)
	@echo Compiling StrongFormFactor.cc; g++ -c StrongFormFactor.cc \
	`root-config --cflags --glibs` -lMathMore -I${HALLD_RECON_HOME}/src/libraries/include -I${HALLD_SIM_HOME}/src/libraries/UTILITIES
#
calc_cs.o: calc_cs.cc Makefile $(wildcard ./*.h)
	@echo Compiling calc_cs.cc; g++ -c calc_cs.cc \
	`root-config --cflags --glibs` -lMathMore -I${HALLD_RECON_HOME}/src/libraries/include -I${HALLD_SIM_HOME}/src/libraries/UTILITIES
#
clean:
	rm -f calc_cs.o PrimExCS.o CoulombFormFactor.o StrongFormFactor.o
