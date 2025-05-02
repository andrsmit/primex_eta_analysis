#!/usr/bin/tcsh
#
mkdir -p log

set vetoOption = $1

./ana_trees -p3 -a6 -ceta_ana_phase3.config -v${vetoOption} -t0 -b0 > log/p3_v${vetoOption}_t0_b0.log;\
./ana_trees -p3 -a6 -ceta_ana_phase3.config -v${vetoOption} -t1 -b0 > log/p3_v${vetoOption}_t1_b0.log;\
./ana_trees -p3 -a6 -ceta_ana_phase3.config -v${vetoOption} -t0 -b1 > log/p3_v${vetoOption}_t0_b1.log;\
./ana_trees -p3 -a6 -ceta_ana_phase3.config -v${vetoOption} -t1 -b1 > log/p3_v${vetoOption}_t1_b1.log &
