#!/usr/bin/tcsh
#
mkdir -p log

set vetoOption = 1
set config = eta_ana_phase1.config

#submit three sets per job:

./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-proton-v3895 > log/Helium-proton-v3895-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-proton-v3912 > log/Helium-proton-v3912-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-proton-v3913 > log/Helium-proton-v3913-v${vetoOption}.log &

./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-proton-v3931 > log/Helium-proton-v3931-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-proton-v3936 > log/Helium-proton-v3936-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-proton-v3937 > log/Helium-proton-v3937-v${vetoOption}.log &

./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-proton-v3938 > log/Helium-proton-v3938-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-proton-v3939 > log/Helium-proton-v3939-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-proton-v3940 > log/Helium-proton-v3940-v${vetoOption}.log &
