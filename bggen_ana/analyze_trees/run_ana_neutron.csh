#!/usr/bin/tcsh
#
mkdir -p log

set vetoOption = 1
set config = eta_ana_phase1.config

#submit three sets per job:

./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-neutron-v3894 > log/Helium-neutron-v3894-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-neutron-v3941 > log/Helium-neutron-v3941-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-neutron-v3942 > log/Helium-neutron-v3942-v${vetoOption}.log &

./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-neutron-v3943 > log/Helium-neutron-v3943-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-neutron-v3944 > log/Helium-neutron-v3944-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-neutron-v3945 > log/Helium-neutron-v3945-v${vetoOption}.log &

./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-neutron-v3946 > log/Helium-neutron-v3946-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-neutron-v3947 > log/Helium-neutron-v3947-v${vetoOption}.log;\
./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-neutron-v3948 > log/Helium-neutron-v3948-v${vetoOption}.log &

./ana_trees -p1 -a6 -v${vetoOption} -c${config} -xHelium-neutron-v3949 > log/Helium-neutron-v3949-v${vetoOption}.log &

