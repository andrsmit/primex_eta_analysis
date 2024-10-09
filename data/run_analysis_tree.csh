#!/usr/bin/tcsh
#

# command line inputs:

set skim_type = $1
set runnumber = $2
set tree_file = $3
set root_file = $4
set  cfg_file = $5
set jsub_file = $6

# source GlueX setup scripts:

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.csh
gxenv /home/andrsmit/version.xml

setenv CCDB_CONNECTION mysql://ccdb_user@hallddb-farm.jlab.org/ccdb
setenv JANA_CALIB_URL mysql://ccdb_user@hallddb-farm.jlab.org/ccdb

setenv HALLD_MY /home/andrsmit/halld_my

echo "pwd"
pwd

echo "ls -a"
ls -a

# copy JANA config file:

echo "cp ${cfg_file} jana.config"
cp ${cfg_file} jana.config

# set up and run hd_root:

set bin = ${HALLD_RECON_HOME}/Linux_Alma9-x86_64-gcc11.4.1/bin/hd_root

echo "$bin --config=jana.config ${skim_type}_${runnumber}_*.evio"
$bin --config=jana.config ${skim_type}_${runnumber}_*.evio

# move output root tree to appropriate directory:
if ( -f eta_gg.root ) then
	echo "mv eta_gg.root ${tree_file}"
	mv eta_gg.root ${tree_file}
endif

# move output root file to appropriate directory:
if ( -f hd_root.root ) then
	echo "mv hd_root.root ${root_file}"
	mv hd_root.root ${root_file}
endif

# delete local copies of input files and jsub file:

echo "rm -f jana.config"
rm -f jana.config

echo "rm -f ${skim_type}_${runnumber}_*.evio"
rm -f ${skim_type}_${runnumber}_*.evio

echo "rm -f ${jsub_file}"
rm -f ${jsub_file}
