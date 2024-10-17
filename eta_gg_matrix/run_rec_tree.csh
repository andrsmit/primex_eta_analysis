#!/usr/bin/tcsh
#

# command line inputs:

set input_file = $1
set  tree_file = $2
set  root_file = $3
set   cfg_file = $4
set  jsub_file = $5

# source GlueX setup scripts:

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.csh
gxenv /home/andrsmit/version.xml

setenv CCDB_CONNECTION mysql://ccdb_user@hallddb-farm.jlab.org/ccdb
setenv JANA_CALIB_URL mysql://ccdb_user@hallddb-farm.jlab.org/ccdb
setenv JANA_CALIB_CONTEXT variation=mc

setenv HALLD_MY /home/andrsmit/halld_my

echo "pwd"
pwd

echo "ls -a"
ls -a

# extract files from .tar file:
echo "tar xvf ${input_file}"
tar xvf ${input_file}
ls -lrth
ls -lrth */*/*/*/*/*/hddm/*.hddm

# copy JANA config file:

echo "cp ${cfg_file} jana.conf"
cp ${cfg_file} jana.conf

# set up hd_root:

set bin = ${HALLD_RECON_HOME}/Linux_Alma9-x86_64-gcc11.4.1/bin/hd_root

set local_fname = `basename ${root_file}`
echo "$bin --config=jana.conf -o ${local_fname} *.hddm"
$bin --config=jana.conf -o ${local_fname} *.hddm

# move output root tree to appropriate directory:
if ( -f eta_gg.root ) then
	echo "mv eta_gg.root ${tree_file}"
	mv eta_gg.root ${tree_file}
endif

# move output root file to appropriate directory:
if ( -f ${local_fname} ) then
	echo "mv ${local_fname} ${root_file}"
	mv ${local_fname} ${root_file}
endif

# delete local copies of input files and jsub file:

echo "rm -rf work"
rm -rf work

echo "rm -f jana.conf"
rm -r jana.conf

echo "rm -f ${input_file}"
rm -f ${input_file}

echo "rm -f ${jsub_file}"
rm -f ${jsub_file}
