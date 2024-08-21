#!/usr/bin/tcsh
#

# command line inputs:

set outdir      = $1
set skim_type   = $2
set target_type = $3
set runnumber   = $4
set phase       = $5
set jfile       = $6

# source GlueX setup scripts:

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.csh
gxenv /home/andrsmit/version.xml

setenv CCDB_CONNECTION mysql://ccdb_user@hallddb-farm.jlab.org/ccdb
setenv JANA_CALIB_URL mysql://ccdb_user@hallddb-farm.jlab.org/ccdb

setenv HALLD_MY /home/andrsmit/halld_my

echo "ls -a"
ls -a

# copy JANA config file:

echo "cp /work/halld/home/andrsmit/primex_eta_analysis/data/config/eta_ana.conf jana.conf"
cp /work/halld/home/andrsmit/primex_eta_analysis/data/config/eta_ana.conf jana.conf

# remove symlink created by swif and replace with local file:

#if ( -f /cache$mss_fname ) then
#	echo "copying file from cache..."
#	rm $loc_fname
#	cp /cache$mss_fname $loc_fname
#endif

# set up and run hd_root:

set bin = ${HALLD_RECON_HOME}/Linux_Alma9-x86_64-gcc11.4.1/bin/hd_root

echo "$bin --config=jana.conf -o ${runnumber}.root ${skim}_${runnumber}_*.evio"
$bin --config=jana.conf -o ${runnumber}.root ${skim}_${runnumber}_*.evio

# move output root file to appropriate directory:

echo "mv ${runnumber}.root $outdir/rootFiles/phase${phase}/${target}/${runnumber}.root"
mv ${runnumber}.root $outdir/rootFiles/phase${phase}/${target}/${runnumber}.root

# delete local copies of input files and jsub file:

echo "rm -f ${skim}_${runnumber}_*.evio"
rm -f ${skim}_${runnumber}_*.evio

echo "rm -f $jfile"
rm -f $jfile
