#!/usr/bin/sh
#
pass=1

outdir=/work/halld/home/andrsmit/primex_eta_analysis/bggen_ana
script=${outdir}/run_rec.csh
ij=0

submit_runs=go

#=========================#
# Choose data set:

evtgen="bggen-upd"
bkgd=""
run_period="2019-01"
nucleus="Helium"
nucleon="proton"
version="3895"

# Set primex phase (1, 2, or 3) or SRC from run period string:
phase=phase1
if [[ $run_period == "2021-08" ]]; then
	phase=phase2
elif [[ $run_period == "2022-08" ]]; then
	phase=phase3
elif [[ $run_period == "2021-11" ]]; then
	phase=src
fi

#=========================#
#Job Resources:

workflow=bggen_upd_primex_eta_analysis_$phase
account=halld
partition=production
ram=12GB
cores=5
run_time=2hr
disk=18GB
constraint=el9

#=========================#

dir_mss=/mss/halld/gluex_simulations/REQUESTED_MC/${evtgen}${bkgd}-rp-${run_period}-nucleus-${nucleus}-nucleon-${nucleon}_${version}
if  [ ! -d $dir_mss ]; then
	echo "directory does not exist."
	exit 1
fi
dir_mss_length=${#dir_mss}

#=========================#

# make sure config file exists:
cfg_file=${outdir}/pass${pass}/eta_ana.conf
if [ ! -f $cfg_file ]; then
	echo "JANA config file does not exist."
	exit 1
fi

dir_output=${outdir}/pass${pass}/${phase}/${nucleus}-${nucleon}-v${version}${bkgd}

# Loop over all files in directory and submit one job per file (per run):
for file in $dir_mss/bggen_upd_*_decay_evtgen_geant4_smeared.hddm.tar; do
	
	# Setup output directory where rootFiles will be written:
	dir_root=${dir_output}/rootFiles
	mkdir -p $dir_root
	
	# Setup output directory where jsub will be written:
	dir_jsub=${dir_output}/jsub
	mkdir -p $dir_jsub
	
	# Setup output directory where log files will be written:
	dir_log=${dir_output}/log
	mkdir -p $dir_log
	
	#
	# Get run number from filename:
	# this code relies on the fact that the filename starts with "bggen_upd_${run_number}..." immediately after $dir_mss/
	# So the run number is given by characters 12-17 after the end of $dir_mss/
	#
	cut_c1=$(($dir_mss_length+12))
	cut_c2=$(($dir_mss_length+17))
	run_number=`echo "$file" | cut -c ${cut_c1}-${cut_c2}`
	#echo "$run_number"
	
	# setup output file name:
	root_file=${dir_root}/${run_number}.root
	if [ -f $root_file ]; then
		continue
	#else
		#echo "Saving output to: ${root_file}"
	fi
	
	# setup jsub file:
	jsub_file=${dir_jsub}/${run_number}.jsub
	if [ -f $jsub_file ]; then
		continue
	fi
	
	# setup log files:
	log_file=${dir_log}/${run_number}.out
	err_file=${dir_log}/${run_number}.err
	
	# Get the file size of the .tar file so we know how much disk space to request:
	sum=0
	while IFS= read -r line; do
	if [[ $line =~ size=([0-9]+) ]]; then
	    size=${BASH_REMATCH[1]}
	    # Add the number to the sum
	    sum=$((sum + size))
	fi
    done < "$file"
	sum_gb=$(($sum / 1000000000))
	#echo "$sum_gb"
	disk="$(($sum_gb+2))GB"
	
	# Set job name:
	job_name="${nucleus}_${nucleon}_v${version}_${run_number}_pass${pass}"
	echo "submitting job ${job_name}"
	
	command="swif2 add-job -workflow ${workflow}"
	command="$command -name ${job_name} -account ${account} -partition ${partition}"
	command="$command -cores ${cores} -ram ${ram} -time ${run_time} -disk ${disk}"
	command="$command -constraint ${constraint}"
	command="$command -input `basename $file` mss:${file}"
	command="$command -stdout $log_file -stderr $err_file"
	command="$command $script `basename $file` $root_file $cfg_file $jsub_file"
	
	if [ $submit_runs == "go" ]; then
		echo "$command" > $jsub_file
		$command
		ij=$(($ij+1))
	fi
done
echo "submitted $ij jobs"
if [ $submit_runs == "go" ]; then
	swif2 run $workflow
fi
