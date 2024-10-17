#!/usr/bin/sh
#
pass=2

outdir=/work/halld/home/andrsmit/primex_eta_analysis/eta_gg_matrix
script=${outdir}/run_rec_tree.csh
ij=0

submit_runs=go

#=========================#
# Choose data set:

dir_mss=/mss/halld/home/gxproj2/simulations/rp-2019-01-04092024-rdm-bkg-80M-batch01
if  [ ! -d $dir_mss ]; then
	echo "/mss directory does not exist."
	exit 1
fi
phase=phase1
dir_mss_length=${#dir_mss}
$dir_mss/hddm-primex_eta_he4-to-flat_coh-eta-to-gg-rnb-${run_number}.tar

#=========================#
#Job Resources:

workflow=eta_gg_matrix_batch1
account=halld
partition=production
ram=12GB
cores=5
run_time=2hr
disk=5GB
constraint=el9

#=========================#

# make sure config file exists:
cfg_file=${outdir}/pass${pass}/eta_ana.conf
if [ ! -f $cfg_file ]; then
	echo "JANA config file does not exist."
	exit 1
fi

dir_output=${outdir}/pass${pass}

# Loop over all files in directory and submit one job per file (per run):
for file in ${dir_mss}/hddm-primex_eta_he4-to-flat_coh-eta-to-gg-rnb-*.tar; do

	# Setup output directory where rootTrees will be written:
	dir_tree=${dir_output}/rootTrees/
	mkdir -p $dir_tree
	
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
	# this code relies on the fact that the filename starts with "hddm-primex_eta_he4-to-flat_coh-eta-to-gg-rnb-${run_number}..." 
	# immediately after $dir_mss/
	# So the run number is given by characters 48-53 after the end of $dir_mss/
	#
	cut_c1=$(($dir_mss_length+48))
	cut_c2=$(($dir_mss_length+52))
	run_number=`echo "$file" | cut -c ${cut_c1}-${cut_c2}`
	#echo "$run_number"
	
	# Setup output ROOT tree name:
	tree_file=${dir_tree}/${run_number}.root
	if [ -f $tree_file ]; then
		continue
	fi
	
	# Setup output ROOT file name:
	root_file=${dir_root}/${run_number}.root
	if [ -f $root_file ]; then
		continue
	#else
		#echo "Saving output to: ${root_file}"
	fi
	
	# Setup job submission file:
	jsub_file=${dir_jsub}/${run_number}.jsub
	if [ -f $jsub_file ]; then
		continue
	fi
	
	# Setup log files:
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
	job_name="eta_gg_matrix_batch1_pass${pass}_${run_number}"
	echo "submitting job ${job_name}"
	
	command="swif2 add-job -workflow ${workflow}"
	command="$command -name ${job_name} -account ${account} -partition ${partition}"
	command="$command -cores ${cores} -ram ${ram} -time ${run_time} -disk ${disk}"
	command="$command -constraint ${constraint}"
	command="$command -input `basename $file` mss:${file}"
	command="$command -stdout $log_file -stderr $err_file"
	command="$command $script `basename $file` $tree_file $root_file $cfg_file $jsub_file"
	
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
