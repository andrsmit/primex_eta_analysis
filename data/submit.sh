#!/usr/bin/sh
#
pass=2

# set which run period you want to analyze:
phase=$1

#directory where the output will be stored:
outdir=/work/halld/home/andrsmit/primex_eta_analysis/data

# script that will execute the commands for each job:
script=${outdir}/run_analysis.csh

# set this to "go" to actually submit the jobs:
submit_runs=go

# directory where skim files are stored on tape:
mss_base_dir=/mss/halld/home/gxproj4/skim/phase${phase}

# this variable will just count how many jobs are submitted:
ij=0

#=========================#
#Job Resources:

workflow=primex_eta2g-skim_analysis
account=halld
partition=production
ram=12GB
cores=5
run_time=2hr
disk=12GB
constraint=el9

#=========================#

# look at the 2-photon skims only:
skim="eta2g-skim"

declare -a target_list=("full" "empty")
declare -a field_list=("bfield" "nobfield")

# directory where runs lists are stored:
run_list_dir=/work/halld/home/andrsmit/run_lists/primex_phase${phase}

#=========================#

# make sure config file exists:
cfg_file=${outdir}/pass${pass}/eta_ana.config
if [ ! -f $cfg_file ]; then
	echo "JANA config file does not exist."
	exit 1
fi

for target in "${target_list[@]}"
do
	for field in "${field_list[@]}"
	do
		echo "$target target, $field runs:"
		
		# get run list:
		run_list=${run_list_dir}/he_${target}_${field}.txt
		
		# skip this configuration if run list doesn't exist:
		if [ ! -f $run_list ]; then
			continue
		fi
		
		# set up output directory where hd_root files will be written:
		dir_rfile=${outdir}/pass${pass}/rootFiles/phase${phase}/${target}_target_${field}
		mkdir -p $dir_rfile
		
		# set up output directory where jsub will be written:
		dir_jsub=${outdir}/pass${pass}/jsub/phase${phase}/${target}_target_${field}
		mkdir -p $dir_jsub
		
		# loop over all runs in the list:
		while read run; do
			
			# Convert run number into 6 digit variable:
			
			run_number=$run
			if [[ $run_number -lt "100000" ]]; then
				run_number="0$run"
			fi
			
			# set up output hd_root file name:
			root_file=${dir_rfile}/${run_number}.root
			if [ -f $root_file ]; then
				continue
			fi
			
			# set up job submission file name:
			jsub_file=${dir_jsub}/${run_number}.jsub
			if [ -f $jsub_file ]; then
				continue
			fi
			
			dir_mss="${mss_base_dir}/ver01/phase${phase}-${field}-primex-eta-${target}/${skim}"
			
			#
			# For phases 1 and 2, all skim files were stored in "ver01".
			# However, for phase 3, most skim files are in 'ver00', while some runs had to be 
			# reprocessed and stored in 'ver01'.
			# To automate the procedure, I first check if the directory for a given run exists in 'ver01'.
			# If so, I use that version. If not, I check if it exists in 'ver00' and use that one instead.
			#
			
			if [ $phase -eq 3 ]; then
				dir_mss="${mss_base_dir}/ver01/phase${phase}-primex-eta-${target}/$skim"
				if [ ! -d ${dir_mss}/${run_number} ]; then
					dir_mss="${mss_base_dir}/ver00/phase${phase}-primex-eta-${target}/$skim"
					if [ ! -d ${dir_mss}/${run_number} ]; then
						echo "Skim files do not exist for Run${run_number}."
						continue
					fi
				fi
			fi
			
			echo "Submitting job for Run${run_number}"
			
			# set job name:
			job_name="${skim}_${run_number}_pass${pass}"
			
			command="swif2 add-job -workflow ${workflow}"
			command="$command -name ${job_name} -account ${account} -partition ${partition}"
			command="$command -cores ${cores} -ram ${ram} -time ${run_time} -disk ${disk}"
			command="$command -constraint ${constraint}"
			for file in ${dir_mss}/${run_number}/${skim}_${run_number}_*.evio; do
				if [ -f $file ]; then
					command="$command -input `basename $file` mss:${file}"
					#echo "  `basename $file`"
				fi
			done
			command="$command $script $skim $run_number $root_file $cfg_file $jsub_file"
			
			if [ $submit_runs == "go" ]; then
				echo "$command" > $jsub_file
				$command
			fi
			ij=$(($ij+1))
			
		done < $run_list # end loop over run numbers
	done # end loop over bfield configs
done # end loop over target configs

echo "submitted $ij jobs"
if [ $submit_runs == "go" ]; then
	swif2 run $workflow
fi
