#!/usr/bin/tcsh
#

# set which run period you want to analyze:
set phase = $1

# directory where the output will be stored:
set outdir = /work/halld/home/andrsmit/primex_eta_analysis/data

# script that will execute the commands for each job:
set script = ${outdir}/run_analysis.csh

# set this to "go" to actually submit the jobs:
set submit_runs = no

# directory where skim files are stored on tape:
set mss_base_dir = /mss/halld/home/gxproj4/skim/phase${phase}

# this variable will just count how many jobs are submitted:
set ij = 0

#=========================#
#Job Resources:

set workflow   = primex_eta2g-skim_analysis
set account    = halld
set partition  = production
set ram        = 12GB
set cores      = 4
set time       = 120min
set disk       = 23GB
set constraint = el9

#=========================#

# look at the 2-photon skims only:
set skim = "eta2g-skim"

set target_list = {"full","empty"}
set  field_list = {"bfield","nobfield"}

# directory where runs lists are stored:
set run_list_dir = /work/halld/home/andrsmit/run_lists/primex_phase${phase}

foreach target ($target_list)
	foreach field ($field_list)
		
		echo "$target target, $field runs:"
		
		# get run list:
		set run_list = ${run_list_dir}/he_${target}_${field}.txt
		
		# skip this configuration if run list doesn't exist:
		if( ! -f $run_list ) then
			continue
		endif
		
		# loop over all runs in the list:
		foreach run ( `cat $run_list`)
			
			# convert 'run' into a six-digit number:
			set runnumber = $run
			if ($runnumber<"100000") then
				set runnumber = "0$run"
			endif
			
			# skip runs that have already been processed:
			if( -f $outdir/rootFiles/phase${phase}/${target}/${runnumber}.root) then
				continue
			endif
			
			# skip runs that are currently being processed:
			set jfile = ${outdir}/jsub/${skim}_${runnumber}.jsub
			if( -f $jfile) then
				continue
			endif
			
			set dir_mss = "${mss_base_dir}/ver01/phase${phase}-${field}-primex-eta-${target}/$skim"
			
			#
			# For phases 1 and 2, all skim files were stored in "ver01".
			# However, for phase 3, most skim files are in 'ver00', while some runs had to be 
			# reprocessed and stored in 'ver01'.
			# To automate the procedure, I first check if the directory for a given run exists in 'ver01'.
			# If so, I use that version. If not, I check if it exists in 'ver00' and use that one instead.
			#
			
			if( $phase == "3" ) then
				set dir_mss = "${mss_base_dir}/ver01/phase${phase}-primex-eta-${target}/$skim"
				if( ! -d ${dir_mss}/${runnumber} ) then
					set dir_mss = "${mss_base_dir}/ver00/phase${phase}-primex-eta-${target}/$skim"
					if( ! -d ${dir_mss}/${runnumber} ) then
						echo "Skim files do not exist for Run${runnumber}."
						continue
					endif
				endif
			endif
			
			set ij = `expr $ij + 1`
			echo "Submitting job for Run${runnumber}"
			
			set command = "swif2 add-job -workflow ${workflow}"
			set command = "$command -name ${skim}_${runnumber}"
			set command = "$command -account ${account} -partition ${partition}"
			set command = "$command -cores ${cores} -ram ${ram} -time ${time} -disk ${disk}"
			set command = "$command -constraint ${constraint}"
			foreach file (${dir_mss}/${runnumber}/${skim}_${runnumber}_*.evio)
				set command = "$command -input `basename $file` mss:${file}"
			end
			set command = "$command $script $outdir $skim $target $runnumber $phase $jfile"
			
			echo "$command" > $jfile
			
			if($submit_runs == "go") then
				$command
			endif
		end # end loop over runs
	end # end loop over bfield configs
end # end loop over target configs

echo "submitted $ij jobs"
