#!/usr/bin/tcsh
#
set phase = $1

set runListDirectory = ""

#check if "PRIMEXDIR" exists:
if ($?PRIMEXDIR) then
	set runListDirectory = "${PRIMEXDIR}/run_lists"
else
	echo "Environment variable PRIMEXDIR is not set. Looking for run lists in this directory:"
	set runListDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/run_lists"
	echo "  ${runListDirectory}"
endif

set wd = `pwd`

set prog = ${wd}/primex_flux.py
if ( ! -f $prog ) then
	echo "python script not found"
	exit
endif

set target_list = {"full","empty"}
set  field_list = {"bfield","nobfield"}

foreach target ($target_list)
	foreach field ($field_list)
		set run_list = ${runListDirectory}/phase${phase}/he_${target}_${field}.txt
		if ( -f $run_list ) then
			
			set outdir_root = ${wd}/rootFiles/phase${phase}/${target}
			set outdir_text = ${wd}/txtFiles/phase${phase}/${target}
			
			# make sure directories exist to store rootFiles and txtFiles:
			mkdir -p ${outdir_root}
			mkdir -p ${outdir_text}
			
			# run python script to get flux for each run:
			echo " "
			echo "Getting flux for runs in: ${run_list}"
			foreach run (`cat ${run_list}`)
				
				# if rootFile already exists for this run, skip it:
				if ( -f ${outdir_root}/${run}_flux.root ) then
					continue
				endif
				
				echo " "
				echo "${run}"
				python ${prog} -b ${run} -e ${run}
				mv ${run}_flux.root ${outdir_root}/
				mv ${run}_tagh_ps_acc_cor.txt ${outdir_text}/
				mv ${run}_tagm_ps_acc_cor.txt ${outdir_text}/
			end
		endif
	end
end
