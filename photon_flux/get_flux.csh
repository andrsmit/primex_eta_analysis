#!/usr/bin/tcsh
#

set phase = $1

set outdir = `pwd`

set prog = ${outdir}/primex_flux.py
if ( ! -f $prog ) then
	echo "python script not found"
	exit
endif

set target_list = {"empty"}
set  field_list = {"bfield","nobfield"}

foreach target ($target_list)
	foreach field ($field_list)
		set run_list = /work/halld/home/andrsmit/run_lists/primex_phase${phase}/he_${target}_${field}.txt
		if ( -f $run_list ) then
			echo " "
			echo "Getting flux for runs in: ${run_list}"
			foreach run (`cat ${run_list}`)
				echo " "
				echo "${run}"
				python ${prog} -b ${run} -e ${run}
				mv ${run}_flux.root rootFiles/phase${phase}/
				mv ${run}_tagh_ps_acc_cor.txt txtFiles/phase${phase}/
				mv ${run}_tagm_ps_acc_cor.txt txtFiles/phase${phase}/
			end
		endif
	end
end
