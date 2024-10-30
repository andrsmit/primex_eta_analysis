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

set target_list = {"full","empty"}
set  field_list = {"bfield","nobfield"}

foreach target ($target_list)
	foreach field ($field_list)
		set run_list = ${runListDirectory}/phase${phase}/he_${target}_${field}.txt
		echo "$run_list"
		if ( -f $run_list ) then
			echo " "
			echo "Getting flux for runs in: ${run_list}"
			foreach run (`cat ${run_list}`)
				
				set textFile = ${wd}/endpointEnergies/phase${phase}/${run}.txt
				
				# if text file already exists for this run, skip it:
				if ( -f $textFile ) then
					continue
				endif
				
				echo " "
				echo "${run}"
				ccdb dump /PHOTON_BEAM/endpoint_energy -r ${run} > $textFile
			end
		endif
	end
end
