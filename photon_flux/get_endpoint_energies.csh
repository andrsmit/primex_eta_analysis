#!/usr/bin/tcsh
#

set phase = $1

set target_list = {"full","empty"}
set  field_list = {"bfield","nobfield"}

foreach target ($target_list)
	foreach field ($field_list)
		set run_list = /work/halld/home/andrsmit/run_lists/primex_phase${phase}/he_${target}_${field}.txt
		echo "$run_list"
		if ( -f $run_list ) then
			echo " "
			echo "Getting flux for runs in: ${run_list}"
			foreach run (`cat ${run_list}`)
				echo " "
				echo "${run}"
				ccdb dump /PHOTON_BEAM/endpoint_energy -r ${run} > endpointEnergies/phase${phase}/${run}.txt
			end
		endif
	end
end
