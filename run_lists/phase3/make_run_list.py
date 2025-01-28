#!/usr/bin/env python 

import rcdb

minRun = 110657
maxRun = 119999

db = rcdb.RCDBProvider('mysql://rcdb@hallddb.jlab.org/rcdb2')

productionRunList = db.select_runs("@is_primex_production and status==1", minRun, maxRun)
with open("he_all_bfield.txt", "w") as file:
    for run in productionRunList:
        file.write(str(run.number)+"\n")

fullTargetRunList = db.select_runs("@is_primex_production and status==1 and target_type=='FULL & Ready'", minRun, maxRun)
with open("he_full_bfield.txt", "w") as file:
    for run in fullTargetRunList:
        file.write(str(run.number)+"\n")

emptyTargetRunList = db.select_runs("@is_primex_production and status==1 and target_type=='EMPTY & Ready'", minRun, maxRun)
with open("he_empty_bfield.txt", "w") as file:
    for run in emptyTargetRunList:
        file.write(str(run.number)+"\n")

#
# Note: The first list "he_all_bfield.txt" should ideally be the sum of the next two lists. 
# However, there is one extra one contained in it that has target_status='Filling'. 
# We should consider changing the 'status' for this run to 3 to indicate that it is a special run in the future.
#
