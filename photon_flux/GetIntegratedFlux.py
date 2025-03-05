#!/usr/bin/env python

# Tool for creating flux histograms from CCDB (ver 0.2)
# Modified version of the flux originally used by GlueX

import os,sys
import rcdb
from ROOT import TFile,TGraph,TH1F,TF1,gRandom, FILE
from optparse import OptionParser
from array import array
from datetime import datetime
import pprint
import math
import MySQLdb

import ccdb
from ccdb import Directory, TypeTable, Assignment, ConstantSet

def LoadCCDB():
    sqlite_connect_str = "mysql://ccdb_user@hallddb.jlab.org/ccdb"
    provider = ccdb.AlchemyProvider()                           # this class has all CCDB manipulation functions
    provider.connect(sqlite_connect_str)                        # use usual connection string to connect to database
    provider.authentication.current_user_name = "andrsmit"      # to have a name in logs

    return provider

def loadCCDBContextList(runPeriod, restVer):
    dbhost = "hallddb.jlab.org"
    dbuser = 'datmon'
    dbpass = ''
    dbname = 'data_monitoring'
    
    conn=MySQLdb.connect(host=dbhost, user=dbuser, db=dbname)
    curs=conn.cursor()    

    cmd = "SELECT revision,ccdb_context FROM version_info WHERE run_period=%s AND data_type='recon' AND revision<=%s ORDER BY revision DESC"
    curs.execute(cmd, [runPeriod, restVer])
    rows=curs.fetchall()
    return rows

def PSAcceptance(x, par):
    
    min = par[1]
    max = par[2]
    
    if x[0] > 2*min and x[0] < min + max:
        return par[0]*(1-2*min/x[0])
    elif x[0] >= min + max:
        return par[0]*(2*max/x[0] - 1)
    
    return 0.

def main():
    
    VARIATION   =  "default"
    
    conv_length   =  750e-6
    Be_rl         =  35.28e-2
    conv_rl       =  conv_length/Be_rl;   # 2.1259 10-3   Used in the PS acceptance determination
    
    ps_scale = 1./((7/9.) * conv_rl)
    
    
    parser = OptionParser(usage = "GetIntegrateFlux_PrimEx.py --run-list fileName --min-energy minEnergy --max-energy maxEnergy")
    
    parser.add_option("-r","--run-list", dest="run_list",
                      help="Run list for integrating photon flux")
    
    parser.add_option("-b","--min-energy", dest="min_energy",
                      help="Minimum beam energy to integrate photon flux over")
    
    parser.add_option("-e","--max_energy", dest="max_energy",
                      help="Maximum beam energy to integrate photon flux over")
    
    (options, args) = parser.parse_args(sys.argv)
    
    run_list_file = options.run_list
    min_energy    = float(options.min_energy)
    max_energy    = float(options.max_energy)
    
    ccdb_conn = LoadCCDB()
    
    fPSAcceptance = TF1("PSAcceptance", PSAcceptance, 2.0, 12.0, 3);
    
    integratedFlux    = 0.0
    integratedFluxErr = 0.0
    
    run_list = open(run_list_file,'r').readlines()
    for run in run_list:
        
        run_number = int(run)
        
        #---------------------------------------------------------#
        # Get endpoint_energy and endpoint_calib:
        
        photon_endpoint_assignment       = ccdb_conn.get_assignment("/PHOTON_BEAM/endpoint_energy",          run_number, VARIATION)
        photon_endpoint_calib_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/hodoscope/endpoint_calib", run_number, VARIATION)
        
        photon_endpoint       = float(photon_endpoint_assignment.constant_set.data_table[0][0])
        photon_endpoint_calib = float(photon_endpoint_calib_assignment.constant_set.data_table[0][0])
        
        #---------------------------------------------------------#
        # TAGH and TAGM scaled_energy_range:
        
        tagh_scaled_energy_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/hodoscope/scaled_energy_range",  run_number, VARIATION)
        tagm_scaled_energy_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/microscope/scaled_energy_range", run_number, VARIATION)
        
        tagh_scaled_energy = tagh_scaled_energy_assignment.constant_set.data_table
        tagm_scaled_energy = tagm_scaled_energy_assignment.constant_set.data_table
        
        #---------------------------------------------------------#
        # Get TAGH and TAGM flux:
        
        tagh_tagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagh_tagged", run_number, VARIATION)
        tagm_tagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagm_tagged", run_number, VARIATION)
        
        tagh_tagged_flux = tagh_tagged_flux_assignment.constant_set.data_table
        tagm_tagged_flux = tagm_tagged_flux_assignment.constant_set.data_table
        
        #---------------------------------------------------------#
        # Get parameters for PS acceptance correction:
        
        PS_accept_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/PS_accept", run_number, VARIATION)
        
        PS_accept = PS_accept_assignment.constant_set.data_table
        
        fPSAcceptance.SetParameters(float(PS_accept[0][0]), float(PS_accept[0][1]), float(PS_accept[0][2]))
        
        #---------------------------------------------------------#
        
        print("  Run",run_number,": Endpoint energy =",photon_endpoint)
        
        locDeltaE = photon_endpoint - photon_endpoint_calib
        
        # Loop over TAGH counters:
        
        for tagh_counter in range(274):
            
            locEnergyLo = photon_endpoint_calib*float(tagh_scaled_energy[tagh_counter][1]) + locDeltaE
            locEnergyHi = photon_endpoint_calib*float(tagh_scaled_energy[tagh_counter][2]) + locDeltaE
            locEnergy   = 0.5 * (locEnergyLo + locEnergyHi)
            
            locPSAcc = fPSAcceptance(locEnergy)
            
            if locPSAcc <= 0:
                locPSAcc = 10
            
            # correct flux in database by PS acceptance:
            
            locFluxCorr    = float(tagh_tagged_flux[tagh_counter][1])*ps_scale/locPSAcc
            locFluxCorrErr = float(tagh_tagged_flux[tagh_counter][2])*ps_scale/locPSAcc
            
            # If the energy of this counter is within our specified range, include it's flux into the total
            
            if locEnergy >= min_energy and locEnergy < max_energy:
                integratedFlux    += locFluxCorr
                integratedFluxErr += (locFluxCorrErr**2)
        
        # Do the same for the TAGM counters:
        
        for tagm_counter in range(102):
            
            locEnergyLo = photon_endpoint_calib*float(tagm_scaled_energy[tagm_counter][1]) + locDeltaE
            locEnergyHi = photon_endpoint_calib*float(tagm_scaled_energy[tagm_counter][2]) + locDeltaE
            locEnergy   = 0.5 * (locEnergyLo + locEnergyHi)
            
            locPSAcc = fPSAcceptance(locEnergy)
            
            if locPSAcc <= 0:
                locPSAcc = 10
            
            # correct flux in database by PS acceptance:
            
            locFluxCorr    = float(tagm_tagged_flux[tagm_counter][1])*ps_scale/locPSAcc
            locFluxCorrErr = float(tagm_tagged_flux[tagm_counter][2])*ps_scale/locPSAcc
            
            # If the energy of this counter is within our specified range, include it's flux into the total
            
            if locEnergy >= min_energy and locEnergy < max_energy:
                integratedFlux    += locFluxCorr
                integratedFluxErr += (locFluxCorrErr**2)
            
    print("\n\n")
    print("Integrated Photon Flux (",min_energy,"GeV -",max_energy,"GeV ):",
        "{0:.5E}".format(integratedFlux),"+/-","{0:.5E}".format(math.sqrt(integratedFlux)))
    print("\n\n")

## main function  
if __name__ == "__main__":
    main()
