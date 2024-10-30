# primex_eta_analysis
Software used for the eta photoproduction cross section measurement and 2-photon decay width extraction for the PrimEx-eta experiment in Hall D at Jefferson Lab.

Directory Structure:

1. run_lists:\
	Defines the list of runs used in the analysis for all three phases.

2. photon_flux:\
	Parses the CCDB for the tagged photon flux for each run (split amongst all three phases)
	
3. plugins:\
	All hd_root plugins which were used in the analysis of the data and MC simulations.
	Of particular importance is the plugin "eta_gg_tree" which is used to produce ROOT trees
		based on very minimal selection criteria. These trees are later further processed
		by the c++ code in the 'analyze_trees' directory (see below) where the full event selection
		criteria was applied, and 2-D histograms of invariant mass vs. production angle were filled
		for the eta->2gamma decay channel.

4. analyze_trees:\
	c++ code to analyze ROOT trees and apply selection criteria for exclusive eta(->2gamma) events
	
