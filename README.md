# ClassicSVfit
Latest version of SVfit_standalone algorithm

# Installation instructions
The SVfitPerformanceStudies package has been developed and tested with CMSSW 7_6_x.
It depends on the following other packages:
	TauAnalysis/SVfitTF

In order to install the code, execute:
	cmsrel CMSSW_7_6_3
	cd CMSSW_7_6_3/src
	cmsenv
	git clone https://github.com/veelken/ClassicSVfit TauAnalysis/ClassicSVfit
	git clone https://github.com/veelken/SVfitTF TauAnalysis/SVfitTF
	cd $CMSSW_BASE/src
	scram b -j 4

In case of compilation problems, please contact christian.veelken AT cern.ch .

