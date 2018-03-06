# ClassicSVfit
Latest version of SVfit_standalone algorithm

# Installation instructions
The SVfitPerformanceStudies package has been tested with CMSSW 9_4_4.
It depends on the following other packages:
	TauAnalysis/SVfitTF

In order to install the code, execute:

```
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
git clone https://github.com/SVfit/ClassicSVfit TauAnalysis/ClassicSVfit
git clone https://github.com/SVfit/SVfitTF TauAnalysis/SVfitTF
cd $CMSSW_BASE/src
scram b -j 4
```

In case of compilation problems, please sutmitt an issue on
https://github.com/SVfit/ClassicSVfit/issues
