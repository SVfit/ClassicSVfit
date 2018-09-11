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
git clone https://github.com/SVfit/ClassicSVfit TauAnalysis/ClassicSVfit -b release_2018Mar20
git clone https://github.com/SVfit/SVfitTF TauAnalysis/SVfitTF
cd $CMSSW_BASE/src
scram b -j 4
```

In case of compilation problems, please sutmitt an issue on
https://github.com/SVfit/ClassicSVfit/issues

# Without CMSSW software

It is possible to build the software without CMSSW framework if the following prerequisites are satisfied (oldest software version the instructions were tested with):
- ROOT (6.10/3 or newer)
- GCC (6.3 or newer)

In order to build the software, please execute the following lines in any directory with write access:
```bash
git clone https://github.com/SVfit/ClassicSVfit TauAnalysis/ClassicSVfit
make -f TauAnalysis/ClassicSVfit/Makefile -j4
```

The test executables will be placed to `$PWD/TauAnalysis/ClassicSVfit/exec`. In order to use them:
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib

# either enter the full path to the executable, e.g.
./TauAnalysis/ClassicSVfit/exec/testClassicSVfit

# or make the executable available globally
export PATH=$PATH:$PWD/TauAnalysis/ClassicSVfit/exec
testClassicSVfit # run anywhere
```
You can add the export statements to your `$HOME/.bashrc` to make their effect permanent.

# Running instructions

- [Presentation, slides 2+3](https://indico.cern.ch/event/684622/contributions/2807248/attachments/1575090/2487044/presentation_tmuller.pdf)
- [Example(s)](https://github.com/SVfit/ClassicSVfit/blob/master/bin/testClassicSVfit.cc)

# Reference

If you use this code, please cite:                                                                                                    
```
   L. Bianchini, B. Calpas, J. Conway, A. Fowlie, L. Marzola, L. Perrini, C. Veelken,                                                  
   "Reconstruction of the Higgs mass in events with Higgs bosons decaying into a pair of tau leptons using matrix element techniques", 
   Nucl. Instrum. Meth. A 862 (2017) 54                                                                                                
```
