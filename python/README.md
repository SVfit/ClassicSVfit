# FastMTT

Author: Wiktor Matyszkiewicz
Last update: 27.02.2025

Documentation for FastMTT implementation in python. It's performance is of the order of C++ version, but probably worse (one should expect around 3.5 times slower calculations).

Presentation: https://indico.cern.ch/event/1467095/

Version is standalone and do not need any installations apart from standard libraries (numpy, pandas, os, scipy, matplotlib, pyplot, argparse).

To see example usage of the code, you can see Short_tutorial.py file with batch division. You can also try Long_tutorial.py, containing example of how to use FastMTT components. You can execute it with:

```
python3 FastMTT_test.py example_data.csv
```

# Basic usage

FastMTT has simple structure and is written in the FastMTT.py file only.
Moreover FastMTT is written as separate python class called FastMTT and to start using it, one have to make an instance of this class:

 ```
import FastMTT

fMTT = FastMTT.FastMTT()
 ```

The function responsible for communication with the program is .run and one can use the programm with:
```
fMTT.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET)
```

Then importing the results is made with .mass component:

```
masses = fMTT.mass
```

Input and output should be contained in numpy arrays with the structure (N, ...), where N is the number of events. For inputs one should need:

1) measuredTauLeptons -- (N, 2, 6), with the structure (number_of_events, taon for reconstruction, kinematic parameters)

Kinematic parameters should be:

lepton[0]: decay_type:
1 - TauToHad
2 - TauToElec
3 - TauToMu

lepton[1]: pt
lepton[2]: eta
lepton[3]: phi
lepton[4]: mass
lepton[5]: hadron decay mode (-1 for non-hadrons)

2) measuredMETx -- (N,) array with the information of x component of reconstructed MET.
3) measuredMETy -- (N,) array with the information of y component of reconstructed MET.
4) covMET -- (N, 2, 2) array with the covariance matrix, containting information of transfer function (Gaussian) between true and reconstructed MET. Event-by-event information will work at best.

For the output one will obtain one array of the size (N,), containing estimated invariant masses.

For another outpyts one can also get the reconstructed pT of the Higgs/Z and the momenta of reconstructed taons with:

'''
ptFast = fMTT.pt
p4_fast_1 = fMTT.tau1P4
p4_fast_2 = fMTT.tau2P4
'''

However one should be carefull, as the resolutions of these results are not perfect -- FastMTT was mainly invented for fast mass reconstruction.

# Batching system

We prepared additional function, that divide input data into batches. The simple function process_FastMTT is located in FastMTT_utils.py.

This approach is recommended, especially for the number of events > 10 thousands (memory allocation is very big in these cases without batching).

Additionally using it makes execution of the code faster -- probably due to the fact, that only the first batch is needed for highly-consuming memory allocation. Using batches allows for at least 2 times better time performance.

Inputs: measuredTauLeptons, METx, METy, covMET (as in previous section). Optional -- batch_size (5 000 by default) and log_interval (1 by default, so the code will raport upon its progress with each batch).

Outputs: mass, pT. Other can be added by hand in FastMTT_utlis.py file.

# Additional User Interface components

1) In case one want to see the likelihood of mass, one could plot it with the functions of FastMTT.
To set it, one should set the value of parameter .WhichLikelihoodPlot. It contains the information which event (of the number of N) should be plotted and saved into images/fastMTT directory.

Example usage:

```
fMTT.WhichLikelihoodPlot = 5
```

WhichLikelihoodPlot = -1 means, that no image will be plotted and is a default option.

2) In order to estimate event-by-event uncertainty, one can set:

```
fMTT.CalculateUncertainties = True
```

It calculates the uncertainty of the mass by estimating the contour, in which there should be masses with the probability in 1 sigma interval (according to the chi^2 test). Then the masses are calculated for the contour and highest and lowest masses give the interval for 1 sigma uncertainty. Additional arbitrary factor is used to adjust the results for chi^2 test.

The procedure produces long tails, but apart from that calculates uncertainties event by event quite ok ~ after some cuts results are aprox. Gaussian. It is also a bit time consuming -- doubles the time of calculation -- so it is disabled by default.

3) There are set two mass constraints, similar to each other (both disabled by default):

a) We modify the likelihood by the normal distribution, by setting:

'''
fMTT.myLikelihood.enable_mass_constraint = True
'''
We set the standard deviation to be equal to 10GeV, as this value seems to bring the best pT resolution. However one should avoid using it in the case of searching for heavy resonances.

b) Suggested by ICL team -- hard cut constraint on possible likelihood. To enable it one can set it and modify the range of window by:

'''
fMTT.myLikelihood.enable_window = True
fMTT.myLikelihood.window = [123, 127]
'''

Idea was already proved to improve the results in CP H->tau tau measurements, if used in proper way.

# Likelihood components

Likelihood.py contains is responsible for calculating likelihood for each point on the grid. These likelihood is crucial to the algorithm, as it is a basis for MLE method.

(Detailed description of physical work will be provided probably when the article will be written)

# Python wrapper

One can also use old C++ code with python wrapper.

This wrapper has been developed with the purpose of using it with the [ColumnFlow](https://columnflow.readthedocs.io/en/latest/) columnar Python-based framework. In order to use `ClassicSVFit` in Python, the [pybind11](https://pybind11.readthedocs.io/en/stable/basics.html) wrapper has been used. The wrapper for the different classes can be found here : [pybind_wrapper.cpp](https://github.com/oponcet/ClassicSVfit/blob/fastMTT_19_02_2019/wrapper/pybind_wrapper.cpp).

The cloned `ClassicSVFit` already contains the `pybind11` and the `wrapper` itself. A few more things to modify [hard-coded, please be patient !!!]
```bash
open TauAnalysis/ClassicSVfit/wrapper/CMakeList.txt
change L10 and L33
```
Now, this wrapper needs to be compiled with :
```bash
export LIBRARY_PATH=$LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/
cmake -S TauAnalysis/ClassicSVfit/wrapper/ -B TauAnalysis/ClassicSVfit/wrapper/
make -C TauAnalysis/ClassicSVfit/wrapper/
```
It should produce a `.so` file which can be used as a module in Python. For example you can import it like :
```py
from modules.extern.TauAnalysis.ClassicSVfit.wrapper.pybind_wrapper import *
```






