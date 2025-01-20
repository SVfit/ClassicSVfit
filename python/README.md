# FastMTT

Documentation for FastMTT implementation in python. It's performance is of the order of C++ version, but probably worse (one should expect around 3.5 times slower calculations).

Version is standalone and do not need any installations apart from standard libraries (numpy, pandas, os, scipy, matplotlib, pyplot, argparse).

To see example usage of the code, please use FastMTT_test.py file and execute it with:

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








