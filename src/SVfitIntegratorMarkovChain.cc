#include "TauAnalysis/ClassicSVfit/interface/SVfitIntegratorMarkovChain.h"

#include "TauAnalysis/ClassicSVfit/interface/svFitAuxFunctions.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <assert.h>

enum { kUniform, kGaus, kNone };

namespace
{
  template <typename T>
  std::string format_vT(const std::vector<T>& vT)
  {
    std::ostringstream os;

    os << "{ ";

    unsigned numEntries = vT.size();
    for ( unsigned iEntry = 0; iEntry < numEntries; ++iEntry ) {
      os << vT[iEntry];
      if ( iEntry < (numEntries - 1) ) os << ", ";
    }

    os << " }";

    return os.str();
  }

  std::string format_vdouble(const std::vector<double>& vd)
  {
    return format_vT(vd);
  }
}

using namespace classic_svFit;

SVfitIntegratorMarkovChain::SVfitIntegratorMarkovChain(const std::string& initMode,
                   unsigned numIterBurnin, unsigned numIterSampling, unsigned numIterSimAnnealingPhase1, unsigned numIterSimAnnealingPhase2,
                   double T0, double alpha,
                   unsigned numChains, unsigned numBatches,
                   double epsilon0, double nu,
                   const std::string& treeFileName, int verbosity)
  : integrand_(0),
    x_(0),
    numIntegrationCalls_(0),
    numMovesTotal_accepted_(0),
    numMovesTotal_rejected_(0),
    errorFlag_(0),
    treeFileName_(treeFileName),
    treeFile_(0),
    tree_(0)
{
  if      ( initMode == "uniform" ) initMode_ = kUniform;
  else if ( initMode == "Gaus"    ) initMode_ = kGaus;
  else if ( initMode == "none"    ) initMode_ = kNone;
  else {
    std::cerr << "<SVfitIntegratorMarkovChain>:"
        << "Invalid Configuration Parameter 'initMode' = " << initMode << ","
        << " expected to be either \"uniform\", \"Gaus\" or \"none\" --> ABORTING !!\n";
    assert(0);
  }

//--- get parameters defining number of "stochastic moves" performed per integration
  numIterBurnin_ = numIterBurnin;
  numIterSampling_ = numIterSampling;

//--- get parameters defining maximum number of attempts to find a valid starting-position for the Markov Chain
  maxCallsStartingPos_ = 1000000;

//--- get parameters defining "simulated annealing" stage at beginning of integration
  numIterSimAnnealingPhase1_ = numIterSimAnnealingPhase1;
  numIterSimAnnealingPhase2_ = numIterSimAnnealingPhase2;
  numIterSimAnnealingPhase1plus2_ = numIterSimAnnealingPhase1_ + numIterSimAnnealingPhase2_;
  if ( numIterSimAnnealingPhase1plus2_ > numIterBurnin_ ) {
    std::cerr << "<SVfitIntegratorMarkovChain>:"
              << "Invalid Configuration Parameters 'numIterSimAnnealingPhase1' = " << numIterSimAnnealingPhase1_ << ","
              << " 'numIterSimAnnealingPhase2' = " << numIterSimAnnealingPhase2_ << ","
              << " sim. Annealing and Sampling stages must not overlap --> ABORTING !!\n";
    assert(0);
  }
  T0_ = T0;
  sqrtT0_ = TMath::Sqrt(T0_);
  alpha_ = alpha;
  if ( !(alpha_ > 0. && alpha_ < 1.) ) {
    std::cerr << "<SVfitIntegratorMarkovChain>:"
              << "Invalid Configuration Parameter 'alpha' = " << alpha_ << ","
              << " value within interval ]0..1[ expected --> ABORTING !!\n";
    assert(0);
  }
  alpha2_ = square(alpha_);

//--- get parameter specifying how many Markov Chains are run in parallel
  numChains_ = numChains;
  if ( numChains_ == 0 ) {
    std::cerr << "<SVfitIntegratorMarkovChain>:"
        << "Invalid Configuration Parameter 'numChains' = " << numChains_ << ","
        << " value greater 0 expected --> ABORTING !!\n";
    assert(0);
  }

  numBatches_ = numBatches;
  if ( numBatches_ == 0 ) {
    std::cerr << "<SVfitIntegratorMarkovChain>:"
              << "Invalid Configuration Parameter 'numBatches' = " << numBatches_ << ","
              << " value greater 0 expected --> ABORTING !!\n";
    assert(0);
  }
  if ( (numIterSampling_ % numBatches_) != 0 ) {
    std::cerr << "<SVfitIntegratorMarkovChain>:"
              << "Invalid Configuration Parameter 'numBatches' = " << numBatches_ << ","
              << " factor of numIterSampling = " << numIterSampling_ << " expected --> ABORTING !!\n";
    assert(0);
  }

  epsilon0_ = epsilon0;
  nu_ = nu;

  verbosity_ = verbosity;
}

SVfitIntegratorMarkovChain::~SVfitIntegratorMarkovChain()
{
  if ( verbosity_ >= 1 ) {
    std::cout << "<SVfitIntegratorMarkovChain::~SVfitIntegratorMarkovChain>:" << std::endl;
    std::cout << " integration calls = " << numIntegrationCalls_ << std::endl;
    std::cout << " moves: accepted = " << numMovesTotal_accepted_ << ", rejected = " << numMovesTotal_rejected_
              << " (fraction = " << (double)numMovesTotal_accepted_/(numMovesTotal_accepted_ + numMovesTotal_rejected_)*100. << "%)" << std::endl;
  }

  delete [] x_;
}

void SVfitIntegratorMarkovChain::setIntegrand(gPtr_C g, const double* xl, const double* xu, unsigned d)
{
  numDimensions_ = d;

  delete [] x_;
  x_ = new double[numDimensions_];

  xMin_.resize(numDimensions_);
  xMax_.resize(numDimensions_);
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    xMin_[iDimension] = xl[iDimension];
    xMax_[iDimension] = xu[iDimension];
  }

  epsilon0s_.resize(numDimensions_);
  epsilon_.resize(numDimensions_);
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    epsilon0s_[iDimension] = epsilon0_;
  }

  p_.resize(2*numDimensions_);   // first N entries = "significant" components, last N entries = "dummy" components
  q_.resize(numDimensions_);     // "potential energy" E(q) depends in the first N "significant" components only
  prob_ = 0.;

  u_.resize(2*numDimensions_);   // first N entries = "significant" components, last N entries = "dummy" components
  pProposal_.resize(numDimensions_);
  qProposal_.resize(numDimensions_);

  probSum_.resize(numChains_*numBatches_);
  for ( vdouble::iterator probSum_i = probSum_.begin();
  probSum_i != probSum_.end(); ++probSum_i ) {
    (*probSum_i) = 0.;
  }
  integral_.resize(numChains_*numBatches_);

  integrand_ = g;
}

void SVfitIntegratorMarkovChain::registerCallBackFunction(const ROOT::Math::Functor& function)
{
  callBackFunctions_.push_back(&function);
}

void SVfitIntegratorMarkovChain::integrate(gPtr_C g, const double* xl, const double* xu, unsigned d, double& integral, double& integralErr)
{
  setIntegrand(g, xl, xu, d);

  if ( !integrand_ ) {
    std::cerr << "<SVfitIntegratorMarkovChain>:"
              << "No integrand function has been set yet --> ABORTING !!\n";
    assert(0);
  }

  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    xMin_[iDimension] = xl[iDimension];
    xMax_[iDimension] = xu[iDimension];
    if ( verbosity_ >= 1 ) {
      std::cout << "dimension #" << iDimension << ": min = " << xMin_[iDimension] << ", max = " << xMax_[iDimension] << std::endl;
    }
  }

//--- CV: set random number generator used to initialize starting-position
//        for each integration, in order to make integration results independent of processing history
rnd_.SetSeed(12345);

  numMoves_accepted_ = 0;
  numMoves_rejected_ = 0;

  unsigned k = numChains_*numBatches_;
  unsigned m = numIterSampling_/numBatches_;

  numChainsRun_ = 0;

  if ( treeFileName_ != "" ) {
    treeFile_ = new TFile(treeFileName_.data(), "RECREATE");
    tree_ = new TTree("tree", "Markov Chain transitions");
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      std::string branchName = Form("x%u", iDimension);
      tree_->Branch(branchName.data(), &x_[iDimension]);
    }
    tree_->Branch("move", &treeMove_);
    tree_->Branch("integrand", &treeIntegrand_);
  }

  for ( unsigned iChain = 0; iChain < numChains_; ++iChain ) {
    bool isValidStartPos = false;
    if ( initMode_ == kNone ) {
      prob_ = evalProb(q_);
      if ( prob_ > 0. ) {
      bool isWithinBounds = true;
      for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
        double q_i = q_[iDimension];
        if ( !(q_i > 0. && q_i < 1.) ) isWithinBounds = false;
      }
      if ( isWithinBounds ) {
        isValidStartPos = true;
      } else {
        if ( verbosity_ >= 1 ) {
          std::cerr << "<SVfitIntegratorMarkovChain>:"
                    << "Warning: Requested start-position = " << format_vdouble(q_) << " not within interval ]0..1[ --> searching for valid alternative !!\n";
        }
      }
      } else {
        if ( verbosity_ >= 1 ) {
          std::cerr << "<SVfitIntegratorMarkovChain>:"
                    << "Warning: Requested start-position = " << format_vdouble(q_) << " returned probability zero --> searching for valid alternative !!";
        }
      }
    }
    unsigned iTry = 0;
    while ( !isValidStartPos && iTry < maxCallsStartingPos_ ) {
      initializeStartPosition_and_Momentum();
      prob_ = evalProb(q_);
      if ( prob_ > 0. ) {
        isValidStartPos = true;
      } else {
        if ( iTry > 0 && (iTry % 100000) == 0 ) {
          if ( iTry == 100000 ) std::cout << "<SVfitIntegratorMarkovChain::integrate>:" << std::endl;
          std::cout << "try #" << iTry << ": did not find valid start-position yet." << std::endl;
        }
      }
      ++iTry;
    }
    if ( !isValidStartPos ) continue;

    for ( unsigned iMove = 0; iMove < numIterBurnin_; ++iMove ) {
//--- propose Markov Chain transition to new, randomly chosen, point
      bool isAccepted = false;
      bool isValid = true;
      do {
  makeStochasticMove(iMove, isAccepted, isValid);
      } while ( !isValid );
    }

    unsigned idxBatch = iChain*numBatches_;

    for ( unsigned iMove = 0; iMove < numIterSampling_; ++iMove ) {
//--- propose Markov Chain transition to new, randomly chosen, point;
//    evaluate "call-back" functions at this point
      bool isAccepted = false;
      bool isValid = true;
      do {
        makeStochasticMove(numIterBurnin_ + iMove, isAccepted, isValid);
      } while ( !isValid );
      if ( isAccepted ) {
        ++numMoves_accepted_;
      } else {
        ++numMoves_rejected_;
      }

      updateX(q_);
      for ( std::vector<const ROOT::Math::Functor*>::const_iterator callBackFunction = callBackFunctions_.begin();
            callBackFunction != callBackFunctions_.end(); ++callBackFunction ) {
        (**callBackFunction)(x_);
      }

      if ( tree_ ) {
        treeMove_ = iMove;
        treeIntegrand_ = prob_;
        tree_->Fill();
      }

      if ( iMove > 0 && (iMove % m) == 0 ) ++idxBatch;
      assert(idxBatch < (numChains_*numBatches_));
      probSum_[idxBatch] += prob_;
    }

    ++numChainsRun_;
  }

  for ( unsigned idxBatch = 0; idxBatch < probSum_.size(); ++idxBatch ) {
    integral_[idxBatch] = probSum_[idxBatch]/m;
    if ( verbosity_ >= 1 ) std::cout << "integral[" << idxBatch << "] = " << integral_[idxBatch] << std::endl;
  }

//--- compute integral value and uncertainty
//   (eqs. (6.39) and (6.40) in [1])
  integral = 0.;
  for ( unsigned i = 0; i < k; ++i ) {
    integral += integral_[i];
  }
  integral /= k;

  integralErr = 0.;
  for ( unsigned i = 0; i < k; ++i ) {
    integralErr += square(integral_[i] - integral);
  }
  if ( k >= 2 ) integralErr /= (k*(k - 1));
  integralErr = TMath::Sqrt(integralErr);

  if ( verbosity_ >= 1 ) std::cout << "--> returning integral = " << integral << " +/- " << integralErr << std::endl;

  errorFlag_ = ( numChainsRun_ >= 0.5*numChains_ ) ? 0 : 1;

  ++numIntegrationCalls_;
  numMovesTotal_accepted_ += numMoves_accepted_;
  numMovesTotal_rejected_ += numMoves_rejected_;

  if ( tree_ ) {
    tree_->Write();
  }
  delete treeFile_;
  treeFile_ = 0;
  //delete tree_;
  tree_ = 0;

  if ( verbosity_ >= 1 ) print(std::cout);
}

void SVfitIntegratorMarkovChain::print(std::ostream& stream) const
{
  stream << "<SVfitIntegratorMarkovChain::print>:" << std::endl;
  for ( unsigned iChain = 0; iChain < numChains_; ++iChain ) {
    double integral = 0.;
    for ( unsigned iBatch = 0; iBatch < numBatches_; ++iBatch ) {
      double integral_i = integral_[iChain*numBatches_ + iBatch];
      //std::cout << "batch #" << iBatch << ": integral = " << integral_i << std::endl;
      integral += integral_i;
    }
    integral /= numBatches_;
    //std::cout << "<integral> = " << integral << std::endl;

    double integralErr = 0.;
    for ( unsigned iBatch = 0; iBatch < numBatches_; ++iBatch ) {
      double integral_i = integral_[iChain*numBatches_ + iBatch];
      integralErr += square(integral_i - integral);
    }
    if ( numBatches_ >= 2 ) integralErr /= (numBatches_*(numBatches_ - 1));
    integralErr = TMath::Sqrt(integralErr);

    std::cout << " chain #" << iChain << ": integral = " << integral << " +/- " << integralErr << std::endl;
  }
  std::cout << "moves: accepted = " << numMoves_accepted_ << ", rejected = " << numMoves_rejected_
            << " (fraction = " << (double)numMoves_accepted_/(numMoves_accepted_ + numMoves_rejected_)*100.
            << "%)" << std::endl;
}

//
//-------------------------------------------------------------------------------
//

void SVfitIntegratorMarkovChain::initializeStartPosition_and_Momentum()
{
//--- randomly choose start position of Markov Chain in N-dimensional space
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    bool isInitialized = false;
    while ( !isInitialized ) {
      double q0 = 0.;
      if ( initMode_ == kGaus ) q0 = rnd_.Gaus(0.5, 0.5);
      else q0 = rnd_.Uniform(0., 1.);
      if ( q0 > 0. && q0 < 1. ) {
  q_[iDimension] = q0;
  isInitialized = true;
      }
    }
  }
  if ( verbosity_ >= 2 ) {
    std::cout << "<SVfitIntegratorMarkovChain::initializeStartPosition_and_Momentum>:" << std::endl;
    std::cout << " q = " << format_vdouble(q_) << std::endl;
  }
}

void SVfitIntegratorMarkovChain::sampleSphericallyRandom()
{
//--- compute vector of unit length
//    pointing in random direction in N-dimensional space
//
//    NOTE: the algorithm implemented in this function
//          uses the fact that a N-dimensional Gaussian is spherically symmetric
//         (u is uniformly distributed over the surface of an N-dimensional hypersphere)
//
  double uMag2 = 0.;
  for ( unsigned iDimension = 0; iDimension < 2*numDimensions_; ++iDimension ) {
    double u_i = rnd_.Gaus(0., 1.);
    u_[iDimension] = u_i;
    uMag2 += (u_i*u_i);
  }
  double uMag = TMath::Sqrt(uMag2);
  for ( unsigned iDimension = 0; iDimension < 2*numDimensions_; ++iDimension ) {
    u_[iDimension] /= uMag;
  }
}

void SVfitIntegratorMarkovChain::makeStochasticMove(unsigned idxMove, bool& isAccepted, bool& isValid)
{
//--- perform "stochastic" move
//    (eq. 24 in [2])

//--- perform random updates of momentum components
  if ( idxMove < numIterSimAnnealingPhase1_ ) {
    for ( unsigned iDimension = 0; iDimension < 2*numDimensions_; ++iDimension ) {
      p_[iDimension] = sqrtT0_*rnd_.Gaus(0., 1.);
    }
  } else if ( idxMove < numIterSimAnnealingPhase1plus2_ ) {
    double pMag2 = 0.;
    for ( unsigned iDimension = 0; iDimension < 2*numDimensions_; ++iDimension ) {
      double p_i = p_[iDimension];
      pMag2 += p_i*p_i;
    }
    double pMag = TMath::Sqrt(pMag2);
    sampleSphericallyRandom();
    for ( unsigned iDimension = 0; iDimension < 2*numDimensions_; ++iDimension ) {
      p_[iDimension] = alpha_*pMag*u_[iDimension] + (1. - alpha2_)*rnd_.Gaus(0., 1.);
    }
  } else {
    for ( unsigned iDimension = 0; iDimension < 2*numDimensions_; ++iDimension ) {
      p_[iDimension] = rnd_.Gaus(0., 1.);
    }
  }

//--- choose random step size
  double exp_nu_times_C = 0.;
  do {
    double C = rnd_.BreitWigner(0., 1.);
    exp_nu_times_C = TMath::Exp(nu_*C);
  } while ( TMath::IsNaN(exp_nu_times_C) || !TMath::Finite(exp_nu_times_C) || exp_nu_times_C > 1.e+6 );
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    epsilon_[iDimension] = epsilon0s_[iDimension]*exp_nu_times_C;
  }

  // Metropolis algorithm: move according to eq. (27) in [2]

//--- update position components
//    by single step of chosen size in direction of the momentum components
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    qProposal_[iDimension] = q_[iDimension] + epsilon_[iDimension]*p_[iDimension];
  }

//--- ensure that proposed new point is within integration region
//   (take integration region to be "cyclic")
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    double q_i = qProposal_[iDimension];
    q_i = q_i - TMath::Floor(q_i);
    assert(q_i >= 0. && q_i <= 1.);
    qProposal_[iDimension] = q_i;
  }

//--- check if proposed move of Markov Chain to new position is accepted or not:
//    compute change in phase-space volume for "dummy" momentum components
//   (eqs. 25 in [2])
  double probProposal = evalProb(qProposal_);

  double deltaE = 0.;
  if      ( probProposal > 0. && prob_ > 0. ) deltaE = -TMath::Log(probProposal/prob_);
  else if ( probProposal > 0.               ) deltaE = -std::numeric_limits<double>::max();
  else if (                      prob_ > 0. ) deltaE = +std::numeric_limits<double>::max();
  else assert(0);

  // Metropolis algorithm: move according to eq. (13) in [2]
  double pAccept = TMath::Exp(-deltaE);

  double u = rnd_.Uniform(0., 1.);

  if ( u < pAccept ) {
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      q_[iDimension] = qProposal_[iDimension];
    }
    prob_ = probProposal;
    isAccepted = true;
  } else {
    isAccepted = false;
  }
}

void SVfitIntegratorMarkovChain::updateX(const std::vector<double>& q)
{
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    const double & q_i = q[iDimension];
    x_[iDimension] = (1. - q_i)*xMin_[iDimension] + q_i*xMax_[iDimension];
  }
}

double SVfitIntegratorMarkovChain::evalProb(const std::vector<double>& q)
{
  double prob = (*integrand_)(q.data(), numDimensions_, 0);
  return prob;
}
