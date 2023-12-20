#include "TauAnalysis/ClassicSVfit/interface/MarkovChainRecorder.h"

#include <assert.h> // assert()

namespace classic_svFit
{

MarkovChainRecorder::MarkovChainRecorder(unsigned int numDimensions)
  : numDimensions_(numDimensions)
{
  assert(numDimensions > 0);
  x_.resize(numDimensions_);
}

MarkovChainRecorder::~MarkovChainRecorder()
{}

unsigned int
MarkovChainRecorder::getNumPoints()
{
  return points_.size();
}

const double*
MarkovChainRecorder::getPoint(unsigned int iPoint)
{
  assert(iPoint < points_.size());
  const std::vector<double>& point = points_[iPoint];
  for ( unsigned int iDimension = 0; iDimension < numDimensions_; ++iDimension )
  {
    x_[iDimension] = point[iDimension];
  }
  return x_.data();
}

double
MarkovChainRecorder::getValue(unsigned int iPoint)
{
  assert(iPoint < values_.size());
  return values_[iPoint];
}

double
MarkovChainRecorder::DoEval(const double* x) const
{
  std::vector<double> point(numDimensions_);
  for ( unsigned int iDimension = 0; iDimension < numDimensions_; ++iDimension )
  {
    point[iDimension] = x[iDimension];
  }
  points_.push_back(point);
  values_.push_back(x[numDimensions_]);
  return 0.;
}

}
